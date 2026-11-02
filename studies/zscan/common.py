"""Shared code of the z-scan analysis: datasets, event selection, in-time hits, fits.

In-time hits (as the test-fibre analysis, analysis/intime_analysis.py): CFD time (relative to the MCP) within +-W_IN of the channel's
arrival peak and peak amplitude above A_HIT; random coincidences (dark counts: the CFD integral is anchored on the largest pulse of the
205 ns record) estimated from sidebands and subtracted.
  light per hit   = (sum of in-time integrals - N_rand <E>_sideband) / (N_in - N_rand)
  light per event = (sum of in-time integrals - N_rand <E>_sideband) / N_events   (events without a hit count as zero)
Integral: long (50 samples) for scintillating, short (15 samples) for the prompt Cherenkov pulses (the short window holds only ~40% of
a large scintillation pulse, checks/pulse_containment.py)."""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from datasets import DATASETS, DATA, REPO, MCP_SELECT, xpos, drs_dump
os.chdir(REPO); sys.path[:0] = [REPO, os.path.join(REPO, "scripts/compare")]
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
from configs.selection_config import get_particle_selection, get_service_drs_cut
from channels.maps.services import get_mcp_channels
from channels.maps.paths import load_channel_map

TYPES = ("Sci", "CerQuartz", "CerPlastic")
TLAB = {"Sci": "scintillating", "CerQuartz": "Cherenkov quartz", "CerPlastic": "Cherenkov plastic"}
TS_NS = 0.2
T_SEARCH = (380.0, 720.0)                         # in-time peak looked for here [TS, CFD time - MCP + 500]
W_IN = {"Sci": 20.0}                              # in-time window: +-20 TS (4 ns) scintillating, +-6 TS (1.2 ns) Cherenkov
W_IN_CER = 6.0
SIDEBANDS = ((-420.0, -120.0), (180.0, 480.0))    # random-coincidence sidebands, relative to the peak [TS]
A_HIT = {False: 40.0, True: 25.0}                 # hit threshold [ADC peak], 3 / 6 mm SiPM (single p.e. ~130 / ~45 ADC high)
MIN_SIGNAL = 100.0                                # light per hit [ADC] needed in every run for a per-channel fit
AMP_LEVELS = (1000.0, 1500.0, 2000.0)             # saturation check: fraction of events above these peak amplitudes
SAT_FRAC = 0.01                                   # channels with more than 1% of events above 1500 ADC (any run) are dropped
MAX_RESID = 0.06                                  # per-channel fits averaged only with run-to-run residuals below 6%
HF_ALL = 0.5                                      # all-event means (no MCP needed): channels hit in more than half of the events
MIN_FERS = 0.05                                   # FERS channels: mean above this (calibrated units) in every run
SPE_WIN = {False: (80.0, 200.0), True: (34.0, 70.0)}  # single-p.e. amplitude peak looked for here [ADC], 3 / 6 mm SiPM


class Dataset:
    """one z-scan dataset: runs with dumps, positions, groups, run sets, DRS channels: the calorimeter channel map (MCP channels and test
    fibres excluded) or the dataset's own channel map (test fibres: label -> channel; scintillating if the label says so, else Cherenkov)"""
    def __init__(self, name):
        c = DATASETS[name]; self.name, self.cfg = name, c
        self.runs = [r for r in c["runs"] if os.path.exists(drs_dump(name, r)) and xpos(r) <= c.get("xmax", 1.0)]
        self.x = {r: xpos(r) + c.get("x_offset", lambda r: 0.0)(r) for r in self.runs}
        self.group = {r: c["groups"](r) for r in self.runs}
        self.runsets = {k: [r for r in v if r in self.runs] for k, v in c.get("runsets", {"main": self.runs}).items()}
        self.daq = lambda r: c.get("daq", {}).get(r, r)
        if c.get("channel_map"):
            self.meta = {ch: dict(label=l, type="Sci" if "scint" in l.lower() else "Cer", is6mm=False)
                         for l, ch in load_channel_map(c["channel_map"]).items()}
        else:
            order = P._build_channel_order(P.build_drs_boards(run_number=c["map_run"]))
            excl = set().union(*(get_mcp_channels(self.daq(r)).values() for r in self.runs))
            if c.get("test_fibres"): excl |= set(load_channel_map("channels/data/testingfibers.json").values())
            self.meta = {ch: dict(type=v, tx=float(tx), ty=float(ty), is6mm=bool(s6)) for ch, b, g, cc, v, tx, ty, s6 in order if ch not in excl}
        self.shift, self.yc = {}, 0.0            # position matching (FERS datasets): hodoscope-y shift per run, window centre

    def w_in(self, ch):
        return W_IN.get(self.meta[ch]["type"], W_IN_CER)


def chip(ch):
    """DRS chip of a channel: the name without the channel number (DRS_Brg0_Board3_Group3_Channel5 -> DRS_Brg0_Board3_Group3)"""
    return ch.rsplit("_", 1)[0]


# ---------------- event selection ----------------
def pid(d, run):
    """positrons: beam-PID flags of the run's period"""
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for det, fire in get_particle_selection("electron", run).items():
        sel &= d[f"is_{det}_fired" if fire else f"is_{det}_vetoed"].astype(bool)
    return sel


def pid_hodo(d, run):
    """positrons (beam-PID flags of the run's period) and 1-3 hodoscope fibres per axis"""
    nx, ny = d["TTU_Hodo_nHitX"], d["TTU_Hodo_nHitY"]
    return pid(d, run) & (nx >= 1) & (nx <= 3) & (ny >= 1) & (ny <= 3)


def mcp_ok(d, run):
    """a pulse in MCP_SELECT above its configured cut (the dumped <MCP>_peak_value is the maximum in the configured window)"""
    return d[f"{MCP_SELECT}_peak_value"] > get_service_drs_cut(MCP_SELECT, run)[4]


def beam_window(d, sel, ds, run, mode):
    """beam within 1 cm of the median of the selected events; mode "matched": vertically within 1 cm of a common centre after
    shifting each run by the measured drift of the shower in the module (matching()), so all runs keep showers at the same place"""
    hx, hy = d["TTU_Hodo_X"], d["TTU_Hodo_Y"]
    w = np.abs(hx - np.median(hx[sel])) < 1.0
    if mode == "matched": return w & (np.abs(hy + ds.shift[run] - ds.yc) < 1.0)
    return w & (np.abs(hy - np.median(hy[sel])) < 1.0)


def load_drs(ds, run, mode):
    """per-event arrays of one run (channels of the dataset): E (integral per fibre type), A, T and the hodoscope for the events with an
    MCP pulse and a valid time reference; Eall for all selected positrons (no MCP: all-event means). None if too few.
    mode "spot" / "matched": hodoscope requirements and beam window (beam_window); "loose" (test-fibre check): positrons only, no
    hodoscope requirement (as the framework mcp_clean, which has no particle ID either)."""
    d = np.load(drs_dump(ds.name, run)); names = list(d["channels"]); daq = ds.daq(run)
    idx = [i for i, ch in enumerate(names) if ch in ds.meta]; chans = [names[i] for i in idx]
    pos = pid(d, daq) if mode == "loose" else pid_hodo(d, daq)
    s1 = pos & mcp_ok(d, daq) & (np.abs(d["T"]) < 3000).any(1)
    if s1.sum() < ds.cfg["min_events"]: return None
    if mode == "loose": sel, sel_all = s1, pos
    else: sel = s1 & beam_window(d, s1, ds, run, mode); sel_all = pos & beam_window(d, pos, ds, run, mode)
    short = np.array([ds.meta[ch]["type"] != "Sci" for ch in chans])[None, :]
    pick = lambda k, s: d[k][s][:, idx]
    return dict(chans=chans, n=int(sel.sum()), n_pos=int(sel_all.sum()),
                E=np.where(short, pick("Es", sel), pick("E", sel)), A=pick("A", sel), T=pick("T", sel),
                Eall=np.where(short, pick("Es", sel_all), pick("E", sel_all)),
                hx=d["TTU_Hodo_X"][sel], hy=d["TTU_Hodo_Y"][sel])


def intime(E, T, A, win, thr, sidebands=SIDEBANDS):
    """in-time hits of one channel in one run: arrival peak t0 [TS], hit mask, randoms expected in time and their mean integral"""
    valid = np.abs(T) < 3000
    w = T[valid & (T > T_SEARCH[0]) & (T < T_SEARCH[1])]
    if len(w) < 100: return None
    h, ed = np.histogram(w, bins=np.arange(T_SEARCH[0], T_SEARCH[1] + 1, 4.0)); k = int(np.argmax(h))
    t0 = float(np.median(w[(w >= ed[max(k - 3, 0)]) & (w < ed[min(k + 4, len(ed) - 1)])]))
    big = valid & (A > thr); inn = big & (np.abs(T - t0) < win)
    side = np.zeros(len(T), bool); span = 0.0
    for lo, hi in sidebands:
        side |= (T >= t0 + lo) & (T < t0 + hi); span += hi - lo
    side &= big
    n_rand = side.sum() * 2 * win / span
    return dict(t0=t0, inn=inn, n_rand=float(n_rand), e_side=float(E[side].mean()) if side.any() else 0.0)


def spe_of(a_out, e_out, is6):
    """single p.e. from out-of-time pulses (dark counts): amplitude peak in SPE_WIN, median integral of the pulses within 20% of it"""
    lo, hi = SPE_WIN[is6]; step = 2.0 if is6 else 4.0
    h, ed = np.histogram(a_out, bins=np.arange(lo, hi + step, step)); k = int(np.argmax(h))
    if h[k] < 20: return np.nan, np.nan
    a1 = ed[k] + step / 2; m = np.abs(a_out - a1) < 0.2 * a1
    return float(a1), float(np.median(e_out[m]))


def gradients(v, hx, hy):
    """per column: linear regression on the hodoscope position, (constant, d/dx, d/dy) [per cm]"""
    Xm = np.column_stack([np.ones(len(hx)), hx - np.median(hx), hy - np.median(hy)])
    C, *_ = np.linalg.lstsq(Xm, v, rcond=None)
    return C


# ---------------- fits ----------------
def fit_channels(xs, groups, sig, err, channels, n_iter=4, per_chip=True):
    """ln(signal) = a_group + b x + delta_(key, run); returns {ch: (b, eb, resid_rms)}, {key: delta}.
    delta: common run-to-run offset of the channels of one DRS chip (per_chip) or of all channels, taken as the median residual and
    kept free of a constant and a linear trend in x within each group, so it cannot absorb the normalisation or the attenuation.
    n_iter=1: no offsets."""
    G = sorted(set(groups)); gi = np.array([G.index(g) for g in groups])
    key = chip if per_chip else (lambda ch: "all")
    delta = {key(ch): np.zeros(len(xs)) for ch in channels}; res = {}
    Xm = np.column_stack([(gi == k).astype(float) for k in range(len(G))] + [xs])
    for _ in range(n_iter):
        resid = {}
        for ch in channels:
            y = np.log(sig[ch]) - delta[key(ch)]; w = (sig[ch] / err[ch]) ** 2; W = np.sqrt(w)
            c, *_ = np.linalg.lstsq(Xm * W[:, None], y * W, rcond=None)
            r = y - Xm @ c; chi2 = (r ** 2 * w).sum() / max(len(y) - Xm.shape[1], 1)
            cov = np.linalg.inv((Xm * w[:, None]).T @ Xm)
            res[ch] = (float(c[-1]), float(np.sqrt(cov[-1, -1]) * max(1.0, np.sqrt(chi2))), float(r.std()))
            resid.setdefault(key(ch), []).append(r)
        for k_, rs in resid.items():
            if len(rs) < 2 and per_chip: continue                    # a lone channel cannot define its chip's offset
            delta[k_] = detrend(delta[k_] + np.median(np.array(rs), axis=0), xs, gi)
    return res, delta


def detrend(dl, xs, gi):
    """remove a constant and a linear trend in x within each group"""
    dl = dl.copy()
    for k in set(gi):
        m = gi == k
        dl[m] -= np.polyval(np.polyfit(xs[m], dl[m], 1), xs[m]) if m.sum() >= 3 else dl[m].mean()
    return dl


def fit_sum(xs, groups, y, ey):
    """ln(y) = a_group + b x, weighted; b, error (scaled with sqrt(chi2/ndf) if > 1), chi2/ndf, residuals, coefficients"""
    G = sorted(set(groups)); gi = np.array([G.index(g) for g in groups])
    l = np.log(y); w = (y / ey) ** 2
    Xm = np.column_stack([(gi == k).astype(float) for k in range(len(G))] + [xs])
    c, *_ = np.linalg.lstsq(Xm * np.sqrt(w)[:, None], l * np.sqrt(w), rcond=None)
    r = l - Xm @ c; chi2 = float((r ** 2 * w).sum() / max(len(l) - Xm.shape[1], 1))
    cov = np.linalg.inv((Xm * w[:, None]).T @ Xm)
    return dict(inv_lambda=float(c[-1]), err=float(np.sqrt(cov[-1, -1]) * max(1.0, np.sqrt(chi2))), chi2ndf=chi2,
                resid=r.tolist(), resid_rms=float(r.std()), coef=c.tolist(), stat=float(np.median(ey / y)))


def fit_time(xs, groups, t, et):
    """time = t_group + slope x; slope [ns/m], its error (scaled with sqrt(chi2/ndf)), residual rms [ns]"""
    G = sorted(set(groups)); gi = np.array([G.index(g) for g in groups]); ok = np.isfinite(t) & np.isfinite(et) & (et > 0)
    if ok.sum() < 5: return np.nan, np.nan, np.nan
    Xm = np.column_stack([(gi[ok] == k).astype(float) for k in range(len(G))] + [xs[ok]])
    keep = Xm[:, :-1].sum(0) > 0; Xm = np.column_stack([Xm[:, :-1][:, keep], Xm[:, -1]])
    w = 1 / et[ok] ** 2; W = np.sqrt(w); c, *_ = np.linalg.lstsq(Xm * W[:, None], t[ok] * W, rcond=None)
    r = t[ok] - Xm @ c; chi2 = (r ** 2 * w).sum() / max(len(r) - Xm.shape[1], 1)
    cov = np.linalg.inv((Xm * w[:, None]).T @ Xm)
    return float(c[-1]), float(np.sqrt(cov[-1, -1]) * max(1.0, np.sqrt(chi2))), float(r.std())


def wmean(v, e):
    """weighted mean, its error scaled with sqrt(chi2/ndf) if > 1, chi2/ndf"""
    v, e = np.asarray(v, float), np.asarray(e, float)
    if not len(v): return np.nan, np.nan, np.nan
    w = 1 / e ** 2; m = (w * v).sum() / w.sum(); chi2 = (((v - m) / e) ** 2).sum() / max(len(v) - 1, 1)
    return float(m), float(np.sqrt(1 / w.sum() * max(1.0, chi2))), float(chi2)


def drift_fit(chs, types, il, eil, gy):
    """lateral drift of the shower relative to the module: 1/lambda_ch = 1/lambda_type + ky gy_ch over the given channels (gy: the
    channel's relative change per hodoscope cm of beam position); robust (3 sigma), 0.03 /m added to the errors for the channel-level
    scatter. ky: hodoscope cm per m of table travel."""
    types, il, eil, gy = map(np.asarray, (types, il, eil, gy))
    present = [t for t in TYPES if (types == t).any()]
    if len(chs) < 10: return None
    Am = np.column_stack([(types == t).astype(float) for t in present] + [gy]); m = np.ones(len(il), bool)
    w = 1 / (eil ** 2 + 0.03 ** 2)
    for _ in range(5):
        c, *_ = np.linalg.lstsq(Am[m] * np.sqrt(w[m])[:, None], il[m] * np.sqrt(w[m]), rcond=None)
        r = il - Am @ c; s = 1.4826 * np.median(np.abs(r[m])); m = np.abs(r) < 3 * max(s, 1e-6)
    cov = np.linalg.inv((Am[m] * w[m][:, None]).T @ Am[m]); chi2 = ((r[m] ** 2) * w[m]).sum() / max(m.sum() - Am.shape[1], 1)
    err = np.sqrt(np.diag(cov)) * max(1.0, np.sqrt(chi2))
    return dict(ky=float(c[-1]), eky=float(err[-1]), n_out=int((~m).sum()),
                types={t: dict(inv_lambda=float(c[i]), err=float(err[i]), n=int(((types == t) & m).sum())) for i, t in enumerate(present)})


# ---------------- FERS (datasets with calorimeter FERS) ----------------
def fers_type(m):
    return "Sci" if not m[2] else ("CerQuartz" if m[3] else "CerPlastic")


def sci_centroid(d, sel):
    """per event: scintillation-weighted tower row (FERS), the shower position in the module"""
    m = d["meta"]; sci = m[:, 2] == 0; v = np.clip(d["vals"][sel][:, sci], 0, None)
    return (v * m[sci, 5]).sum(1) / np.maximum(v.sum(1), 1e-9)


def matching(ds):
    """vertical drift of the shower in the module with the table position: run medians of the FERS shower centroid (tower rows) vs X,
    straight line over the main run set; converted to hodoscope cm with d(centroid)/d(hodoscope y) from the beam spread within the runs.
    Sets ds.shift (hodoscope cm added to hodoscope y per run) and ds.yc (common window centre)."""
    cy, k, hys = {}, [], {}
    for r in ds.runs:
        d = np.load(f"{DATA}/fers/fers_{r}.npz"); sel = pid_hodo(d, ds.daq(r))
        c = sci_centroid(d, sel); hy = d["TTU_Hodo_Y"][sel]
        cy[r] = float(np.median(c)); k.append(np.polyfit(hy, c, 1)[0]); hys[r] = hy
    use = ds.runsets["main"]; xs = np.array([ds.x[r] for r in use])
    b = np.polyfit(xs, [cy[r] for r in use], 1)[0]; k = float(np.median(k)); xref = float(xs.mean())
    ds.shift = {r: float(b * (ds.x[r] - xref) / k) for r in ds.runs}
    ds.yc = float(np.median([np.median(hys[r] + ds.shift[r]) for r in use]))
    return dict(centroid_median=[cy[r] for r in ds.runs], centroid_slope_rows_per_m=float(b), rows_per_hodo_cm=k,
                hodo_cm_per_m=float(b / k), shift_cm=[ds.shift[r] for r in ds.runs], yc=ds.yc)


def load_fers(ds, run, mode, with_mcp):
    """FERS calibrated values (mean = light per event; no in-time selection needed: charge over a long gate), positrons"""
    d = np.load(f"{DATA}/fers/fers_{run}.npz"); daq = ds.daq(run)
    s = pid_hodo(d, daq)
    if with_mcp: s &= mcp_ok(d, daq)
    sel = s & beam_window(d, s, ds, run, mode)
    v = d["vals"][sel]
    C = gradients(v, d["TTU_Hodo_X"][sel], d["TTU_Hodo_Y"][sel])
    names = [f"FERS_B{int(m[0])}_{int(m[1])}" for m in d["meta"]]
    meta = {n: dict(type=fers_type(m), tx=float(m[4]), ty=float(m[5])) for n, m in zip(names, d["meta"])}
    with np.errstate(divide="ignore", invalid="ignore"):
        gx, gy = C[1] / C[0], C[2] / C[0]
    return dict(names=names, meta=meta, n=int(sel.sum()), mean=v.mean(0), err=v.std(0) / np.sqrt(len(v)), gx=gx, gy=gy,
                lg_max=np.percentile(d["lg"][sel], 99.9, axis=0), cy=float(np.median(sci_centroid(d, sel))), vals=v)
