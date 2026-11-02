"""Electronics part of the timing resolution, from the reference dumps (dump_refs.py) of REF_RUNS:
  mcp_same_chip   width of MCP_US_1 - MCP_DS_1: two MCPs on the same DRS chip (no group alignment); in 2025 both MCPs give ~125 ADC,
                  so one MCP (with its digitization) is the width / sqrt(2); in 2026 the upstream MCP gives ~80 ADC against ~160
  same_chip_boards  the same-chip MCP pair on every board that has both MCPs (2025: all four boards)
  copies          every pair of digitized copies of the same MCP signal on different DRS boards (after the group alignment): no MCP
                  jitter; what the alignment of two groups adds to real pulses, plus the digitization of two large pulses. 2025: each
                  MCP on all four boards (6 pairs each); 2026: two copies
  per_board       with three or more copies (2025): one term per board from sigma^2(i - j) = a_i^2 + a_j^2
  groups          reference (Channel8) times of two groups on the same board, per group = rms of the difference / sqrt(2), for the
                  integer sample (former definition), the interpolated crossing (framework) and a line fit; larger than what reaches
                  a channel time, since a pulse and the reference of its own group share part of the cell-timing error of the chip
  cells           group-to-group difference vs the DRS cell of the edge (non-uniform cell widths), one pair of CELL_RUN, and the width
                  after a per-cell correction learned on half of the events and applied to the other half
  channels        channel - MCP_DS_1 width of every channel with integer and interpolated references, by the group relation to the MCP
Widths: gaussian core (2-sigma clipping) for continuous distributions; robust rms where the integer reference makes them discrete.
Usage (repository root): python3 studies/timing/electronics.py      Output: $TIMING_DATA/electronics.json"""
import json, sys
import numpy as np
from config import DATA, REF_RUNS, CELL_RUN, TS_PS
sys.path.insert(0, ".")
from configs.selection_config import get_particle_selection

REF_OFFSET = 790.0                     # framework: TS_cfd_ref = 790 + TS_cfd - ref_TS
VARIANTS = ("ref_int", "ref_fw", "ref_fit")
VLAB = {"ref_int": "integer", "ref_fw": "interpolated", "ref_fit": "line fit"}


def core(x, nit=8, nmin=50):
    """gaussian core width (iterative 2-sigma clipping, x 1.137), its error and the fraction beyond 3 sigma"""
    x = x[np.isfinite(x)]
    if len(x) < nmin: return np.nan, np.nan, np.nan
    m, s = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
    for _ in range(nit):
        k = np.abs(x - m) < 2 * s
        if k.sum() < nmin or s <= 0: return np.nan, np.nan, np.nan
        m, s = x[k].mean(), x[k].std() * 1.137
    return float(s), float(s / np.sqrt(2 * k.sum())), float(np.mean(np.abs(x - m) > 3 * s))


def rrms(x):
    """robust rms: within +-5 x max(MAD sigma, 100 ps) of the median (works for the discrete integer differences)"""
    m = np.median(x); s = 1.4826 * np.median(np.abs(x - m)); k = np.abs(x - m) < 5 * max(s, 100.0)
    return float(x[k].std())


def gname(ch): return ch.rsplit("_Channel", 1)[0]
def bname(ch): return gname(ch).rsplit("_Group", 1)[0]


def run_analysis(run):
    d = np.load(f"{DATA}/refs/refs_{run}.npz"); groups = list(d["groups"]); gi = {g: i for i, g in enumerate(groups)}
    mcps = dict(x.split("=") for x in d["mcps"]); mk = list(mcps); phys = list(d["phys"])
    pid = np.ones(len(d["mcp_T"]), bool)
    for det, fire in get_particle_selection("electron", run).items():
        k = f"is_{det}_fired" if fire else f"is_{det}_vetoed"
        if k in d: pid &= d[k].astype(bool)
    T = {m: d["mcp_T"][:, i] for i, m in enumerate(mk)}; A = {m: d["mcp_A"][:, i] for i, m in enumerate(mk)}
    ref = {v: d[v] for v in VARIANTS}
    def aligned(t, ch, v, s):                    # channel or MCP time after the alignment of its group [ps], events s
        r = ref[v][s, gi[gname(ch)]]; return np.where((t > 0) & (r > 0), (REF_OFFSET + t - r) * TS_PS, np.nan)
    sel = pid & (A["MCP_DS_1"] > 50)          # the z-scan MCP requirement
    out = dict(run=run, mcps=mcps, amp={m: float(np.median(A[m][sel])) for m in mk})
    # two MCPs, same chip
    s = pid & (T["MCP_US_1"] > 0) & (T["MCP_DS_1"] > 0) & (A["MCP_US_1"] > 50) & (A["MCP_DS_1"] > 50)
    x = (T["MCP_US_1"][s] - T["MCP_DS_1"][s]) * TS_PS; w, e, tail = core(x)
    out["mcp_same_chip"] = dict(sigma=w, err=e, tail=tail, n=int(s.sum()), hist=np.histogram(x - np.median(x), np.arange(-600, 601, 10))[0].tolist())
    # two MCPs on the same chip, on every board that has both (2025: copies on all four boards)
    out["same_chip_boards"] = {}
    for k in sorted({m.rsplit("_", 1)[1] for m in mk}):
        u, w_ = f"MCP_US_{k}", f"MCP_DS_{k}"
        if u in T and w_ in T:
            s = pid & (T[u] > 0) & (T[w_] > 0) & (A[u] > 50) & (A[w_] > 50)
            out["same_chip_boards"][bname(mcps[w_])] = dict(sigma=core((T[u][s] - T[w_][s]) * TS_PS)[0], n=int(s.sum()))
    # copies of the same MCP on two boards, after alignment: every pair of copies
    out["copies"] = {}
    for det in ("MCP_DS", "MCP_US"):
        cs = sorted(m for m in mk if m.startswith(det + "_"))
        for i, a in enumerate(cs):
            for b in cs[i + 1:]:
                s = pid & (A[a] > 50) & (A[b] > 50); rec = dict(boards=[bname(mcps[a]), bname(mcps[b])], n=int(s.sum()))
                if s.sum() < 100: continue
                for v in VARIANTS:
                    x = aligned(T[a][s], mcps[a], v, s) - aligned(T[b][s], mcps[b], v, s); x = x[np.isfinite(x)]
                    rec[v] = dict(core=core(x)[0], rms=rrms(x) if len(x) > 50 else np.nan, tail=core(x)[2])
                    if v == "ref_fw": rec["hist"] = np.histogram(x - np.median(x), np.arange(-600, 601, 10))[0].tolist()
                out["copies"][f"{a}-{b}"] = rec
    # with three or more copies: one term per board, sigma^2(i - j) = a_i^2 + a_j^2 (least squares over the pairs of both MCPs)
    out["per_board"] = {}
    for v in ("ref_int", "ref_fw"):
        prs = [(c["boards"], c[v]["core"]) for c in out["copies"].values() if v in c and np.isfinite(c[v]["core"])]
        bl = sorted({b for p_, _ in prs for b in p_})
        if len(bl) >= 3 and len(prs) >= len(bl):
            M = np.array([[1.0 if b in p_ else 0.0 for b in bl] for p_, _ in prs]); y = np.array([w_ ** 2 for _, w_ in prs])
            a2 = np.linalg.lstsq(M, y, rcond=None)[0]
            out["per_board"][v] = {b: float(np.sqrt(max(x, 0))) for b, x in zip(bl, a2)}
    # group alignment, groups on the same board
    pairs = [(a, b) for i, a in enumerate(groups) for b in groups[i + 1:] if a.rsplit("_Group", 1)[0] == b.rsplit("_Group", 1)[0]]
    gr = {}
    for v in VARIANTS:
        vals = []
        for a, b in pairs:
            ra, rb = ref[v][:, gi[a]], ref[v][:, gi[b]]; ok = (ra > 0) & (rb > 0)
            if ok.sum() > 100: vals.append(rrms((ra[ok] - rb[ok]) * TS_PS))
        gr[v] = dict(pair=float(np.median(vals)), per_group=float(np.median(vals) / np.sqrt(2)), p10=float(np.percentile(vals, 10)), p90=float(np.percentile(vals, 90)), n=len(vals))
    a, b = pairs[0]; ok = (ref["ref_int"][:, gi[a]] > 0) & (ref["ref_int"][:, gi[b]] > 0)
    gr["example"] = dict(pair=[a, b], **{v: np.histogram((ref[v][ok, gi[a]] - ref[v][ok, gi[b]]) * TS_PS - np.median((ref[v][ok, gi[a]] - ref[v][ok, gi[b]]) * TS_PS),
                                                          np.arange(-500, 501, 10))[0].tolist() for v in VARIANTS})
    out["groups"] = gr
    # DRS cell pattern of the group-to-group difference
    if run == CELL_RUN:
        sc = d["start_cell"].astype(int); ta, tb = ref["ref_fw"][ok, gi[a]], ref["ref_fw"][ok, gi[b]]
        dt = (ta - tb) * TS_PS; dt -= np.median(dt)
        ca = (sc[ok, gi[a]] + np.floor(ta).astype(int)) % 1024; cb = (sc[ok, gi[b]] + np.floor(tb).astype(int)) % 1024
        prof = [float(np.median(dt[ca == c])) if np.sum(ca == c) > 5 else None for c in range(1024)]
        idx = np.arange(len(dt)); tr, te = idx % 2 == 0, idx % 2 == 1
        pa = np.array([np.median(dt[tr & (ca == c)]) if np.sum(tr & (ca == c)) > 3 else 0 for c in range(1024)])
        pb = np.array([np.median((dt - pa[ca])[tr & (cb == c)]) if np.sum(tr & (cb == c)) > 3 else 0 for c in range(1024)])
        out["cells"] = dict(pair=[a, b], profile=prof, n=int(ok.sum()), raw=rrms(dt[te]), corrected=rrms(dt[te] - pa[ca[te]] - pb[cb[te]]))
    # channel widths with integer / interpolated references, by group relation to the reference MCP
    rel = {}
    tM = T["MCP_DS_1"]; okM = pid & (tM > 0) & (A["MCP_DS_1"] > 50); gM = gname(mcps["MCP_DS_1"])
    for j, ch in enumerate(phys):
        t, a = d["phys_T"][:, j], d["phys_A"][:, j]; s = okM & (t > 0) & (a > 40)
        if s.sum() < 300: continue
        tv = {v: aligned(t[s], ch, v, s) - aligned(tM[s], mcps["MCP_DS_1"], v, s) for v in ("ref_int", "ref_fw")}
        x0 = tv["ref_int"][np.isfinite(tv["ref_int"])]
        if len(x0) < 300: continue
        h, e = np.histogram(x0, np.arange(x0.min(), x0.max() + 400, 400)); t0 = e[np.argmax(h)] + 200
        inn = np.abs(tv["ref_int"] - t0) < 3000
        if inn.sum() < 200: continue
        r = "same group" if gname(ch) == gM else ("same board" if bname(ch) == bname(mcps["MCP_DS_1"]) else "other board")
        w = {v: core(tv[v][inn])[0] for v in tv}
        if all(np.isfinite(list(w.values()))): rel.setdefault(r, []).append((ch, w["ref_int"], w["ref_fw"]))
    out["channels"] = {r: dict(n=len(v), int=float(np.median([x[1] for x in v])), interp=float(np.median([x[2] for x in v])),
                               gain=float(np.median([np.sqrt(max(x[1] ** 2 - x[2] ** 2, 0)) for x in v]))) for r, v in rel.items()}
    return out


def main():
    res = {}
    for period, runs in REF_RUNS.items():
        rr = [run_analysis(r) for r in runs]; res[period] = dict(runs=rr)
        r0 = rr[0]
        print(f"{period}: MCP amplitudes {r0['amp']}")
        for r in rr:
            c = r["copies"].get("MCP_DS_0-MCP_DS_1", {})
            print("    copies: " + "; ".join(f"{k} ({'/'.join(b.replace('DRS_', '') for b in x['boards'])}) {x['ref_int']['core']:.0f}->{x['ref_fw']['core']:.0f}" for k, x in r["copies"].items()))
            if r["per_board"]: print("    per board (interpolated): " + ", ".join(f"{b.replace('DRS_', '')} {x:.0f}" for b, x in r["per_board"].get("ref_fw", {}).items()))
            print(f"  run {r['run']}: same-chip MCP pair {r['mcp_same_chip']['sigma']:.0f} ps (>3σ {100 * r['mcp_same_chip']['tail']:.0f}%);"
                  f" DS copies on boards {c.get('boards')}: " + ", ".join(f"{VLAB[v]} {c[v]['core']:.0f}" for v in VARIANTS if v in c) +
                  "; groups on one board, per group: " + ", ".join(f"{VLAB[v]} {r['groups'][v]['per_group']:.0f}" for v in VARIANTS) +
                  "; channels int->interp: " + ", ".join(f"{k} {v['int']:.0f}->{v['interp']:.0f} ({v['n']})" for k, v in r["channels"].items()))
            if "cells" in r: print(f"    cells ({r['cells']['pair']}): {r['cells']['raw']:.0f} ps -> {r['cells']['corrected']:.0f} ps with a per-cell correction")
    json.dump(res, open(f"{DATA}/electronics.json", "w"))


if __name__ == "__main__":
    main()
