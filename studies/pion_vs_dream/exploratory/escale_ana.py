import numpy as np
runs = {1423: ("e+ 20", "electron"), 1416: ("e+ 30", "electron"), 1514: ("e+ 40", "electron"),
        1526: ("e+ 60", "electron"), 1527: ("e+ 100", "electron"), 1465: ("pi+ 80", "pion")}
need = {"electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"],
        "pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]}
edges = np.array([0.05, 0.1, 0.2, 0.4, 0.7, 1, 1.5, 2, 3, 4.5, 6.5, 9, 13, 18, 25, 35, 50, 70])
data = {}
for run, (lab, part) in runs.items():
    d = np.load(f"escale_{run}.npz")
    sel = np.ones(len(d["t0_S"]) if "t0_S" in d.files else len(d[d.files[3]]), bool)
    for f in need[part]:
        sel &= d[f].astype(bool)
    ids = [int(r[0]) for r in d["meta"]]
    S = np.stack([d[f"t{i}_S"][sel] for i in ids], 1); C = np.stack([d[f"t{i}_C"][sel] for i in ids], 1)
    SLG = np.stack([d[f"t{i}_S_lg"][sel] for i in ids], 1); CLG = np.stack([d[f"t{i}_C_lg"][sel] for i in ids], 1)
    data[run] = (lab, S, C, SLG, CLG)
    tot_S, tot_C = S.sum(1), C.sum(1)
    hot = np.argmax(S, 1); hotS = S[np.arange(len(S)), hot]
    lg_frac = (S * SLG).sum() / S.sum()
    print(f"{lab:7s} N={sel.sum():6d}  <S>={np.median(tot_S):6.1f} <C>={np.median(tot_C):6.1f}  hottest S tower: median {np.median(hotS):5.1f} GeV "
          f"({np.median(hotS / tot_S):.2f} of S) | S energy in low-gain channels {lg_frac:.2f}, C {(C * CLG).sum() / C.sum():.2f}")
print()
print("Sum S / sum C per tower, in bins of the tower signal (S+C)/2 [GeV]:")
hdr = "bin [GeV]      " + " ".join(f"{data[r][0]:>14s}" for r in runs)
print(hdr)
for lo, hi in zip(edges[:-1], edges[1:]):
    cells = []
    for run in runs:
        lab, S, C, SLG, CLG = data[run]
        a = 0.5 * (S + C); b = (a >= lo) & (a < hi)
        n = b.sum()
        cells.append(f"{S[b].sum() / C[b].sum():6.3f} ({n:6d})" if n >= 50 else " " * 14)
    print(f"{lo:5.2f}-{hi:5.2f}    " + " ".join(cells))
print()
print("Same, split by whether the S channel is on low gain (switch at 7500 HG counts):")
for run in runs:
    lab, S, C, SLG, CLG = data[run]
    a = 0.5 * (S + C)
    out = []
    for lo, hi in ((3, 6.5), (6.5, 13), (13, 25), (25, 50)):
        b = (a >= lo) & (a < hi)
        hg, lg = b & ~SLG & ~CLG, b & SLG & ~CLG
        f = lambda m: f"{S[m].sum() / C[m].sum():.3f}/{m.sum()}" if m.sum() >= 30 else "-"
        out.append(f"{lo}-{hi}: HG {f(hg)} LG {f(lg)}")
    print(f"  {lab:7s} " + " | ".join(out))
