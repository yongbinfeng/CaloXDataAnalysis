import json, numpy as np
T = json.load(open("test12.json")); pos = {int(k): v for k, v in T["pos"].items()}; ref = pos[1409]
SEL = {"electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"],
       "pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]}
cen = {r: (v["ex"], v["ey"]) for r, v in pos.items() if r != 1399}; x0, y0 = cen[1409]
offr = [r for r in cen if r != 1409]; rad_d = np.mean([np.hypot(cen[r][0]-x0, cen[r][1]-y0) for r in offr])
for key, isc in (("S", 0), ("C", 1)):
    out = {}
    for run, part in ((1432, "pion"), (1409, "electron")):
        d = np.load(f"chan_{run}.npz"); meta = d["meta"]; E = float(d["energy"]); sel = np.ones(len(d["vals"]), bool)
        for f in SEL[part]: sel &= d[f].astype(bool)
        V = np.nan_to_num(d["vals"][sel]); Ssum = V[:, meta[:, 2] == 0].sum(1)
        V = V[(Ssum > 0.85*E) & (Ssum < 1.25*E)] if part == "electron" else V[Ssum > 0.3*E]
        m = meta[:, 2] == isc; x, y = meta[m, 4]*1.2, meta[m, 5]*1.6; v = V[:, m].mean(0)
        xs, ys = (x*v).sum()/v.sum(), (y*v).sum()/v.sum(); dx, dy = x-xs+x0, y-ys+y0
        R_nn = np.array([pos[min(cen, key=lambda r: (cen[r][0]-a)**2 + (cen[r][1]-b)**2)][key]/ref[key] for a, b in zip(dx, dy)])
        rr = np.hypot(dx-x0, dy-y0); rad_off = np.mean([pos[r][key]/ref[key] for r in offr])
        R_rad = 1 + (rad_off-1)*np.clip(rr/rad_d, 0, 1.5)
        out[part] = (len(V), V[:, m].sum(1).mean(), (v*R_nn).sum()/v.sum(), (v*R_rad).sum()/v.sum(), [v[rr > R].sum()/v.sum() for R in (3, 6, 10)])
    p, e = out["pion"], out["electron"]
    print(f"{key}: N pion {p[0]}, e+ {e[0]} (e+ {key} sum {e[1]:.1f} GeV) | off-centre rel. response {rad_off:.3f} at r~{rad_d:.1f} cm")
    print(f"     energy fraction beyond 3/6/10 cm: pion {np.round(p[4], 3)}, e+ {np.round(e[4], 3)}")
    print(f"     weighted response, nearest run: pion {p[2]:.4f}, e+ {e[2]:.4f} -> pion/e+ {p[2]/e[2]:.4f};  radial-linear: pion {p[3]:.4f}, e+ {e[3]:.4f} -> {p[3]/e[3]:.4f}")
