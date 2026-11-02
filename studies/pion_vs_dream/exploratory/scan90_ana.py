import numpy as np, json
rl = json.load(open("/data1/yofeng/CaloX/CaloXDataAnalysis/data/Runlist.json"))
RUNS = [1499, 1500, 1501, 1502, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1513]
def num(s): return float(str(s).replace("mm", "").strip())
def peak(x):
    m = np.median(x); s = 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 1.5 * s
        if w.sum() < 50: break
        m2, s2 = x[w].mean(), x[w].std() / 0.74
        if abs(m2 - m) < 1e-4: break
        m, s = m2, s2
    return m, s / np.sqrt(max(np.sum(np.abs(x - m) < 1.5 * s), 1))
rows = []
for r in RUNS:
    try: d = np.load(f"scan90_{r}.npz")
    except FileNotFoundError: continue
    S, Cq, Cp = d["FERS_S"], d["FERS_Cq"], d["FERS_Cp"]
    if len(S) < 2000: continue
    sh = S > 2.0
    core = S > 0.5 * np.median(S[sh])                   # main peak region, for the peak fit only
    row = {"run": r, "X": num(rl[str(r)]["Table position X"]), "Y": num(rl[str(r)]["Table position Y"]),
           "deg": rl[str(r)]["Table Degree"], "N": int(sh.sum())}
    for k, x in (("S", S), ("Cq", Cq), ("Cp", Cp)):
        row[f"{k}_med"] = float(np.median(x[sh])); row[f"{k}_pk"], row[f"{k}_pkerr"] = map(float, peak(x[core]))
    rows.append(row)
print(f"{'run':>5s} {'X':>7s} {'Y':>6s} {'deg':>5s} {'N':>6s} | {'S med':>6s} {'S peak':>7s} | {'Cq med':>6s} {'Cq peak':>7s} | {'Cp med':>6s} {'Cp peak':>7s}")
for r in sorted(rows, key=lambda t: (t["Y"], t["X"])):
    print(f"{r['run']:5d} {r['X']:7.1f} {r['Y']:6.0f} {r['deg']:>5s} {r['N']:6d} | {r['S_med']:6.2f} {r['S_pk']:7.2f} | {r['Cq_med']:6.2f} {r['Cq_pk']:7.2f} | {r['Cp_med']:6.2f} {r['Cp_pk']:7.2f}")
print("\nexponential fit per Y (runs at X = -168, -218, -268 mm): signal change per 100 mm of X, and |lambda| = 100 mm / |ln ratio|")
for Y in sorted({r["Y"] for r in rows}):
    g = [r for r in rows if r["Y"] == Y and r["X"] in (-168.0, -218.0, -268.0)]
    if len(g) < 2: continue
    xs = np.array([r["X"] for r in g])
    out = []
    for k in ("S", "Cq", "Cp"):
        for kind in ("pk", "med"):
            y = np.log([r[f"{k}_{kind}"] for r in g]); sl = np.polyfit(xs, y, 1)[0] * 100
            out.append(f"{k} {kind} {np.expm1(sl) * 100:+5.1f}%/100mm (|lam| {abs(100 / sl) / 1000 if sl else np.inf:5.2f} m)")
    print(f"  Y={Y:5.0f} ({len(g)} runs, X {sorted(xs.tolist())}): " + " | ".join(out[i] for i in (0, 2, 4)))
    print(f"  {'':12s}  medians:            " + " | ".join(out[i] for i in (1, 3, 5)))
json.dump(rows, open("scan90_2025.json", "w"), indent=1)
