"""Raw FERS HG of every calorimeter channel (no pedestal subtraction, no calibration) for the gain fit, in the channel order of the
framework map.
Usage (repository root): python3 studies/pe_yield/dump_raw.py RUN NEVENTS      Output: $PE_YIELD_DATA/raw/raw_<RUN>.npz"""
import json, os, sys
sys.path.insert(0, ".")
import numpy as np
import ROOT
from configs.run_config import jsonFile
from channels.maps.fers import build_fers_boards
from config import DATA

ROOT.gROOT.SetBatch(True)
run, n = int(sys.argv[1]), int(sys.argv[2])
files = json.load(open(jsonFile))[str(run)]
chain = ROOT.TChain("EventTree")
for f in (files if isinstance(files, list) else [files]): chain.Add(f)
boards = build_fers_boards(run_number=run)
meta = [(b.board_no, ch.channel_no, int(ch.isCer), int(ch.isQuartz), ch.i_tower_x, ch.i_tower_y) for b in boards.values() for ch in b]
cols = sorted({f"FERS_Board{m[0]}_energyHG" for m in meta})
d = ROOT.RDataFrame(chain).Range(n).AsNumpy(cols)
arr = {c: np.stack([np.asarray(x) for x in d[c]]) for c in cols}
hg = np.stack([arr[f"FERS_Board{m[0]}_energyHG"][:, m[1]] for m in meta], 1).astype(np.int16)
os.makedirs(f"{DATA}/raw", exist_ok=True)
np.savez(f"{DATA}/raw/raw_{run}.npz", hg=hg, meta=np.array(meta, dtype=np.int16))
print("done", run, hg.shape)
