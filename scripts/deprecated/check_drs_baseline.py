"""
Compare DRS baseline subtraction methods for a single service-DRS detector.

Books the baseline-subtracted waveform (blsub) vs time slice two ways:
  - constant : median(0, 200) subtracted as a flat offset (current default)
  - linear   : two-anchor linear baseline (clean pre-pulse + post-pulse medians),
               which removes the post-pulse pedestal droop seen on big pulses.

Produces a side-by-side 2D comparison PNG plus tail-distribution statistics, so
you can eyeball how much the linear baseline cleans up the right-side drift.

Usage
-----
  python scripts/check_drs_baseline.py --run 1804
  python scripts/check_drs_baseline.py --run 1804 --det HoleVeto
  python scripts/check_drs_baseline.py --run 1804 --det Cer474 --nfiles 3
  python scripts/check_drs_baseline.py --run 1804 --pre 0 120 --post 800 1000
"""

import os
import json
import argparse

import ROOT
from utils.root_setup import setup_root
from channels.channel_map import get_service_drs_channels
from variables.drs import subtract_baseline


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=int, required=True, help="Run number")
    ap.add_argument("--det", default="HoleVeto",
                    help="Service-DRS detector name (default: HoleVeto)")
    ap.add_argument("--json", default="data/datafiles_local.json",
                    help="JSON file mapping run -> file list")
    ap.add_argument("--nfiles", type=int, default=2,
                    help="Number of segment files to chain (default: 2)")
    ap.add_argument("--pre", type=int, nargs=2, default=(0, 120),
                    metavar=("LO", "HI"), help="Pre-pulse anchor window")
    ap.add_argument("--post", type=int, nargs=2, default=(800, 1000),
                    metavar=("LO", "HI"), help="Post-pulse anchor window")
    ap.add_argument("--outdir", default="results/baseline_check",
                    help="Output directory for the comparison PNG")
    args = ap.parse_args()

    setup_root(n_threads=10, batch_mode=True, load_functions=True)

    channel = get_service_drs_channels(args.run).get(args.det)
    if channel is None:
        raise SystemExit(
            f"Detector '{args.det}' has no service-DRS channel for run {args.run}. "
            f"Available: {list(get_service_drs_channels(args.run))}")
    print(f"Run {args.run}  detector {args.det}  channel {channel}")

    files = json.load(open(args.json))[str(args.run)][:args.nfiles]
    chain = ROOT.TChain("EventTree")
    for f in files:
        chain.Add(f)
    print(f"Chained {len(files)} file(s), {chain.GetEntries()} entries")

    # Service channels are flipped (negative-going pulse -> positive blsub).
    # Build both baselines in one node so the two blsub columns coexist:
    #   blsub_const : constant median(0,200) via subtract_baseline(linear=False)
    #   blsub_lin   : linear two-anchor line, defined inline below
    flip = [channel]
    rdf = ROOT.RDataFrame(chain)
    node = subtract_baseline(rdf, [channel], drs_channels_to_flip=flip, linear=False)
    node = node.Define("blsub_const", f"{channel}_blsub")
    node = node.Define("bl_pre",  f"compute_baseline_median({channel}, {args.pre[0]}, {args.pre[1]})")
    node = node.Define("bl_post", f"compute_baseline_median({channel}, {args.post[0]}, {args.post[1]})")
    pre_c = 0.5 * (args.pre[0] + args.pre[1])
    post_c = 0.5 * (args.post[0] + args.post[1])
    node = node.Define("bl_line",
                       f"bl_pre + ((bl_post - bl_pre) / double({post_c} - {pre_c})) * (ts - {pre_c})")
    node = node.Define("blsub_lin", f"-({channel} - bl_line)")
    node = node.Define("v900_const", "blsub_const[900]")
    node = node.Define("v900_lin", "blsub_lin[900]")

    YN, YLO, YHI = 500, -500, 3000
    h2c = node.Histo2D(("h2c", f"{args.det} constant baseline (current);Time Slice;blsub ADC",
                        1024, 0, 1024, YN, YLO, YHI), "ts", "blsub_const")
    h2l = node.Histo2D(("h2l", f"{args.det} linear two-anchor baseline;Time Slice;blsub ADC",
                        1024, 0, 1024, YN, YLO, YHI), "ts", "blsub_lin")
    hc = node.Histo1D(("hc", "blsub[900] const", 400, -300, 1200), "v900_const")
    hl = node.Histo1D(("hl", "blsub[900] linear", 400, -300, 1200), "v900_lin")
    total = node.Count()
    band_c = node.Filter("v900_const > 400").Count()
    band_l = node.Filter("v900_lin > 400").Count()

    ROOT.RDF.RunGraphs([h2c, h2l, hc, hl, total, band_c, band_l])

    T = total.GetValue()
    print(f"\nTotal events: {T}")
    print(f"Tail-band events (blsub[900] > 400):")
    print(f"  constant: {band_c.GetValue()} ({100*band_c.GetValue()/T:.3f}%)")
    print(f"  linear  : {band_l.GetValue()} ({100*band_l.GetValue()/T:.3f}%)")
    print(f"blsub[900]  constant: mean={hc.GetMean():7.2f} rms={hc.GetRMS():6.2f}")
    print(f"blsub[900]  linear  : mean={hl.GetMean():7.2f} rms={hl.GetRMS():6.2f}")

    os.makedirs(args.outdir, exist_ok=True)
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "c", 1500, 650)
    c.Divide(2, 1)
    c.cd(1); ROOT.gPad.SetLogz(); h2c.Draw("colz")
    c.cd(2); ROOT.gPad.SetLogz(); h2l.Draw("colz")
    out = os.path.join(args.outdir, f"baseline_compare_run{args.run}_{args.det}.png")
    c.SaveAs(out)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
