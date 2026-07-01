import argparse
from configs.run_config import run_number, firstEvent, lastEvent, jsonFile
from utils.data_loader import getRunInfo


def get_args():
    parser = argparse.ArgumentParser(description="prepare CaloX data analysis")
    parser.add_argument("--run", default=run_number,
                        type=int, help="Run number")
    parser.add_argument(
        "--first-event", default=firstEvent, type=int, help="First event")
    parser.add_argument(
        "--last-event", default=lastEvent, type=int, help="Last event")
    parser.add_argument(
        "--json-file", default=jsonFile, type=str, help="JSON file with data files")
    parser.add_argument(
        "--sequences", default=None, nargs="+", metavar="SEQ",
        help="DQM sequences to run (default: all enabled-by-default sequences). "
             "Example: --sequences monitor_conditions drs_waveforms drs_stats")
    parser.add_argument(
        "--particle", default=None, type=str, metavar="TYPE",
        help="Select a specific particle type before analysis "
             "(e.g. pion, electron, proton, muon). Default: no selection.")
    parser.add_argument(
        "--jsroot", action="store_true", default=False,
        help="Save canvases to ROOT file and generate JSROOT-powered HTML (no PNG files)")
    parser.add_argument(
        "--hole-veto", action="store_true", default=False,
        help="Apply the hole veto selection before analysis.")
    parser.add_argument(
        "--mcp-clean", action="store_true", default=False,
        help="Require the mcp_clean selection (all MCPs fire) before analysis.")
    parser.add_argument(
        "--channels", default=None, type=str, metavar="FILE",
        help="JSON file mapping labels to channel names "
             "(e.g. data/channel_maps/testingfibers.json). "
             "When given, only the listed channels are histogrammed.")
    args = parser.parse_args()

    print(f"Run number: {args.run}")
    print(f"First event: {args.first_event}")
    print(f"Last event: {args.last_event}")
    print(f"JSON file: {args.json_file}")

    btype, benergy = getRunInfo(args.run)
    print(f"Beam type: {btype}")
    print(f"Beam energy: {benergy} GeV")
    return args
