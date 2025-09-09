import argparse
from runconfig import runNumber, firstEvent, lastEvent
from utils.utils import getRunInfo


def get_args():
    parser = argparse.ArgumentParser(description="prepare CaloX data analysis")
    parser.add_argument("--run", default=runNumber,
                        type=int, help="Run number")
    parser.add_argument(
        "--first-event", default=firstEvent, type=int, help="First event")
    parser.add_argument(
        "--last-event", default=lastEvent, type=int, help="Last event")
    args = parser.parse_args()

    print(f"Run number: {args.run}")
    print(f"First event: {args.first_event}")
    print(f"Last event: {args.last_event}")

    btype, benergy = getRunInfo(args.run)
    print(f"Beam type: {btype}")
    print(f"Beam energy: {benergy} GeV")
    return args
