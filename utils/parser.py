import argparse
from runconfig import runNumber, firstEvent, lastEvent


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
    return args
