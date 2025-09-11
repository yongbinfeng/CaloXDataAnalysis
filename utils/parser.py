import argparse
from runconfig import runNumber, firstEvent, lastEvent, jsonFile


def get_args():
    parser = argparse.ArgumentParser(description="prepare CaloX data analysis")
    parser.add_argument("--run", default=runNumber,
                        type=int, help="Run number")
    parser.add_argument(
        "--first-event", default=firstEvent, type=int, help="First event")
    parser.add_argument(
        "--last-event", default=lastEvent, type=int, help="Last event")
    parser.add_argument(
        "--jsonFile", default=jsonFile, type=str, help="JSON file with data files")
    args = parser.parse_args()

    print(f"Run number: {args.run}")
    print(f"First event: {args.first_event}")
    print(f"Last event: {args.last_event}")
    print(f"JSON file: {args.jsonFile}")
    return args
