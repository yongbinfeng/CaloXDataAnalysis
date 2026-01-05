# read csv file
import csv
import sys
import os
import re
import argparse
import json
from datetime import datetime


def read_csv(file_path):
    with open(file_path, mode='r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        data = [row for row in reader]
        print(f"Read {len(data)} rows from {file_path}")
        print(data)
    return data


def write_json(data, file_path):
    # instead of a list, make data a dictionary with RunNumber as key
    # if runNumber can not be converted to int, drop that row
    data = {row['RunNumber']: row for row in data if row['RunNumber'].isdigit()}
    with open(file_path, mode='w', encoding='utf-8') as jsonfile:
        json.dump(data, jsonfile, indent=4)


def testGetRunInfo(runNum=None):
    from utils.utils import getRunInfo
    if runNum == None:
        from runconfig import runNumber
        runNum = runNumber
    runInfo = getRunInfo(runNumber)
    print(runInfo)


if __name__ == "__main__":
    # datafile = "/Users/yongbinfeng/Desktop/CaloX/CaloXDataAnalysis/data/RunlistAugust.csv"
    datafile = "/Users/yongbinfeng/Desktop/runlist_sep.csv"
    data = read_csv(datafile)
    jsonfile = "/Users/yongbinfeng/Desktop/runlist_sep.json"
    write_json(data, jsonfile)
    # testGetRunInfo()
