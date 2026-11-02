"""Locations for the pion-vs-DREAM study: the scripts (this folder), the repository, and the data directory with the
event dumps and intermediate results (override with the environment variable PION_DREAM_DATA)."""
import os

STUDY_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.abspath(os.path.join(STUDY_DIR, "..", ".."))
DATA_DIR = os.environ.get("PION_DREAM_DATA", "/data1/yofeng/CaloX/pion_vs_dream_dumps")


def code(name):
    """Path of a script of this study, top level first, then exploratory/ (scripts load each other with exec)."""
    for d in (STUDY_DIR, os.path.join(STUDY_DIR, "exploratory")):
        if os.path.exists(os.path.join(d, name)):
            return os.path.join(d, name)
    raise FileNotFoundError(name)
