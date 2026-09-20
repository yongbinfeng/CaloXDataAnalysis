"""Location of the channel tables.

Everything the channel maps read from disk (DRS->FERS CSVs, frozen per-era DRS
maps, trigger delay and trigger map, channel selections such as
testingfibers.json) lives under channels/data/, next to the code that
interprets it, so it is versioned with the maps and found from any working
directory.
"""
import os

TABLE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def table_path(*parts):
    """Absolute path of a file (or directory) under channels/data/."""
    return os.path.join(TABLE_DIR, *parts)
