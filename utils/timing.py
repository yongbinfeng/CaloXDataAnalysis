from contextlib import contextmanager
import functools
import atexit
import os
import sys
import time
from dataclasses import dataclass
from typing import Optional
try:
    import resource  # Unix-only
except Exception:
    resource = None


@dataclass
class _Report:
    label: str
    t0: float


def _maxrss_mib() -> Optional[float]:
    if resource is None:
        return None
    r = resource.getrusage(resource.RUSAGE_SELF)
    rss = r.ru_maxrss
    # Linux: KiB, macOS: bytes
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _cpu_seconds() -> Optional[float]:
    if resource is None:
        return None
    u = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (u.ru_utime + u.ru_stime) + (c.ru_utime + c.ru_stime)


# In yongbinfeng/caloxdataanalysis/CaloXDataAnalysis-main/utils/timing.py

_global_mgr = None  # Global registry for the manager


def register_manager(mgr):
    """Register the CaloXAnalysisManager to be reported at exit."""
    global _global_mgr
    _global_mgr = mgr


def _print_report(rep: _Report):
    elapsed = time.perf_counter() - rep.t0
    cpu = _cpu_seconds()
    rss = _maxrss_mib()

    print(f"\n{'='*50}", file=sys.stderr)
    print(f"{'FINAL ANALYSIS REPORT':^50}", file=sys.stderr)
    print(f"{'='*50}", file=sys.stderr)

    # Header: Script/Label Name
    print(f" Target: {rep.label}", file=sys.stderr)
    print(f"{'-'*50}", file=sys.stderr)

    # Resource Metrics aligned with f-strings
    print(f" {'Wall Time:':<15} {elapsed:>10.3f} s", file=sys.stderr)
    if cpu is not None:
        print(f" {'CPU Time:':<15} {cpu:>10.3f} s", file=sys.stderr)
    if rss is not None:
        print(f" {'Max RSS:':<15} {rss:>10.2f} MiB", file=sys.stderr)

    # Global Manager Metrics (Loop count and Cutflow)
    global _global_mgr
    if _global_mgr:
        print(f"{'-'*50}", file=sys.stderr)
        print(f" {'RDF Loops:':<15} {_global_mgr.rdf.GetNRuns():>10}",
              file=sys.stderr)

        # SelectionManager usually handles its own header/formatting
        if hasattr(_global_mgr, 'sel_mgr'):
            _global_mgr.sel_mgr.print_cutflow()

        _global_mgr = None  # Clean up to prevent duplicate reports

    print(f"{'='*50}\n", file=sys.stderr)


def auto_timer(label: Optional[str] = None):
    """
    Start a process-wide timer that prints at program exit.
    Call this once near the top of your script.
    Disable via env: TIMING_OFF=1
    """
    if os.environ.get("TIMING_OFF") == "1":
        return None
    rep = _Report(label or os.path.basename(sys.argv[0]), time.perf_counter())
    atexit.register(_print_report, rep)
    return rep


@contextmanager
def timed_block(label: str = "block"):
    t0 = time.perf_counter()
    try:
        yield
    finally:
        rep = _Report(label, t0)
        _print_report(rep)


def timed(label: Optional[str] = None):
    """@timed() decorator for a function (e.g., your main())."""
    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*a, **k):
            with timed_block(label or fn.__name__):
                return fn(*a, **k)
        return wrapper
    return deco
