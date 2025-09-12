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


def _print_report(rep: _Report):
    elapsed = time.perf_counter() - rep.t0
    cpu = _cpu_seconds()
    rss = _maxrss_mib()
    parts = [f"wall={elapsed:.3f}s"]
    if cpu is not None:
        parts.append(f"cpu={cpu:.3f}s")
    if rss is not None:
        parts.append(f"rss={rss:.2f}MiB")
    msg = f"[time] {rep.label}: " + "  ".join(parts)
    print(msg, file=sys.stderr)


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
