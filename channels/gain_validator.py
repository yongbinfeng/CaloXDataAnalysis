import functools

ALLOWED_GAINS = {"HG", "LG", "Mix"}


def enforce_gain(func):
    """Decorator to enforce that gain is one of HG, LG, Mix (default=HG)."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if "gain" in kwargs:
            gain = kwargs["gain"]
        elif args:
            gain = args[0]
        else:
            gain = "HG"  # default

        if gain not in ALLOWED_GAINS:
            raise ValueError(
                f"gain must be one of {ALLOWED_GAINS}, got {gain}")

        return func(*args, **kwargs)
    return wrapper
