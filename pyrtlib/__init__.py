def __getattr__(name):
    """Pyrtlib package initialization file."""
    if name == "__version__":
        from .version import __version__
        return __version__
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
