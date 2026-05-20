"""Package version information."""

from importlib.metadata import PackageNotFoundError, version


def get_version() -> str:
    try:
        return version("telomererepeatloci")
    except PackageNotFoundError:
        return "unknown"


__version__ = get_version()
