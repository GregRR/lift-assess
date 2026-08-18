"""Installed package version for liftAssess."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("liftassess")
except PackageNotFoundError:
    # Raw source-tree imports are not the supported execution path. Installed and
    # editable environments obtain the authoritative version from package metadata.
    __version__ = "0+unknown"
