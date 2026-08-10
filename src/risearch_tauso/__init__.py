"""RIsearch1 (tauso fork) packaged for Python.

The C binary is bundled as package data. Use `executable_path()` to locate it
and invoke it via `subprocess`. RIsearch is a standalone command-line tool;
this package does not wrap its CLI semantics.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from importlib.resources import files
from pathlib import Path

__all__ = ["executable_path", "__version__"]

# Read back from the installed distribution, whose version comes from
# pyproject.toml -- the same place CMake reads it for the binary's banner. A
# literal here is a second copy that nothing checks against the first.
try:
    __version__ = _installed_version("risearch-tauso")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.0.0+unknown"

_BINARY_NAME = "RIsearch"


def executable_path() -> str:
    """Absolute path to the bundled RIsearch binary, as a string.

    Raises FileNotFoundError if the binary is missing — that means the wheel
    was built incorrectly or an editable install hasn't compiled it yet.
    """
    p = files(__name__) / "bin" / _BINARY_NAME
    if not Path(str(p)).is_file():
        raise FileNotFoundError(
            f"Bundled RIsearch binary not found at {p}. "
            "If installing in editable mode, ensure `make` succeeded during install."
        )
    return str(p)
