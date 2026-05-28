"""`python -m risearch_tauso` and the `risearch-tauso` console-script entry point.

Both forward all arguments straight to the bundled RIsearch binary.
"""

from __future__ import annotations

import os
import sys

from . import executable_path


def main() -> "int | None":
    binary = executable_path()
    # execv replaces this Python process with RIsearch — the user sees its exit
    # code and signals directly, no Python wrapper overhead.
    os.execv(binary, [binary, *sys.argv[1:]])


if __name__ == "__main__":
    main()
