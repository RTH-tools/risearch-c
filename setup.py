"""
Build hook for the `risearch-tauso` Python package.

The C source for RIsearch1 lives in `RIsearch1/src/` and is compiled by its own
Makefile to `RIsearch1/bin/RIsearch`. This script wraps that build into a normal
setuptools step: at wheel-build time it runs `make`, then copies the binary into
the package as `risearch_tauso/bin/RIsearch` so it ships as package data.

All static metadata is in `pyproject.toml`.
"""

import logging
import os
import shutil
import subprocess

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py
from setuptools.command.egg_info import egg_info

try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:  # `wheel` may be absent for sdist-only builds
    bdist_wheel = None

PACKAGE_NAME = "risearch_tauso"
BINARY_NAME = "RIsearch"


class BinaryDistribution(Distribution):
    """The wheel ships a compiled binary, so it is platform-specific.

    Without this, setuptools tags the wheel `py3-none-any` and pip on a
    Windows or macOS machine would happily install a Linux ELF binary.
    """

    def has_ext_modules(self):
        return True


if bdist_wheel is not None:

    class GenericWheel(bdist_wheel):
        """Tag the wheel `py3-none-<plat>` instead of `cp3X-cp3X-<plat>`.

        The RIsearch binary is invoked as a subprocess and doesn't link against
        Python's ABI, so the wheel is compatible with every CPython 3 on a given
        platform. Without this override, setuptools would pin each wheel to the
        builder's exact Python minor version.
        """

        def get_tag(self):
            _python, _abi, plat = super().get_tag()
            return "py3", "none", plat
else:
    GenericWheel = None

logger = logging.getLogger(__name__)


def _build_risearch_binary(force_target_dir=None):
    """Compile RIsearch1 and copy the binary into the package.

    For `pip install .` / wheel builds, `force_target_dir` is the staged
    `build/lib/risearch` directory. For `pip install -e .`, it is None and the
    binary lands in `src/risearch/bin/` so the editable install can find it.
    """
    root_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(root_dir, "RIsearch1")
    src_dir = os.path.join(base_dir, "src")
    bin_dir = os.path.join(base_dir, "bin")

    if force_target_dir:
        dest_dir = os.path.join(force_target_dir, "bin")
    else:
        dest_dir = os.path.join(root_dir, "src", PACKAGE_NAME, "bin")

    dest_path = os.path.join(dest_dir, BINARY_NAME)

    if not force_target_dir and os.path.exists(dest_path):
        print(f"Binary already exists at {dest_path}, skipping rebuild.")
        return

    print("--- BUILDING RIsearch ---")

    if not os.path.exists(src_dir):
        raise FileNotFoundError(
            f"C source not found at {src_dir}. The sdist should include RIsearch1/."
        )

    os.makedirs(bin_dir, exist_ok=True)
    try:
        subprocess.check_call(["make", "-C", src_dir])
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"make failed to compile RIsearch: {e}")

    compiled_bin = os.path.join(bin_dir, BINARY_NAME)
    if not os.path.exists(compiled_bin):
        raise FileNotFoundError(
            f"make ran successfully, but {compiled_bin} was not produced."
        )

    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy(compiled_bin, dest_path)
    os.chmod(dest_path, 0o755)
    print(f"Installed {dest_path}")


class CustomBuild(build_py):
    """Triggered by `pip install .` and wheel builds."""

    def run(self):
        target_dir = os.path.join(self.build_lib, PACKAGE_NAME)
        _build_risearch_binary(target_dir)
        super().run()


class CustomEggInfo(egg_info):
    """Triggered by both `pip install .` and `pip install -e .`."""

    def run(self):
        _build_risearch_binary(force_target_dir=None)
        super().run()


cmdclass = {
    "build_py": CustomBuild,
    "egg_info": CustomEggInfo,
}
if GenericWheel is not None:
    cmdclass["bdist_wheel"] = GenericWheel

setup(
    distclass=BinaryDistribution,
    cmdclass=cmdclass,
)
