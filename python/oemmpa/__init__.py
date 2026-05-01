"""
OEMMPA - Enhanced matched molecular pair capabilities with the OpenEye Toolkits
"""

import os
import re
import warnings

__version__ = "0.1.0"
__version_info__ = (0, 1, 0)


def _ensure_library_compat():
    """Create compatibility symlinks when OpenEye library versions differ from build time.

    When this package is built with shared OpenEye libraries, the compiled extension
    records the exact versioned library filenames (e.g., liboechem-4.3.0.1.dylib).
    If the user upgrades openeye-toolkits, these filenames change and the dynamic
    linker fails to load the extension.

    This function detects version mismatches and creates symlinks from the expected
    (build-time) library names to the actual (runtime) library files.
    """
    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return False

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    pkg_dir = os.path.dirname(__file__)
    created_any = False

    for expected_name in expected_libs:
        if os.path.exists(os.path.join(oe_lib_dir, expected_name)):
            continue

        symlink_path = os.path.join(pkg_dir, expected_name)
        if os.path.islink(symlink_path):
            if os.path.exists(symlink_path):
                continue
            try:
                os.unlink(symlink_path)
            except OSError:
                continue
        elif os.path.exists(symlink_path):
            continue

        match = re.match(r'(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$', expected_name)
        if not match:
            continue
        base_name = match.group(1)

        actual_path = None
        for f in os.listdir(oe_lib_dir):
            if f.startswith(base_name + '-') or f.startswith(base_name + '.'):
                actual_path = os.path.join(oe_lib_dir, f)
                break

        if actual_path:
            try:
                os.symlink(actual_path, os.path.join(pkg_dir, expected_name))
                created_any = True
            except OSError:
                pass

    return created_any


def _preload_shared_libs():
    """Preload OpenEye shared libraries so the C extension can find them.

    On Linux, the extension's RUNPATH (set at build time) normally handles
    dependency resolution, but preloading ensures libraries are available
    even if RUNPATH is stripped (e.g. by certain packaging tools).
    On macOS, @rpath references may not resolve without preloading.

    Only the libraries recorded in ``OPENEYE_EXPECTED_LIBS`` are loaded,
    and they are loaded with ``RTLD_GLOBAL`` so that cross-module C++
    symbol references resolve correctly. Loading the entire OpenEye
    library directory (which can contain 70+ unrelated shared objects)
    would pollute the global symbol namespace and cause segfaults in
    unrelated C extensions such as ``_sqlite3``.
    """
    import ctypes
    import sys
    if sys.platform not in ('linux', 'darwin'):
        return

    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return

    if not os.path.isdir(oe_lib_dir):
        return

    pkg_dir = os.path.dirname(__file__)
    for lib_name in expected_libs:
        # Try the OpenEye lib directory first, then local symlinks
        oe_path = os.path.join(oe_lib_dir, lib_name)
        local_path = os.path.join(pkg_dir, lib_name)
        path = oe_path if os.path.exists(oe_path) else local_path
        if os.path.exists(path) or os.path.islink(path):
            try:
                ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass


def _preload_bundled_libs():
    """Preload libraries bundled by auditwheel from the .libs directory.

    auditwheel repair bundles non-manylinux dependencies (e.g. libraries
    from FetchContent or system packages) into a ``<package>.libs/``
    directory next to the package. The bundled copies have hashed filenames
    and must be loaded before the C extension to satisfy its DT_NEEDED
    entries.

    Libraries may have inter-dependencies, so we do multiple passes
    until no new libraries can be loaded. Libraries are loaded without
    ``RTLD_GLOBAL`` to avoid polluting the global symbol namespace.
    """
    import sys
    if sys.platform != 'linux':
        return

    import ctypes
    pkg_name = __name__
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    site_dir = os.path.dirname(pkg_dir)
    for libs_name in (f'{pkg_name}.libs', f'.{pkg_name}.libs'):
        libs_dir = os.path.join(site_dir, libs_name)
        if not os.path.isdir(libs_dir):
            continue
        remaining = [
            os.path.join(libs_dir, f)
            for f in sorted(os.listdir(libs_dir))
            if '.so' in f
        ]
        while remaining:
            failed = []
            for lib_path in remaining:
                try:
                    ctypes.CDLL(lib_path)
                except OSError:
                    failed.append(lib_path)
            if len(failed) == len(remaining):
                break
            remaining = failed


def _preload_extension_openeye_libs():
    """Preload OpenEye libraries named by the build-tree extension.

    Development builds can produce a shared extension even when generated
    ``_build_info`` reports ``STATIC``. In that case, importing OpenEye first
    masks missing runtime resolution, while a fresh package import fails before
    the extension can load. Read only this extension's linked OpenEye library
    names and load matching files from the OpenEye runtime directory.
    """
    import ctypes
    import importlib.util
    import sys
    from pathlib import Path

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    extension_path = os.path.join(pkg_dir, "_oemmpa.so")
    if not os.path.exists(extension_path):
        return

    try:
        with open(extension_path, "rb") as extension:
            linked = extension.read()
    except OSError:
        return

    pattern = (
        rb"lib(?:oechem|oemath|oesystem|oeplatform|oezstd|zstd)"
        rb"[A-Za-z0-9._+-]*(?:\.dylib|\.so(?:\.\d+)*)"
    )
    lib_names = {
        match.decode("utf-8", errors="ignore")
        for match in re.findall(pattern, linked)
    }
    if not lib_names:
        return

    search_locations = []
    openeye_module = sys.modules.get("openeye")
    openeye_path = getattr(openeye_module, "__path__", None)
    if openeye_path is not None:
        search_locations.extend(openeye_path)

    if not search_locations:
        try:
            openeye_spec = importlib.util.find_spec("openeye")
        except (ImportError, ValueError):
            openeye_spec = None
        if openeye_spec is not None and openeye_spec.submodule_search_locations is not None:
            search_locations.extend(openeye_spec.submodule_search_locations)

    if not search_locations:
        return

    lib_paths = {}
    for package_root in search_locations:
        libs_root = Path(package_root) / "libs"
        if not libs_root.is_dir():
            continue

        # Importing openeye.libs eagerly imports oechem in some environments.
        # The runtime libraries are shipped below openeye/libs, so filesystem
        # discovery preserves the fresh-import condition.
        for root, _, files in os.walk(libs_root):
            for lib_name in files:
                if ".dylib" in lib_name or ".so" in lib_name:
                    lib_paths.setdefault(lib_name, os.path.join(root, lib_name))

    for lib_name in list(lib_names):
        if lib_name.startswith("liboechem-"):
            grid_name = lib_name.replace("liboechem-", "liboegrid-", 1)
            if grid_name not in lib_names and grid_name in lib_paths:
                lib_names.add(grid_name)

    def load_order(lib_name):
        order = (
            "libzstd",
            "liboeplatform",
            "liboesystem",
            "liboegrid",
            "liboemath",
            "liboechem",
        )
        for index, prefix in enumerate(order):
            if lib_name.startswith(prefix):
                return index
        return len(order)

    for lib_name in sorted(lib_names, key=load_order):
        path = lib_paths.get(lib_name)
        if path and os.path.exists(path):
            try:
                ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass


def _check_openeye_version():
    """Check that the OpenEye version matches what was used at build time."""
    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    if not build_version:
        return

    try:
        from openeye import oechem
        runtime_version = oechem.OEToolkitsGetRelease()
        if runtime_version and build_version:
            build_parts = build_version.split('.')[:2]
            runtime_parts = runtime_version.split('.')[:2]
            if build_parts != runtime_parts:
                warnings.warn(
                    f"OpenEye version mismatch: oemmpa was built with "
                    f"OpenEye Toolkits {build_version} but runtime has {runtime_version}. "
                    f"This may cause compatibility issues.",
                    RuntimeWarning
                )
    except ImportError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "Install with: pip install openeye-toolkits",
            ImportWarning
        )


_ensure_library_compat()
_preload_shared_libs()
_preload_bundled_libs()
_preload_extension_openeye_libs()
_check_openeye_version()

from . import _oemmpa  # type: ignore
from . import oemmpa as _swig_proxy

_RAW_BINDING_EXPORTS = (
    "AnalysisMethod",
    "AnalysisStateError",
    "Analyzer",
    "CutBond",
    "CutBondVector",
    "DuplicateIdError",
    "Fragmentation",
    "FragmentationError",
    "FragmentationMethod",
    "FragmentationStrategy",
    "FragmentationVector",
    "Fragmenter",
    "InvalidMoleculeError",
    "InvalidQueryError",
    "LoadError",
    "LoadErrorVector",
    "LoadReport",
    "MatchedPair",
    "MatchedPairVector",
    "MemoryIndex",
    "MissingPropertyError",
    "MoleculeRecord",
    "OEMMPAError",
    "PairScoring",
    "QueryOptions",
    "ScoringMode_FewerCutsThenHeavyAtomChange",
    "ScoringMode_FewerCutsThenHeavyBondChange",
    "ScoringMode_KeepAll",
    "ScoringMode_MinimalHeavyAtomChange",
    "ScoringMode_MinimalHeavyBondChange",
    "ScoringOptions",
    "SmartsFragmentationStrategy",
    "StringVector",
    "Transform",
    "TransformVector",
)

_missing_raw_exports = []
for _name in _RAW_BINDING_EXPORTS:
    if hasattr(_swig_proxy, _name):
        if not hasattr(_oemmpa, _name):
            setattr(_oemmpa, _name, getattr(_swig_proxy, _name))
    else:
        _missing_raw_exports.append(_name)

if _missing_raw_exports:
    raise ImportError(
        "generated oemmpa wrapper is missing raw exports: "
        + ", ".join(_missing_raw_exports)
    )

from .oemmpa import (
    calculate_molecular_weight,
)
from ._facade import Analyzer
from ._loading import LoadReport, RowError
from ._results import (
    PairCollection,
    PairResult,
    TransformCollection,
    TransformResult,
)

del _missing_raw_exports, _name, _swig_proxy

__all__ = [
    "__version__",
    "__version_info__",
    "_oemmpa",
    "Analyzer",
    "LoadReport",
    "PairCollection",
    "PairResult",
    "RowError",
    "TransformCollection",
    "TransformResult",
    "calculate_molecular_weight",
]
