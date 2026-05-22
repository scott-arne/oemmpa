"""Package import behavior tests."""

import importlib
from pathlib import Path
import shutil
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]


_STUB_MODULE = 'class _Stub:\n    def __init__(self, *args, **kwargs):\n        pass\n    def __call__(self, *args, **kwargs):\n        return None\n\ndef __getattr__(name):\n    return _Stub\n'


def test_package_exposes_oemmpa_command_without_legacy_cli_alias():
    """The package should expose the product command name, not the template alias."""
    expected = 'oemmpa = "oemmpa.cli:main"'
    forbidden = "oemmpa-cli"

    for path in (REPO_ROOT / "pyproject.toml", REPO_ROOT / "python" / "pyproject.toml"):
        text = path.read_text(encoding="utf-8")

        assert expected in text
        assert forbidden not in text


def test_package_exports_notebook_friendly_workflow_names():
    import oemmpa

    assert oemmpa.analyze is oemmpa.analyze_dataframe
    assert oemmpa.open is oemmpa.open_store
    assert oemmpa.Objective("pIC50").higher_is_better is True
    assert oemmpa.Selection(property_name="pIC50").property_name == "pIC50"


def test_import_uses_user_cache_for_broken_openeye_runtime_compat_symlink(
    monkeypatch,
    tmp_path,
):
    """Broken OpenEye runtime symlinks should not mutate oemmpa."""
    package = "oemmpa"
    source_dir = tmp_path / package
    shutil.copytree(
        "python/oemmpa",
        source_dir,
        ignore=shutil.ignore_patterns(
            "__pycache__",
            "_*.so",
            "_*.pyd",
            "_*.dylib",
            "lib*.so",
            "lib*.dylib",
            "lib*.a",
        ),
    )
    expected_name = "liboechem-4.3.0.1.so"
    runtime_name = "liboechem-4.3.0.3.so"

    (source_dir / "_build_info.py").write_text(
        "OPENEYE_LIBRARY_TYPE = 'SHARED'\n"
        f"OPENEYE_EXPECTED_LIBS = [{expected_name!r}]\n"
        "OPENEYE_BUILD_VERSION = '2025.2.1'\n"
    )
    (source_dir / "_oemmpa.py").write_text(_STUB_MODULE)
    (source_dir / "oemmpa.py").write_text(_STUB_MODULE)

    fake_openeye = tmp_path / "openeye"
    fake_libs = fake_openeye / "libs"
    fake_runtime = fake_libs / "python3-linux-x64-g++10.x"
    fake_runtime.mkdir(parents=True)
    (fake_openeye / "__init__.py").write_text("")
    marker = tmp_path / "openeye_imported.txt"
    (fake_libs / "__init__.py").write_text(
        f"from pathlib import Path\nPath({str(marker)!r}).write_text('libs')\n"
    )
    (fake_openeye / "oechem.py").write_text(
        f"from pathlib import Path\nPath({str(marker)!r}).write_text('oechem')\n"
    )
    (fake_runtime / runtime_name).write_text("not a real library")
    (fake_runtime / expected_name).symlink_to(fake_runtime / "missing-liboechem.so")
    cache_home = tmp_path / "cache"

    for module_name in list(sys.modules):
        if module_name == package or module_name.startswith(f"{package}."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)
        if module_name == "openeye" or module_name.startswith("openeye."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)

    monkeypatch.setattr(
        sys,
        "meta_path",
        [
            finder
            for finder in sys.meta_path
            if package not in type(finder).__module__
        ],
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_home))
    importlib.invalidate_caches()

    importlib.import_module(package)

    assert not marker.exists()
    assert "openeye.libs" not in sys.modules
    assert "openeye.oechem" not in sys.modules
    assert not (source_dir / expected_name).exists()
    cached_aliases = list(
        cache_home.glob(f"{package}/openeye-libs/**/{expected_name}")
    )
    assert len(cached_aliases) == 1
    assert cached_aliases[0].is_symlink()
    assert cached_aliases[0].resolve().name == runtime_name


def test_import_creates_cache_alias_for_static_named_openeye_runtime_lib(
    monkeypatch,
    tmp_path,
):
    """OpenEye `.a` build records may still point at runtime shared libraries."""
    package = "oemmpa"
    source_dir = tmp_path / package
    shutil.copytree(
        "python/oemmpa",
        source_dir,
        ignore=shutil.ignore_patterns(
            "__pycache__",
            "_*.so",
            "_*.pyd",
            "_*.dylib",
            "lib*.so",
            "lib*.dylib",
            "lib*.a",
        ),
    )
    expected_name = "liboegrid.a"
    runtime_name = "liboegrid-4.3.0.3.dylib"

    (source_dir / "_build_info.py").write_text(
        "OPENEYE_LIBRARY_TYPE = 'SHARED'\n"
        f"OPENEYE_EXPECTED_LIBS = [{expected_name!r}]\n"
        "OPENEYE_BUILD_VERSION = '2025.2.2'\n"
    )
    (source_dir / "_oemmpa.py").write_text(_STUB_MODULE)
    (source_dir / "oemmpa.py").write_text(_STUB_MODULE)

    fake_openeye = tmp_path / "openeye"
    fake_libs = fake_openeye / "libs"
    fake_runtime = fake_libs / "python3-osx-universal-clang++"
    fake_runtime.mkdir(parents=True)
    (fake_openeye / "__init__.py").write_text("")
    (fake_libs / "__init__.py").write_text("")
    (fake_runtime / runtime_name).write_text("not a real library")
    cache_home = tmp_path / "cache"

    for module_name in list(sys.modules):
        if module_name == package or module_name.startswith(f"{package}."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)
        if module_name == "openeye" or module_name.startswith("openeye."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)

    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_home))
    importlib.invalidate_caches()

    importlib.import_module(package)

    assert not (source_dir / expected_name).exists()
    cached_aliases = list(
        cache_home.glob(f"{package}/openeye-libs/**/{expected_name}")
    )
    assert len(cached_aliases) == 1
    assert cached_aliases[0].is_symlink()
    assert cached_aliases[0].resolve().name == runtime_name
