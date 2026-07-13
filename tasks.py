"""Invoke tasks for OEMMPA project management."""

from __future__ import annotations

import os
from pathlib import Path
import shlex
import shutil
import sys

from invoke.tasks import task


PROJECT_ROOT = Path(__file__).parent.absolute()
DOCS_DIR = PROJECT_ROOT / "docs"
BUILD_DIR = DOCS_DIR / "_build"
HTML_DIR = BUILD_DIR / "html"
SPHINXBUILD = f"{sys.executable} -m sphinx.cmd.build"

# Generated artifact locations relative to the project root. Only known build
# outputs are listed; local developer files that happen to be gitignored
# (CMakeUserPresets.json, .venv, .vscode) are deliberately excluded so `clean`
# never removes machine configuration.
_CLEAN_DIRS = (
    "build",
    "build-debug",
    "build-release",
    "dist",
    "docs/_build",
    "docs/_doxygen",
    "docs/cpp-api",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
)
# In-tree artifacts copied/generated into the editable package by the CMake
# build (the compiled extension, generated SWIG wrapper, build info, and the
# bundled OpenEye shared/static libraries).
_CLEAN_PACKAGE_FILES = (
    "python/oemmpa/_oemmpa.so",
    "python/oemmpa/oemmpa.py",
    "python/oemmpa/_build_info.py",
)
_CLEAN_PACKAGE_GLOBS = ("lib*.dylib", "lib*.so", "lib*.a")


@task
def docs(ctx, clean=False):
    """Build Sphinx documentation.

    :param clean: Remove the documentation build directory first.
    """
    os.chdir(DOCS_DIR)

    if clean:
        print("Cleaning build directory...")
        ctx.run("make clean", warn=True)

    print("Building documentation...")
    result = ctx.run(f"make html SPHINXBUILD='{SPHINXBUILD}'", warn=True)

    if result.ok:
        print("\nDocumentation built successfully.")
        print(f"Open: file://{HTML_DIR}/index.html")
    else:
        print("\nDocumentation build failed.", file=sys.stderr)
        sys.exit(1)


@task
def serve_docs(ctx, port=8000, clean=False, watch=False):
    """Build documentation and serve it at localhost.

    :param port: Port to serve on.
    :param clean: Remove the documentation build directory first.
    :param watch: Auto-rebuild and reload on source changes.
    """
    if clean:
        os.chdir(DOCS_DIR)
        print("Cleaning build directory...")
        ctx.run("make clean", warn=True)

    if watch:
        print(f"\nWatching for changes and serving at http://localhost:{port}")
        print("Press Ctrl+C to stop.\n")
        os.chdir(DOCS_DIR)
        ctx.run(
            f"{sys.executable} -m sphinx_autobuild"
            f" . {HTML_DIR}"
            f" --port {int(port)}"
            f" --open-browser"
            f" --re-ignore '/_doxygen/'"
            f" --re-ignore '/cpp-api/'"
            f" --re-ignore '/_build/'"
        )
    else:
        docs(ctx, clean=False)

        print(f"\nServing documentation at http://localhost:{port}")
        print("Press Ctrl+C to stop.\n")

        os.chdir(HTML_DIR)
        ctx.run(f"{sys.executable} -m http.server {int(port)}")


@task
def docs_check(ctx):
    """Build documentation with warnings as errors."""
    os.chdir(DOCS_DIR)

    print("Building documentation with strict checking...")
    result = ctx.run(f"make check SPHINXBUILD='{SPHINXBUILD}'", warn=True)

    if result.ok:
        print("\nDocumentation check passed.")
    else:
        print("\nDocumentation check failed. Fix warnings.", file=sys.stderr)
        sys.exit(1)


@task
def docs_deps(ctx):
    """Install documentation dependencies into the active interpreter via uv."""
    # Use uv and target the same interpreter this task runs under (consistent
    # with SPHINXBUILD) rather than a bare pip install, which would silently
    # install into whichever interpreter happens to be active.
    print(f"Installing documentation dependencies into {sys.executable}...")
    ctx.run(f"uv pip install --python {sys.executable} -r {DOCS_DIR}/requirements.txt")
    print("Done.")


@task
def docs_build(ctx, clean=False):
    """Compatibility alias for :func:`docs`."""
    docs(ctx, clean=clean)


@task
def docs_serve(ctx, port=8000, clean=False, watch=False):
    """Compatibility alias for :func:`serve_docs`."""
    serve_docs(ctx, port=port, clean=clean, watch=watch)


@task
def clean(ctx, pycache=False):
    """Remove generated build, docs, and in-tree package artifacts.

    Removes only known build outputs (CMake build trees, generated docs, the
    compiled extension and bundled OpenEye libraries copied into the editable
    package). Local developer files such as ``CMakeUserPresets.json`` and
    ``.venv`` are never touched.

    :param pycache: Also remove ``__pycache__`` directories under the project.
    """
    removed = []

    for relative in _CLEAN_DIRS:
        target = PROJECT_ROOT / relative
        if target.is_dir():
            shutil.rmtree(target)
            removed.append(relative)

    for relative in _CLEAN_PACKAGE_FILES:
        target = PROJECT_ROOT / relative
        if target.exists():
            target.unlink()
            removed.append(relative)

    package_dir = PROJECT_ROOT / "python" / "oemmpa"
    for pattern in _CLEAN_PACKAGE_GLOBS:
        for target in package_dir.glob(pattern):
            target.unlink()
            removed.append(str(target.relative_to(PROJECT_ROOT)))

    if pycache:
        for cache_dir in PROJECT_ROOT.rglob("__pycache__"):
            shutil.rmtree(cache_dir, ignore_errors=True)
            removed.append(str(cache_dir.relative_to(PROJECT_ROOT)))

    if removed:
        print("Removed:")
        for item in removed:
            print(f"  {item}")
    else:
        print("Nothing to clean.")


# Local fallbacks used only when the corresponding env var is unset, so the
# benchmark task works out-of-the-box here without hardcoding a personal path as
# the sole option.
_BENCHMARK_ENV_DEFAULTS = {
    "OPENEYE_ROOT": "/Users/johnss51/Support/openeye/lib/openeye/toolkits",
    "HTTP_PROXY": "http://proxy-server.bms.com:8080",
    "HTTPS_PROXY": "http://proxy-server.bms.com:8080",
    "NO_PROXY": "s3.amazonaws.com,bms.com,localhost,127.0.0.1,169.254.169.254",
}


@task(
    help={
        "head_to_head": "Run only the flagship three-way head-to-head benchmark.",
        "sizes": "Comma-separated molecule counts for head-to-head.",
        "smiles": "Override SMILES corpus path.",
        "output": "Write benchmark rows to a CSV path.",
        "repeats": "Number of timed repeats.",
    }
)
def benchmark(c, head_to_head=False, sizes=None, smiles=None, output=None, repeats=None):
    """Run the OEMMPA benchmark suite with the environment pre-set."""
    from invoke.exceptions import Exit

    # --sizes and --smiles only have meaning for the head-to-head subcommand
    # (the full suite's per-benchmark datasets are fixed by the §6.1 matrix).
    # Reject them for the full suite instead of building an argv the group CLI
    # does not accept.
    if not head_to_head and (sizes is not None or smiles is not None):
        raise Exit(
            "--sizes/--smiles require --head-to-head (the full suite uses fixed "
            "per-benchmark datasets).",
            code=2,
        )

    env = dict(os.environ)
    for key, fallback in _BENCHMARK_ENV_DEFAULTS.items():
        env.setdefault(key, fallback)
    # Ensure the built worktree extension is importable by the suite and the
    # oemmpa CLI it spawns, and that the env's mmpdb/oemmpa executables resolve.
    python_root = str(PROJECT_ROOT / "python")
    swig_root = str(PROJECT_ROOT / "build-debug" / "swig")
    env["PYTHONPATH"] = os.pathsep.join(
        [python_root, swig_root, env.get("PYTHONPATH", "")]
    )
    # Put the running interpreter's bin dir on PATH so the head-to-head
    # subprocess can find `mmpdb`/`oemmpa` even under a minimal PATH.
    interpreter_bin = str(Path(sys.executable).parent)
    env["PATH"] = os.pathsep.join([interpreter_bin, env.get("PATH", "")])

    script = str(PROJECT_ROOT / "benchmarks" / "benchmark_suite.py")
    argv = [sys.executable, script]
    if head_to_head:
        argv.append("head-to-head")
        if sizes is not None:
            argv += ["--sizes", str(sizes)]
        if smiles is not None:
            argv += ["--smiles", str(smiles)]
    if output is not None:
        argv += ["--output", str(output)]
    if repeats is not None:
        argv += ["--repeats", str(repeats)]

    c.run(shlex.join(argv), env=env, pty=False)


@task(
    help={
        "sizes": "Comma-separated molecule counts for the size sweep.",
        "threads": "Comma-separated worker-thread counts for the parallel sweep.",
        "parquet": "Public parquet to sample corpora from.",
        "parallel_size": "Corpus size for the thread sweep (default: largest size).",
        "json": "Write the record+metadata JSON bundle to this path.",
        "html": "Write the self-contained HTML report to this path.",
        "output": "Write benchmark rows to a CSV path.",
        "repeats": "Warm-timing repeats (auto-reduced at large sizes).",
    }
)
def stage_benchmark(
    c,
    sizes=None,
    threads=None,
    parquet=None,
    parallel_size=None,
    json=None,
    html=None,
    output=None,
    repeats=None,
):
    """Run the per-stage OEMMPA/RDKit/MMPDB benchmark and render the web report."""
    env = dict(os.environ)
    for key, fallback in _BENCHMARK_ENV_DEFAULTS.items():
        env.setdefault(key, fallback)
    python_root = str(PROJECT_ROOT / "python")
    # Prefer the release SWIG build for accurate performance numbers, falling
    # back to the debug build when release has not been built.
    release_swig = PROJECT_ROOT / "build-release" / "swig"
    swig_root = str(
        release_swig if release_swig.is_dir() else PROJECT_ROOT / "build-debug" / "swig"
    )
    env["PYTHONPATH"] = os.pathsep.join(
        [python_root, swig_root, env.get("PYTHONPATH", "")]
    )
    interpreter_bin = str(Path(sys.executable).parent)
    env["PATH"] = os.pathsep.join([interpreter_bin, env.get("PATH", "")])

    script = str(PROJECT_ROOT / "benchmarks" / "benchmark_suite.py")
    argv = [sys.executable, script, "stage-benchmark"]
    for flag, value in (
        ("--sizes", sizes),
        ("--threads", threads),
        ("--parquet", parquet),
        ("--parallel-size", parallel_size),
        ("--json", json),
        ("--html", html),
        ("--output", output),
        ("--repeats", repeats),
    ):
        if value is not None:
            argv += [flag, str(value)]

    c.run(shlex.join(argv), env=env, pty=False)
