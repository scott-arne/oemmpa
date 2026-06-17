"""Invoke tasks for OEMMPA project management."""

from __future__ import annotations

import os
from pathlib import Path
import sys

from invoke.tasks import task


PROJECT_ROOT = Path(__file__).parent.absolute()
DOCS_DIR = PROJECT_ROOT / "docs"
BUILD_DIR = DOCS_DIR / "_build"
HTML_DIR = BUILD_DIR / "html"
SPHINXBUILD = f"{sys.executable} -m sphinx.cmd.build"


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
