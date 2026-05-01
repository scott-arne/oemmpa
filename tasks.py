"""Invoke tasks for documentation workflows."""

import sys

from invoke import task


@task
def docs_build(c):
    """Build Sphinx documentation."""
    c.run(f"{sys.executable} -m sphinx -b html docs docs/_build/html", pty=True)


@task
def docs_check(c):
    """Build Sphinx documentation with warnings treated as errors."""
    c.run(
        f"{sys.executable} -m sphinx -W --keep-going -b html docs docs/_build/html",
        pty=True,
    )


@task
def docs_serve(c, port=8000):
    """Serve the built HTML documentation."""
    c.run(
        f"{sys.executable} -m http.server {int(port)} -d docs/_build/html",
        pty=True,
    )
