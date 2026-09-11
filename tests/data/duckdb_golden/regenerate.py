"""Regenerate the committed DuckDB golden fixtures in this directory.

The golden files are the oracle for
``tests/python/test_duckdb_bulk_equivalence.py::test_bulk_save_matches_golden``.
This script is deliberately a thin wrapper over the same three
``tests/python/_duckdb_dump`` helpers that the test uses, and it takes the
fixture path and the ``with_properties`` argument from that module's
``GOLDEN_FIXTURE`` and ``GOLDEN_WITH_PROPERTIES`` constants, so the generator
and the test cannot drift apart. It never touches the proprietary
``tests/data/dhu_glu_ymin.smi``.

The output is byte-reproducible without pinning a thread count: every dump
query in ``_duckdb_dump`` orders by all of its projected columns, so the file
order is fixed by SQL rather than by the order in which ``Analyzer.analyze()``
happens to finish its work units. Passing ``threads=1`` here would only make
the generator diverge from the call the test makes.

Run it from anywhere; every path is resolved from ``__file__``::

    python tests/data/duckdb_golden/regenerate.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

GOLDEN_DIR = Path(__file__).resolve().parent
REPO_ROOT = GOLDEN_DIR.parents[2]

# The compiled ``oemmpa`` package and the shared dump helpers live outside this
# directory, and this script is not collected by pytest, so neither
# ``tests/python/conftest.py`` nor the project root is on ``sys.path`` yet.
for path in (REPO_ROOT / "python", REPO_ROOT / "tests" / "python"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def regenerate(out_dir: Path | None = None) -> None:
    """Build the store from the public fixture and write one TSV per table.

    :param out_dir: Directory to write the TSV files into. Defaults to the
        directory holding this script, i.e. the committed golden fixtures.
    """
    # Both are resolved through the sys.path entries appended above, which a
    # static checker cannot model; duckdb also ships no stubs.
    import _duckdb_dump as dump  # type: ignore[import-not-found]
    import duckdb  # type: ignore[import-not-found]

    target = GOLDEN_DIR if out_dir is None else Path(out_dir)
    with tempfile.TemporaryDirectory() as scratch:
        database = Path(scratch) / "golden.duckdb"
        # The fixture and the with_properties argument come from
        # _duckdb_dump so the generator cannot drift from the test.
        dump.build_store_from_fixture(
            str(dump.GOLDEN_FIXTURE),
            str(database),
            with_properties=dump.GOLDEN_WITH_PROPERTIES,
        )
        connection = duckdb.connect(str(database), read_only=True)
        try:
            dump.write_golden(dump.dump_all(connection), str(target))
        finally:
            connection.close()


def main() -> None:
    """Regenerate into the directory named on the command line, or in place."""
    if len(sys.argv) > 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} [output-directory]")
    regenerate(Path(sys.argv[1]) if len(sys.argv) == 2 else None)


if __name__ == "__main__":
    main()
