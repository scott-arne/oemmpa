"""Detection of a second DuckDB runtime already loaded into this process.

The pip ``duckdb`` package does not link the shared library oemmpa links: it
carries its own statically linked copy of DuckDB and exports that copy's C++
symbols. C++ weak symbols -- vtables, typeinfo, inline function bodies -- are
coalesced across images by the dynamic loader whatever namespace they were
linked in, so a process holding both copies gives both of them one set. When
the two are different releases, one release's vtables end up dispatching calls
on the other's objects.

On macOS that is not theoretical: importing pip ``duckdb`` 1.5.6 before an
oemmpa linked against 1.5.4 makes every store write fail inside DuckDB with
``INTERNAL Error: Internal system catalog does not have storage``, a message
that names nothing a caller could act on. Importing oemmpa first is safe, since
its own copy is then the one already loaded. Linux is unaffected, because
CPython loads extension modules with ``RTLD_LOCAL`` there and the pip copy's
symbols never reach the global namespace -- so the check is limited to the
platform that has the problem rather than failing an import that works.
"""

import re
import sys

_RELEASE = re.compile(r"^v?(\d+)\.(\d+)\.(\d+)")

_MESSAGE = (
    "oemmpa links DuckDB {ours}, but the pip duckdb package ({theirs}) was "
    "already imported into this process. Both carry their own copy of the "
    "DuckDB C++ runtime, and the loader gives the two one shared set of "
    "vtables, so the mismatched pair makes every store write fail inside "
    "DuckDB with an error that explains nothing.\n"
    "Either import oemmpa before duckdb, which leaves oemmpa's copy in "
    "charge, or install duckdb=={ours} so the two agree."
)


def release_triple(version):
    """The ``(major, minor, patch)`` of a DuckDB version string.

    Accepts both spellings in play: ``duckdb_library_version()`` answers
    ``v1.5.6`` while the pip package's ``__version__`` is ``1.5.6``. Any
    trailing suffix is ignored, so a development build compares on its release.

    :param version: Version string, or None.
    :returns: The release triple, or None when there is no readable one.
    """
    if not version:
        return None
    match = _RELEASE.match(str(version))
    if match is None:
        return None
    return tuple(int(part) for part in match.groups())


def links_shared_duckdb():
    """Whether this build's DuckDB symbols can be coalesced with another copy's.

    A build that absorbed DuckDB from static archives, and exports only its
    module init symbol, has nothing in the dynamic symbol table to coalesce; a
    build linking the shared library does. The answer is recorded at configure
    time because the alternative is reading the extension's load commands on
    every import.

    :returns: True when the linkage is shared, or when this build predates the
        record and cannot say -- the guard costs an informative error at worst,
        and skipping it on an unknown would be the expensive way to be wrong.
    """
    try:
        from . import _build_info
    except ImportError:
        return True
    return getattr(_build_info, "DUCKDB_LINKAGE", "SHARED") == "SHARED"


def conflict(ours, theirs):
    """The message for a mismatched pair of DuckDB runtimes, if there is one.

    :param ours: Version of the DuckDB oemmpa links.
    :param theirs: Version of the already-imported pip duckdb package.
    :returns: The message to raise, or None when the two agree or either
        version cannot be read -- an unreadable version is not evidence of a
        conflict, and refusing the import on it would be a false alarm.
    """
    our_release, their_release = release_triple(ours), release_triple(theirs)
    if our_release is None or their_release is None or our_release == their_release:
        return None
    return _MESSAGE.format(
        ours=".".join(str(part) for part in our_release),
        theirs=".".join(str(part) for part in their_release),
    )


def check(library_version):
    """Refuse the import when a mismatched DuckDB is already loaded.

    Only the already-imported case is checked: a pip duckdb imported after
    oemmpa finds oemmpa's copy in place and is harmless.

    :param library_version: Callable answering the version of the DuckDB
        oemmpa links, as ``duckdb_library_version()`` does. Taken as a callable
        and consulted last, so a package with nothing to compare against never
        pays for the lookup.
    :raises ImportError: When a different pip duckdb release is already loaded
        on a platform where the two collide.
    """
    if sys.platform != "darwin":
        return
    if not links_shared_duckdb():
        return
    loaded = sys.modules.get("duckdb")
    if loaded is None:
        return
    theirs = getattr(loaded, "__version__", None)
    if theirs is None:
        return
    try:
        ours = library_version()
    except Exception:   # noqa: BLE001 - see below
        # A build that cannot name its DuckDB cannot be checked against one.
        # Broad on purpose: this guard exists to stop an import failing
        # obscurely later, so it must never become the reason one fails now.
        return
    message = conflict(ours, theirs)
    if message is not None:
        raise ImportError(message)
