"""Build configuration guardrails."""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tomllib

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]

# Poisoning sys.modules before the package is imported is the only way to show
# that nothing on the p-value path imports SciPy at module scope, and it has to
# happen in a fresh interpreter because pytest has already imported oemmpa by
# the time any test body runs. The package directory is passed in argv so the
# child resolves this worktree rather than whatever editable install is on the
# inherited sys.path.
_SCIPY_FREE_P_VALUE_SCRIPT = """
import json
import sys

sys.modules["scipy"] = None
sys.modules["scipy.stats"] = None
sys.path.insert(0, sys.argv[1])

import oemmpa
from oemmpa import _analytics

print(json.dumps({
    "package": oemmpa.__file__,
    "defined": _analytics._p_value(2.228138851986273, 11, 1.0),
    "undefined": [
        _analytics._p_value(None, 11, 1.0),
        _analytics._p_value(1.0, 1, 1.0),
        _analytics._p_value(1.0, 11, 0.0),
    ],
}))
"""


def _cmake_executable():
    cmake = shutil.which("cmake")
    if cmake is not None:
        return cmake
    return "/opt/homebrew/bin/cmake"


def test_cmake_presets_are_readable_and_expose_debug_build():
    result = subprocess.run(
        [_cmake_executable(), "--list-presets=build"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert '"debug"' in result.stdout


def test_wheel_workflow_uses_installed_package_smoke_tests():
    workflow = (REPO_ROOT / ".github/workflows/build-wheels.yml").read_text()

    assert "pytest tests/python/" not in workflow
    assert "mktemp -d" in workflow
    assert "import oemmpa" in workflow
    assert "python -m oemmpa --help" in workflow


def test_wheel_build_enables_duckdb_storage():
    """Published wheels must be built with DuckDB storage enabled.

    Guards the full chain: pyproject flips the CMake flag for the wheel build,
    each platform job provisions the libduckdb C++ bundle, and every smoke test
    asserts ``duckdb_available()`` so a disabled build fails CI rather than
    silently shipping a wheel that cannot persist.
    """
    pyproject = (REPO_ROOT / "pyproject.toml").read_text()
    workflow = (REPO_ROOT / ".github/workflows/build-wheels.yml").read_text()

    assert 'OEMMPA_BUILD_DUCKDB = "ON"' in pyproject

    # One DuckDB provisioning step per platform job (linux, macos, windows).
    assert workflow.count("Download DuckDB C++ library") == 3
    assert "libduckdb-linux-amd64.zip" in workflow
    assert "libduckdb-osx-universal.zip" in workflow
    assert "libduckdb-windows-amd64.zip" in workflow
    # Windows lacks auditwheel/delocate, so the DLL is vendored explicitly.
    assert "Bundle duckdb.dll into wheel" in workflow
    # Every job's smoke test must fail closed when DuckDB is missing.
    assert "wheel built without DuckDB storage" in workflow


def test_rdkit_license_ships_with_binary_distributions():
    """RDKit's BSD notice has to reach the wheel, not only the source tree.

    ``src/RingPerception.cpp`` carries the notice in its header comment, which
    covers source redistribution (BSD condition 1). Condition 2 covers binary
    forms, and a wheel ships no ``.cpp``: the notice gets there only as a file
    named in ``license-files``, which the build backend copies into the
    dist-info ``licenses/`` directory.
    """
    license_file = REPO_ROOT / "LICENSES" / "RDKit-BSD-3-Clause.txt"

    assert license_file.is_file(), f"missing {license_file}"

    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    assert "LICENSES/RDKit-BSD-3-Clause.txt" in pyproject["project"]["license-files"]

    notice = license_file.read_text()
    assert "Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors" in notice
    assert (
        "Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, "
        "and Julie Penzotti and others" in notice
    )

    # The shipped file must stay byte-for-byte the text carried in the source
    # header, or the two drift apart silently. Both blocks are bounded by their
    # first and last lines, stripped of the comment markers, and required to
    # appear in the license file in order.
    source_lines = _comment_body_lines(REPO_ROOT / "src" / "RingPerception.cpp")
    notice_block = _bounded_block(
        source_lines,
        "Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors",
        "of the RDKit source tree.",
    )
    license_block = _bounded_block(
        source_lines,
        "BSD 3-Clause License",
        "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.",
    )
    shipped = [line.strip() for line in notice.splitlines() if line.strip()]
    for block in (notice_block, license_block):
        assert len(block) > 1, "source header block not found"
        assert _is_contiguous_run(block, shipped), block[0]


def _comment_body_lines(path):
    """Header-comment text of a C++ source with its ``//`` markers removed.

    The RDKit notice is quoted inside the OEMMPA comment as ``//   //  text``,
    so the marker is stripped twice; every line is then whitespace-trimmed.
    """
    lines = []
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped.startswith("//"):
            continue
        body = stripped[2:].strip()
        if body.startswith("//"):
            body = body[2:].strip()
        if body:
            lines.append(body)
    return lines


def _bounded_block(lines, first, last):
    """The slice of ``lines`` from the line equal to ``first`` through ``last``."""
    start = lines.index(first)
    end = lines.index(last, start)
    return lines[start : end + 1]


def _is_contiguous_run(block, haystack):
    """Whether ``block`` occurs as a contiguous run inside ``haystack``."""
    width = len(block)
    return any(haystack[i : i + width] == block for i in range(len(haystack) - width + 1))


#: The install rule the notices depend on, pinned to its exact form: the two
#: file names in either order and nothing else between ``FILES`` and the
#: destination, then the destination and the closing parenthesis. CMake takes
#: newlines as argument separators, so the rule may be written on one line or
#: five without changing meaning, and ``\s+`` covers both. Everything the form
#: excludes is a way to keep the files out of a package while the rule still
#: reads correctly: ``EXCLUDE_FROM_ALL`` leaves the rule out of a full install,
#: ``CONFIGURATIONS`` stages the notices for some build types only, ``OPTIONAL``
#: downgrades a missing notice to silence, and CMake accepts each of them before
#: the destination as readily as after it, which is why the file list is matched
#: exactly rather than searched. ``COMPONENT`` and a third file fail for the
#: same reason, being spellings this rule does not use, and the closing
#: parenthesis right after the destination is what keeps
#: ``share/licenses/oemmpa/extra`` out. ``install`` itself has to open its line:
#: CMake reads one command invocation per line, so a command token can follow
#: nothing but indentation, and that is what turns away ``not_install(...)``
#: and ``message(STATUS install(...))`` -- text that is executable, but not as
#: an install.
#: The pattern says what a rule means, not whether CMake would run it, so it is
#: applied through ``_cmake_stages_license_notices`` -- to text whose comments
#: and whose quoted and bracket arguments have been blanked, and only at the
#: top level of the listfile -- rather than to the raw file. What it has to
#: match is therefore the unquoted, command-level spelling ``CMakeLists.txt``
#: uses: the same characters inside a ``message()`` or a ``set()`` are a string
#: and install nothing.
_INSTALL_LICENSES_RULE = re.compile(
    r"^[ \t]*install\s*\(\s*FILES\s+"
    r"(?:LICENSE\s+LICENSES/RDKit-BSD-3-Clause\.txt"
    r"|LICENSES/RDKit-BSD-3-Clause\.txt\s+LICENSE)"
    r"\s+DESTINATION\s+share/licenses/oemmpa\s*\)",
    re.MULTILINE,
)

#: A CMake bracket delimiter's opening, whether it follows a ``#`` and opens a
#: bracket comment or stands alone and opens a bracket argument. The equals signs
#: are captured because the closing delimiter has to carry the same number of
#: them, which is the entire reason the ``[=[`` spelling exists.
_CMAKE_BRACKET_OPEN = re.compile(r"\[(=*)\[")

#: The opening and closing lines of every CMake block -- conditional, loop,
#: function or macro definition, variable scope -- however the command name is
#: cased or spaced, and the ``return()`` that ends a listfile's processing.
_CMAKE_BLOCK_OPEN = re.compile(
    r"\s*(?:if|foreach|while|function|macro|block)\s*\(", re.IGNORECASE
)
_CMAKE_BLOCK_CLOSE = re.compile(
    r"\s*end(?:if|foreach|while|function|macro|block)\s*\(", re.IGNORECASE
)
_CMAKE_RETURN = re.compile(r"\s*return\s*\(", re.IGNORECASE)


def _blanked(chunk):
    """``chunk`` with every character but its newlines replaced by a space."""
    return re.sub(r"[^\n]", " ", chunk)


def _cmake_bracket_end(text, opening):
    """The index just past the bracket construct that ``opening`` begins.

    :param text: The contents of a CMake listfile.
    :param opening: An ``_CMAKE_BRACKET_OPEN`` match for the opening delimiter.
    :returns: The index after the matching close, or the end of the text when
        the construct is never closed.
    """
    closing = f"]{opening.group(1)}]"
    closed_at = text.find(closing, opening.end())
    return len(text) if closed_at < 0 else closed_at + len(closing)


def _blank_cmake_comments_and_arguments(text):
    """Overwrite everything CMake would not execute as a command with spaces.

    What survives is what CMake runs: command names, parentheses and unquoted
    arguments. Comments go, in both spellings -- ``#`` to the end of the line and
    the ``#[[ ... ]]`` bracket form with any number of equals signs -- and with
    them the contents of quoted and of ``[[ ... ]]`` bracket arguments, which are
    data. A pattern applied to the result therefore cannot be satisfied by a
    listfile that merely mentions a rule in a ``message()`` or a ``set()``.

    Blanking rather than deleting keeps the result line-for-line aligned with
    the input, so a line index taken from the stripped text still names the same
    line of the file. A backslash escape inside a quoted argument is consumed
    whole, two characters for two spaces, so an escaped quote does not end the
    string. The quotes themselves stay, so what is blanked is still recognisably
    one argument of the surrounding command.

    :param text: The contents of a CMake listfile.
    :returns: The same text with everything but the executable commands replaced
        by spaces.
    """
    kept = []
    index = 0
    while index < len(text):
        char = text[index]
        if char == "#":
            bracket = _CMAKE_BRACKET_OPEN.match(text, index + 1)
            if bracket is None:
                end = text.find("\n", index)
                end = len(text) if end < 0 else end
            else:
                end = _cmake_bracket_end(text, bracket)
            kept.append(_blanked(text[index:end]))
            index = end
            continue
        bracket = _CMAKE_BRACKET_OPEN.match(text, index)
        if bracket is not None:
            end = _cmake_bracket_end(text, bracket)
            kept.append(_blanked(text[index:end]))
            index = end
            continue
        if char == '"':
            kept.append(char)
            index += 1
            while index < len(text) and text[index] != '"':
                # Consuming the escape whole is what keeps \" from closing the
                # string, and a # inside it from starting a comment.
                width = 2 if text[index] == "\\" and index + 1 < len(text) else 1
                kept.append(_blanked(text[index : index + width]))
                index += width
            if index < len(text):
                kept.append('"')
                index += 1
            continue
        kept.append(char)
        index += 1
    return "".join(kept)


def _cmake_is_top_level(lines, target_index):
    """Whether the line at ``target_index`` runs unconditionally, once, at load.

    A relative of ``_cmake_duckdb_guard_depth`` below that answers "is this
    line inside anything at all" without asking what. Every block counts --
    ``if()`` and ``while()`` because their condition is not evaluated here,
    ``foreach()`` because a loop body is not a top-level statement,
    ``function()``, ``macro()`` and ``block()`` because their bodies run only
    when something calls or enters them -- and a ``return()`` at the top level
    ends the listfile, so nothing after it runs either. The only caller is
    checking a rule documented as unconditional, where any enclosure is
    already the finding; a ``return()`` inside a block is a conditional early
    exit this pass does not evaluate, which is the staged install's half.
    With comments and argument contents already blanked, every parenthesis
    left is structural, so a block keyword counts only on a line that starts
    at parenthesis depth zero; one written inside another command's argument
    list is that command's data.

    :param lines: The listfile's lines, comments and argument contents already
        blanked.
    :param target_index: Index of the line to report on.
    :returns: True when no block encloses the line and no top-level
        ``return()`` precedes it.
    """
    depth = 0
    open_parentheses = 0
    for line in lines[:target_index]:
        # A block keyword is a command only where it opens a command: at
        # structural parenthesis depth zero. Inside another command's
        # parentheses the same word is an argument, so message( / endfunction()
        # / ) neither closes the function around it nor opens anything. The
        # keyword line is classified before its own parentheses are counted,
        # because the opener's parenthesis belongs to the block it opens.
        if open_parentheses == 0:
            # Command names are case-insensitive to CMake and may be followed
            # by whitespace, so IF (FALSE) guards a rule as surely as if(FALSE).
            if _CMAKE_BLOCK_OPEN.match(line):
                depth += 1
            elif _CMAKE_BLOCK_CLOSE.match(line) and depth:
                depth -= 1
            elif depth == 0 and _CMAKE_RETURN.match(line):
                return False
        open_parentheses = max(0, open_parentheses + line.count("(") - line.count(")"))
    return depth == 0


def _cmake_stages_license_notices(cmake_text):
    """Whether a listfile has a rule that will actually stage both notices.

    Searching the raw text cannot tell a rule from a mention of one: commented
    out, or quoted as the argument of some other command, it matches the same
    characters and installs nothing. So comments and argument contents go first,
    leaving only what CMake would execute, and what survives has to open a
    command of its own -- ``install`` leading its line, at parenthesis depth
    zero, since an unquoted ``install(...)`` written inside another command's
    parentheses is just that command's arguments -- and to sit at the top
    level of the listfile: outside every block, be it a conditional, a loop, a
    function or macro body or a variable scope, and before any top-level
    ``return()``, the rule being documented as unconditional. A rule disabled
    any other way -- an option that is off, a path that does not exist, an
    override further down -- is the staged install's half of the check, not
    this one's.

    :param cmake_text: The contents of a CMake listfile.
    :returns: True when a live, unquoted ``install(FILES ...)`` rule stages
        ``LICENSE`` and the RDKit notice into ``share/licenses/oemmpa``.
    """
    stripped = _blank_cmake_comments_and_arguments(cmake_text)
    lines = stripped.splitlines()
    # Every match, not merely the first: a guarded copy earlier in the file must
    # not mask the live rule below it. With comments and argument contents
    # blanked, every parenthesis left is structural, so a match preceded by as
    # many closing as opening parentheses opens a command of its own rather
    # than sitting, line-leading or not, inside another command's arguments.
    return any(
        stripped[: match.start()].count("(") == stripped[: match.start()].count(")")
        and _cmake_is_top_level(lines, stripped[: match.start()].count("\n"))
        for match in _INSTALL_LICENSES_RULE.finditer(stripped)
    )


def _build_tree_predates_source(build_tree, source):
    """Whether ``build_tree`` was configured before ``source`` last changed.

    ``cmake --install`` replays ``cmake_install.cmake`` exactly as the configure
    step wrote it, so a tree older than the listfile stages the rules the
    project used to have. The caller skips on that rather than failing: a skip
    is visibly not a pass, a developer who has just edited ``CMakeLists.txt``
    reads the reason and reconfigures, and whether CI keeps a build tree around
    at all is CI's decision rather than this test's.

    :param build_tree: The CMake binary directory.
    :param source: The listfile the tree should have been configured from.
    :returns: True when the generated install script is missing or older than
        the source.
    """
    generated = build_tree / "cmake_install.cmake"
    if not generated.is_file():
        return True
    return generated.stat().st_mtime_ns < source.stat().st_mtime_ns


def test_cmake_install_stages_the_license_notices(tmp_path):
    """The CMake install surface has to carry the notices, not just the wheel.

    ``license-files`` only reaches distributions built by the Python backend.
    BSD condition 2 covers every binary form, so a package assembled from
    ``cmake --install`` needs the same two files.

    The check is in two halves because neither one is sufficient. The source
    half reads ``CMakeLists.txt`` through ``_cmake_stages_license_notices``,
    which blanks comments and the contents of quoted and bracket arguments
    before searching and requires the surviving rule to open a command of its
    own at the top level of the listfile, outside every block and before any
    top-level ``return()`` -- so what it accepts is the unquoted, command-level
    spelling the file uses. It catches a rule that was deleted, commented out in
    any of CMake's
    spellings, demoted to a string in some other command's argument, parked
    behind a guard, stripped of one of the two notices, split so they no longer
    travel together, retargeted -- to another directory, or into a subdirectory
    of the right one, where nothing looks for the notices -- or qualified, on
    either side of the destination, by a modifier that keeps it out of an
    install, such as ``EXCLUDE_FROM_ALL``, ``CONFIGURATIONS`` or ``OPTIONAL``.

    The staged half catches what reading the source cannot: a rule that is
    present and unguarded but never runs, because it names a path that does not
    exist or is overridden further down. It is worth no more than the build tree
    it replays, though. ``cmake --install`` re-runs ``cmake_install.cmake`` as
    the configure step last generated it, so a tree older than ``CMakeLists.txt``
    skips instead of reporting a pass on rules the project may no longer have.
    """
    cmake_source = REPO_ROOT / "CMakeLists.txt"
    assert _cmake_stages_license_notices(cmake_source.read_text(encoding="utf-8")), (
        "CMakeLists.txt has no unconditional install(FILES ...) rule staging "
        "LICENSE and LICENSES/RDKit-BSD-3-Clause.txt into share/licenses/oemmpa"
    )

    cmake = shutil.which("cmake")
    if cmake is None:
        pytest.skip("cmake is not on PATH")

    # The build tree is gitignored, so it is absent on a fresh checkout.
    build_tree = REPO_ROOT / "build-release"
    if not (build_tree / "CMakeCache.txt").is_file():
        pytest.skip(f"no configured build tree at {build_tree}")
    if _build_tree_predates_source(build_tree, cmake_source):
        pytest.skip(
            f"{build_tree} predates {cmake_source}: the tree has to be "
            "reconfigured before its install script describes the current rules"
        )

    prefix = tmp_path / "prefix"
    result = subprocess.run(
        [cmake, "--install", str(build_tree), "--prefix", str(prefix)],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"

    staged = prefix / "share" / "licenses" / "oemmpa"
    for installed, original in (
        (staged / "LICENSE", REPO_ROOT / "LICENSE"),
        (
            staged / "RDKit-BSD-3-Clause.txt",
            REPO_ROOT / "LICENSES" / "RDKit-BSD-3-Clause.txt",
        ),
    ):
        assert installed.is_file(), f"missing {installed}"
        assert installed.read_bytes() == original.read_bytes(), f"{installed} differs"


#: The rule exactly as ``CMakeLists.txt`` writes it. Every mutation below is a
#: way the notices could stop being staged while the file still looks right.
_LICENSE_RULE_AS_WRITTEN = (
    "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
    "    DESTINATION share/licenses/oemmpa\n"
    ")\n"
)

#: The same rule folded onto one line and without its newline, so that it can be
#: embedded in the argument of another command.
_LICENSE_RULE_ON_ONE_LINE = (
    "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt"
    " DESTINATION share/licenses/oemmpa)"
)

#: ``(id, CMake text, whether the notices are staged)``.
_LICENSE_RULE_MUTATIONS = (
    ("as-written", _LICENSE_RULE_AS_WRITTEN, True),
    (
        "reordered-on-one-line",
        (
            "install(FILES LICENSES/RDKit-BSD-3-Clause.txt LICENSE"
            " DESTINATION share/licenses/oemmpa)\n"
        ),
        True,
    ),
    (
        "trailing-comment-on-the-closing-line",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa\n"
            ")  # notices\n"
        ),
        True,
    ),
    (
        "hash-inside-a-quoted-argument",
        'set(OEMMPA_NOTE "# not a comment")\n' + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "quoted-mention-before-the-live-rule",
        f'message(STATUS "{_LICENSE_RULE_ON_ONE_LINE}")\n' + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "after-a-closed-if-block",
        "if(OEMMPA_BUILD_TESTS)\n    set(OEMMPA_TEST_NOTE 1)\nendif()\n"
        + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "one-line-form-commented-out",
        (
            "#install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt"
            " DESTINATION share/licenses/oemmpa)\n"
        ),
        False,
    ),
    (
        "every-line-commented-out",
        (
            "#install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "#    DESTINATION share/licenses/oemmpa\n"
            "#)\n"
        ),
        False,
    ),
    (
        "inside-a-bracket-comment",
        "#[[\n" + _LICENSE_RULE_AS_WRITTEN + "]]\n",
        False,
    ),
    (
        "quoted-mention-in-set",
        f'set(EXAMPLE "{_LICENSE_RULE_ON_ONE_LINE}")\n',
        False,
    ),
    (
        "quoted-mention-in-message",
        f'message(STATUS "{_LICENSE_RULE_ON_ONE_LINE}")\n',
        False,
    ),
    (
        "multi-line-quoted-mention",
        'set(EXAMPLE "\n' + _LICENSE_RULE_AS_WRITTEN + '")\n',
        False,
    ),
    (
        "bracket-argument-mention",
        "set(EXAMPLE [[\n" + _LICENSE_RULE_AS_WRITTEN + "]])\n",
        False,
    ),
    (
        "inside-if-false",
        "if(FALSE)\n" + _LICENSE_RULE_AS_WRITTEN + "endif()\n",
        False,
    ),
    (
        "inside-upper-case-if-false",
        "IF(FALSE)\n" + _LICENSE_RULE_AS_WRITTEN + "ENDIF()\n",
        False,
    ),
    (
        "inside-if-false-spaced-before-the-parenthesis",
        "if (FALSE)\n" + _LICENSE_RULE_AS_WRITTEN + "endif ()\n",
        False,
    ),
    (
        "after-a-block-closed-by-upper-case-endif",
        "if(OEMMPA_BUILD_TESTS)\n    set(OEMMPA_TEST_NOTE 1)\nENDIF()\n"
        + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "destination-one-directory-deeper",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa/extra\n"
            ")\n"
        ),
        False,
    ),
    (
        "exclude-from-all",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa EXCLUDE_FROM_ALL\n"
            ")\n"
        ),
        False,
    ),
    (
        "configurations-debug",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa CONFIGURATIONS Debug\n"
            ")\n"
        ),
        False,
    ),
    (
        "optional-modifier",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa OPTIONAL\n"
            ")\n"
        ),
        False,
    ),
    (
        "component-modifier",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa COMPONENT notices\n"
            ")\n"
        ),
        False,
    ),
    (
        "exclude-from-all-before-the-destination",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt EXCLUDE_FROM_ALL\n"
            "    DESTINATION share/licenses/oemmpa\n"
            ")\n"
        ),
        False,
    ),
    (
        "optional-before-the-destination",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt OPTIONAL\n"
            "    DESTINATION share/licenses/oemmpa\n"
            ")\n"
        ),
        False,
    ),
    (
        "third-file-added",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt NOTICE\n"
            "    DESTINATION share/licenses/oemmpa\n"
            ")\n"
        ),
        False,
    ),
    (
        "prefixed-command-name",
        "not_" + _LICENSE_RULE_AS_WRITTEN,
        False,
    ),
    (
        "unquoted-rule-inside-message",
        (
            "message(STATUS install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt"
            " DESTINATION share/licenses/oemmpa))\n"
        ),
        False,
    ),
    (
        "unquoted-rule-inside-set",
        (
            "set(EXAMPLE install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt"
            " DESTINATION share/licenses/oemmpa))\n"
        ),
        False,
    ),
    (
        "unquoted-rule-inside-a-multi-line-message",
        "message(\n" + _LICENSE_RULE_AS_WRITTEN + ")\n",
        False,
    ),
    (
        "indented-rule",
        "    " + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "after-a-nested-parenthesis-that-closed",
        "if(NOT (OEMMPA_BUILD_TESTS AND OEMMPA_BUILD_DUCKDB))\n"
        "    set(OEMMPA_TEST_NOTE 1)\nendif()\n" + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "inside-an-uncalled-function",
        "function(oemmpa_stage_notices)\n" + _LICENSE_RULE_AS_WRITTEN + "endfunction()\n",
        False,
    ),
    (
        "inside-an-uncalled-macro",
        "macro(oemmpa_stage_notices)\n" + _LICENSE_RULE_AS_WRITTEN + "endmacro()\n",
        False,
    ),
    (
        "inside-while-false",
        "while(FALSE)\n" + _LICENSE_RULE_AS_WRITTEN + "endwhile()\n",
        False,
    ),
    (
        "inside-a-foreach-body",
        "foreach(notice IN ITEMS once)\n" + _LICENSE_RULE_AS_WRITTEN + "endforeach()\n",
        False,
    ),
    (
        "inside-a-block-scope",
        "block()\n" + _LICENSE_RULE_AS_WRITTEN + "endblock()\n",
        False,
    ),
    (
        "after-a-top-level-return",
        "return()\n" + _LICENSE_RULE_AS_WRITTEN,
        False,
    ),
    (
        "after-a-function-that-closed",
        "function(oemmpa_note)\n    message(STATUS note)\nendfunction()\n"
        + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "block-closer-forged-by-a-multi-line-argument",
        "function(oemmpa_stage_notices)\n"
        "    message(\n"
        "endfunction()\n"
        "    )\n" + _LICENSE_RULE_AS_WRITTEN + "endfunction()\n",
        False,
    ),
    (
        "block-opener-forged-by-a-multi-line-argument",
        "message(\n"
        "if(FALSE)\n"
        ")\n" + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "return-forged-by-a-multi-line-argument",
        "message(\n"
        "return()\n"
        ")\n" + _LICENSE_RULE_AS_WRITTEN,
        True,
    ),
    (
        "license-dropped",
        (
            "install(FILES LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/licenses/oemmpa\n"
            ")\n"
        ),
        False,
    ),
    (
        "notice-dropped",
        "install(FILES LICENSE\n    DESTINATION share/licenses/oemmpa\n)\n",
        False,
    ),
    (
        "destination-retargeted",
        (
            "install(FILES LICENSE LICENSES/RDKit-BSD-3-Clause.txt\n"
            "    DESTINATION share/doc/oemmpa\n"
            ")\n"
        ),
        False,
    ),
    ("empty-file", "", False),
    (
        "split-into-two-rules",
        (
            "install(FILES LICENSE DESTINATION share/licenses/oemmpa)\n"
            "install(FILES LICENSES/RDKit-BSD-3-Clause.txt"
            " DESTINATION share/licenses/oemmpa)\n"
        ),
        False,
    ),
)


@pytest.mark.parametrize(
    ("cmake_text", "expected"),
    [
        pytest.param(cmake_text, expected, id=identifier)
        for identifier, cmake_text, expected in _LICENSE_RULE_MUTATIONS
    ],
)
def test_license_rule_detection_rejects_a_rule_that_cannot_run(cmake_text, expected):
    """Reading the source has to distinguish a live rule from a disabled one.

    Searching the raw text answers "are these characters present", which is not
    the question: a rule that is commented out, in any of the three ways CMake
    spells a comment, or written -- quoted or not -- inside another command's
    parentheses, or spelled as some other command's name, or parked inside
    ``if(FALSE)``, a ``while(FALSE)`` body, a loop, a function or macro nobody
    calls, a ``block()`` scope, or after a top-level ``return()``, or inside a
    function that only a block closer forged by a multi-line argument appears
    to close, stages nothing and reads as a pass -- and one qualified by a
    modifier runs, but not in the install the notices have to reach.
    The accepted mutations are the rewrites that keep staging both notices --
    CMake takes newlines as argument separators and does not care in which order
    the files are named, and only whitespace may stand between the destination
    and the closing parenthesis -- plus the cases a blanking pass can get wrong
    in the other direction: a real trailing comment, a ``#`` inside a quoted
    argument that is not one, and a quoted mention standing above a rule that is
    live.
    """
    assert _cmake_stages_license_notices(cmake_text) is expected


def test_a_build_tree_older_than_the_cmake_source_is_detected(tmp_path):
    """``cmake --install`` replays the configure step, so its age matters.

    The staged half proves nothing about a ``CMakeLists.txt`` the build tree has
    never seen: ``cmake_install.cmake`` still names whatever the last configure
    generated. Comparing its mtime against the source is how the test tells the
    two apart. Equal mtimes count as current, because a tree configured from
    exactly this text is exactly what the staged half wants.
    """
    source = tmp_path / "CMakeLists.txt"
    source.write_text(_LICENSE_RULE_AS_WRITTEN, encoding="utf-8")
    build_tree = tmp_path / "build-release"
    build_tree.mkdir()
    generated = build_tree / "cmake_install.cmake"
    generated.write_text("# configure output\n", encoding="utf-8")

    source_mtime = 1_700_000_000_000_000_000
    os.utime(source, ns=(source_mtime, source_mtime))

    for offset, expected in (
        (-1_000_000_000, True),
        (1_000_000_000, False),
        (0, False),
    ):
        stamp = source_mtime + offset
        os.utime(generated, ns=(stamp, stamp))
        assert _build_tree_predates_source(build_tree, source) is expected

    generated.unlink()
    assert _build_tree_predates_source(build_tree, source) is True


def test_swig_openeye_grid_dependency_is_explicitly_linked():
    swig_interface = (REPO_ROOT / "swig/oemmpa.i").read_text()
    swig_cmake = (REPO_ROOT / "swig/CMakeLists.txt").read_text()

    uses_grid_bindings = "oegrid.h" in swig_interface or "OEScalarGrid" in swig_interface

    assert not uses_grid_bindings or "OpenEye::OEGrid" in swig_cmake


def _cmake_duckdb_guard_depth(lines, target_index):
    """Count the DuckDB ``if()`` blocks still open at ``target_index``.

    CMake conditionals nest, so a stack of "is this block a DuckDB guard" flags
    is the only reading that stays correct when an unrelated ``if()`` sits in
    between.
    """
    stack = []
    for line in lines[:target_index]:
        stripped = line.strip()
        if stripped.startswith("if("):
            stack.append("DUCKDB" in stripped)
        elif stripped.startswith("endif(") and stack:
            stack.pop()
    return sum(stack)


def _inside_duckdb_preprocessor_guard(lines, target_index):
    """Report whether ``target_index`` sits in an ``#if OEMMPA_HAS_DUCKDB`` block.

    Include guards and the OEMedChem block also start with ``#if``, so the
    stack records which conditional is the DuckDB one instead of counting
    directives.
    """
    stack = []
    for line in lines[:target_index]:
        stripped = line.strip()
        if stripped.startswith("#if"):
            stack.append("OEMMPA_HAS_DUCKDB" in stripped)
        elif stripped.startswith("#endif") and stack:
            stack.pop()
    return any(stack)


def _only_index(lines, wanted):
    matches = [index for index, line in enumerate(lines) if line.strip() == wanted]
    assert len(matches) == 1, f"expected exactly one {wanted!r}, found {matches}"
    return matches[0]


def test_statistics_unit_is_built_and_wrapped_without_duckdb():
    """The Student's t p-value must survive a DuckDB-less build.

    ``python/oemmpa/_analytics.py`` calls the C++ ``two_sided_t_p_value`` and
    has no DuckDB dependency of its own, so all three placements of the
    statistics unit have to sit outside the DuckDB guards: the source in the
    unconditional ``OEMMPA_SOURCES`` list, the umbrella include, and the SWIG
    ``%include``. Each is checked against the guard that would hide it, with
    the DuckDB-only sibling as the positive control that the guard is detected.
    """
    cmake = (REPO_ROOT / "CMakeLists.txt").read_text().splitlines()
    umbrella = (REPO_ROOT / "include/oemmpa/oemmpa.h").read_text().splitlines()
    swig = (REPO_ROOT / "swig/oemmpa.i").read_text().splitlines()

    sources_start = _only_index(cmake, "set(OEMMPA_SOURCES")
    sources_end = next(
        index
        for index, line in enumerate(cmake[sources_start:], sources_start)
        if line.strip() == ")"
    )
    statistics_source = _only_index(cmake, "src/Statistics.cpp")
    assert sources_start < statistics_source < sources_end
    assert _cmake_duckdb_guard_depth(cmake, statistics_source) == 0
    # Positive control: the DuckDB-only source is guarded, so a depth of zero
    # above means the helper found no guard, not that it cannot see one.
    duckdb_source = _only_index(
        cmake, "list(APPEND OEMMPA_SOURCES src/DuckDBStore.cpp)"
    )
    assert _cmake_duckdb_guard_depth(cmake, duckdb_source) == 1

    statistics_include = _only_index(umbrella, '#include "oemmpa/Statistics.h"')
    assert not _inside_duckdb_preprocessor_guard(umbrella, statistics_include)
    duckdb_include = _only_index(umbrella, '#include "oemmpa/DuckDBStore.h"')
    assert _inside_duckdb_preprocessor_guard(umbrella, duckdb_include)

    statistics_wrap = _only_index(
        swig, '%include "../include/oemmpa/Statistics.h"'
    )
    swig_guard = _only_index(swig, "#if OEMMPA_HAS_DUCKDB")
    assert statistics_wrap < swig_guard
    assert not _inside_duckdb_preprocessor_guard(swig, statistics_wrap)


def test_p_value_does_not_require_scipy():
    """The p-value comes from the C++ extension, not SciPy.

    Runs in a subprocess that blocks SciPy before the package is imported, so a
    module-level ``import scipy`` anywhere on the p-value path fails the test
    instead of riding on an interpreter that already imported SciPy for the
    Kramer overlay. Deliberately not in test_statistics_parity.py: that module
    is skipped entirely without DuckDB, which is the build this guarantee is
    for.
    """
    package_root = str(REPO_ROOT / "python")
    result = subprocess.run(
        [sys.executable, "-c", _SCIPY_FREE_P_VALUE_SCRIPT, package_root],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    reported = json.loads(result.stdout.strip().splitlines()[-1])

    assert reported["package"].startswith(package_root), reported["package"]
    assert reported["defined"] == pytest.approx(0.05, abs=1e-9)
    # The undefined cases keep returning None: no pairs, zero variance, no t.
    assert reported["undefined"] == [None, None, None]
