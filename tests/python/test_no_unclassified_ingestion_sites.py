"""Fail if a new caller-molecule ingestion site is added without classification.

Scans the C++ sources and the Python direct-call surface for every
``MoleculeRecord::From*`` call (M1), every ``ApplySmirks`` /
``ApplyVariableTransform`` / ``GenerateProducts`` / ``ComputeQueryEnvironments`` /
``AddMolecule`` / ``AddMoleculesFromSmilesFile`` declaration (M2/M3), and every
``_oemmpa.MoleculeRecord.From*`` Python call, and asserts each occurrence is in
the explicit allowlist below classified ``ingestion`` (desalter threaded) or
``internal`` (raw by design). Adding a new site fails this test until it is
classified.

This is the durable fix for the recurring "enumeration missed a sibling caller
path" defect class: a molecule enters oemmpa via THREE C++ mechanisms plus the
Python query surface, and a new object-form overload that bypasses
``MoleculeRecord`` is exactly the kind of site that silently skips desalting.
See ``docs/hyperpowers/specs/2026-07-05-salt-remover-design.md`` section 5.
"""

import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Each key is (relative_path, matched_text); each value is the classification.
# "ingestion" => the site threads the shared Desalter; "internal" => raw by
# design; "by-caller" => the same helper is ingestion or internal depending on
# its caller (resolved at the call site, see the Python query wrappers). The
# scanners below emit a canonical key for every occurrence; a mismatch (new,
# removed, or renamed site) fails the test until the allowlist is updated.
ALLOWLIST = {
    # M1: MoleculeRecord::From* / .From* call sites (definitions excluded).
    ("src/Analyzer.cpp", "MoleculeRecord::FromSmiles"): "ingestion",
    ("src/Analyzer.cpp", "MoleculeRecord::FromMol"): "ingestion",
    ("src/DuckDBStore.cpp", "MoleculeRecord::FromSmiles"): "ingestion",
    ("src/TransformApplication.cpp", "MoleculeRecord::FromSmiles"): "ingestion",
    ("src/QueryEnvironment.cpp", "MoleculeRecord::FromSmiles"): "ingestion",
    ("src/MemoryIndex.cpp", "MoleculeRecord::FromSmiles"): "internal",
    ("python/oemmpa/_rule_environment.py", "_oemmpa.MoleculeRecord.FromSmiles"): "by-caller",
    # M2/M3: caller-molecule operation entry-point OVERLOADS declared in
    # headers, keyed PER OVERLOAD by (symbol, first-param-type) so the
    # std::string source form and the const OEMolBase& object form are DISTINCT
    # allowlist entries. The object (OEMolBase) overloads are the bypass class
    # this guard exists to catch — each must be classified on its own.
    ("include/oemmpa/TransformApplication.h", "ApplySmirks", "std::string"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "ApplySmirks", "OEMolBase"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "ApplyVariableTransform", "std::string"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "ApplyVariableTransform", "OEMolBase"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "GenerateProducts", "std::string"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "GenerateProducts", "OEMolBase"): "ingestion",
    ("include/oemmpa/TransformApplication.h", "ApplyPairTransform", "MatchedPair"): "internal",
    ("include/oemmpa/QueryEnvironment.h", "ComputeQueryEnvironments", "std::string"): "ingestion",
    ("include/oemmpa/Analyzer.h", "AddMolecule", "std::string"): "ingestion",
    ("include/oemmpa/Analyzer.h", "AddMolecule", "OEMolBase"): "ingestion",
    ("include/oemmpa/DuckDBStore.h", "AddMoleculesFromSmilesFile", "std::string"): "ingestion",
}

# Every generate call site that takes a user --source / query molecule MUST
# pass a `desalter=` argument. Keyed by (file, class_or_None, function) so the
# scanner can disambiguate same-named methods on different classes: _query.py
# has BOTH ObjectiveAnalysis.generate (a delegator) and AnalysisResult.generate
# (the real direct generate_products call). A NEW generate path — a fourth
# _generate_* in cli.py OR a new source method on AnalysisResult — that forgets
# `desalter=` fails the test. Delegating wrappers (ObjectiveAnalysis.generate,
# AnalysisResult.opportunities) are covered transitively and are NOT listed.
DESALTER_CALL_SITES = {
    ("python/oemmpa/cli.py", None, "_generate_no_properties"): "generate_products",
    ("python/oemmpa/cli.py", None, "_generate_stateless"): "generate_products",
    ("python/oemmpa/cli.py", None, "_generate_persisted"): "generate_products_from_rule_environments",
    ("python/oemmpa/_query.py", "AnalysisResult", "generate"): "generate_products",
}

# Operation symbols scanned per header. Overloads are keyed by first-parameter
# type; declarations whose first parameter carries no molecule (e.g.
# BuildVariableTransformSmirks) are skipped by _header_overloads.
_OPERATION_SYMBOLS_BY_HEADER = {
    "include/oemmpa/TransformApplication.h": [
        "ApplySmirks", "ApplyVariableTransform", "GenerateProducts", "ApplyPairTransform",
    ],
    "include/oemmpa/QueryEnvironment.h": ["ComputeQueryEnvironments"],
    "include/oemmpa/Analyzer.h": ["AddMolecule"],
    "include/oemmpa/DuckDBStore.h": ["AddMoleculesFromSmilesFile"],
}


def _header_overloads():
    """Return ``{(rel_path, symbol, first_param_type)}`` for each declared overload.

    ``first_param_type`` is ``'std::string'``, ``'OEMolBase'``, or
    ``'MatchedPair'`` — the classifying feature that distinguishes a
    source-SMILES form from the object form that bypasses ``MoleculeRecord``.
    """
    found = set()
    for rel, symbols in _OPERATION_SYMBOLS_BY_HEADER.items():
        text = (REPO / rel).read_text(encoding="utf-8")
        for symbol in symbols:
            # Match each declaration from the symbol up to the first ',' or ')'
            # so only the first parameter is inspected.
            for match in re.finditer(re.escape(symbol) + r"\s*\(([^,)]*)", text):
                first_param = match.group(1)
                if "OEMolBase" in first_param:
                    param_type = "OEMolBase"
                elif "MatchedPair" in first_param:
                    param_type = "MatchedPair"
                elif "std::string" in first_param:
                    param_type = "std::string"
                else:
                    continue  # not a molecule-bearing overload
                found.add((rel, symbol, param_type))
    return found


def _from_call_occurrences():
    """Return ``{(rel, matched)}`` for ``MoleculeRecord::From*`` / ``.From*``.

    Factory definitions (``MoleculeRecord.cpp``/``.h``/``.py``) are excluded so
    only caller sites are classified.
    """
    found = set()
    pattern = re.compile(r"MoleculeRecord(?:::|\.)(FromSmiles|FromMol)\s*\(")
    for sub in ("src", "python/oemmpa"):
        for path in (REPO / sub).rglob("*"):
            if path.suffix not in (".cpp", ".py") or path.name.startswith("MoleculeRecord"):
                continue
            rel = str(path.relative_to(REPO))
            prefix = "_oemmpa.MoleculeRecord." if rel.endswith(".py") else "MoleculeRecord::"
            for match in pattern.finditer(path.read_text(encoding="utf-8")):
                found.add((rel, prefix + match.group(1)))
    return found


def test_every_from_site_is_classified():
    occurrences = _from_call_occurrences()
    unclassified = sorted(occ for occ in occurrences if occ not in ALLOWLIST)
    missing = sorted(key for key in ALLOWLIST if len(key) == 2 and key not in occurrences)
    assert not missing, (
        "ALLOWLIST From* entries no longer present in source (renamed/removed?):\n"
        + "\n".join(map(str, missing))
    )
    assert not unclassified, (
        "Unclassified MoleculeRecord::From* site(s) - classify in ALLOWLIST "
        "('ingestion' => thread the Desalter; 'internal'/'by-caller' => raw):\n"
        + "\n".join(f"{path}: {matched}" for path, matched in unclassified)
    )


def test_every_operation_overload_is_classified():
    overloads = _header_overloads()
    unclassified = sorted(overload for overload in overloads if overload not in ALLOWLIST)
    missing = sorted(key for key in ALLOWLIST if len(key) == 3 and key not in overloads)
    assert not missing, (
        "ALLOWLIST overload entries no longer present in headers (renamed/removed?):\n"
        + "\n".join(map(str, missing))
    )
    assert not unclassified, (
        "Unclassified caller-molecule operation OVERLOAD(s) - each std::string "
        "AND OEMolBase form must be classified 'ingestion'/'internal':\n"
        + "\n".join(f"{path}: {symbol}({param})" for path, symbol, param in unclassified)
    )


def _function_source(rel, class_name, func_name):
    """Return the source text of a specific ``(class, function)`` via AST.

    ``class_name`` is ``None`` for a top-level function. Disambiguates
    same-named methods on different classes (e.g. two ``generate`` methods in
    ``_query.py``).
    """
    text = (REPO / rel).read_text(encoding="utf-8")
    tree = ast.parse(text)
    if class_name is None:
        nodes = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == func_name
        ]
    else:
        cls = next(
            (
                node for node in tree.body
                if isinstance(node, ast.ClassDef) and node.name == class_name
            ),
            None,
        )
        assert cls is not None, f"{rel}: class {class_name} not found; update DESALTER_CALL_SITES"
        nodes = [
            node for node in cls.body
            if isinstance(node, ast.FunctionDef) and node.name == func_name
        ]
    label = f"{rel}: {class_name or '<module>'}.{func_name}"
    assert nodes, f"{label} not found; update DESALTER_CALL_SITES"
    return ast.get_source_segment(text, nodes[0])


def test_generate_source_call_sites_pass_desalter():
    missing = []
    for (rel, class_name, func_name), call in DESALTER_CALL_SITES.items():
        body = _function_source(rel, class_name, func_name)
        label = f"{rel}: {class_name or '<module>'}.{func_name}"
        call_index = body.find(call + "(")
        assert call_index != -1, f"{label} no longer calls {call}(); update DESALTER_CALL_SITES"
        # Scan a generous window past the call opening: these are multi-line
        # calls with nested parens, so a first-close-paren scan would truncate.
        call_region = body[call_index: call_index + 500]
        if "desalter=" not in call_region:
            missing.append(f"{label}: {call}(...) is missing a desalter= argument")
    assert not missing, (
        "Source-ingesting path(s) generate from a user molecule without threading "
        "the desalter:\n" + "\n".join(missing)
    )
