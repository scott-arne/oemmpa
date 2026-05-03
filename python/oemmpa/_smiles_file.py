"""SMILES file parsing helpers."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path


VALID_DELIMITERS = {"whitespace", "space", "tab", "comma", "to-eol"}


@dataclass(frozen=True)
class SmilesFileRow:
    """Parsed SMILES file row or row-local parse error."""

    row_number: int
    smiles: str | None = None
    molecule_id: str | None = None
    error: Exception | None = None


def iter_smiles_file(path, delimiter="whitespace", has_header=False):
    """Yield parsed rows from a SMILES file.

    :param path: Plain text or ``.gz`` SMILES file.
    :param delimiter: One of ``"whitespace"``, ``"space"``, ``"tab"``,
        ``"comma"``, or ``"to-eol"``.
    :param has_header: Skip the first physical row when true.
    :returns: Iterator of :class:`SmilesFileRow` objects.
    :raises ValueError: If ``delimiter`` is unsupported.
    """
    delimiter = str(delimiter)
    if delimiter not in VALID_DELIMITERS:
        raise ValueError(f"unsupported SMILES file delimiter: {delimiter}")

    with _open_text(path) as handle:
        for row_number, line in enumerate(handle, start=1):
            if has_header and row_number == 1:
                continue
            if not line.strip():
                continue
            stripped = line.rstrip("\n\r")
            try:
                smiles, molecule_id = _parse_smiles_line(stripped, delimiter)
            except ValueError as exc:
                yield SmilesFileRow(row_number=row_number, error=exc)
            else:
                yield SmilesFileRow(
                    row_number=row_number,
                    smiles=smiles,
                    molecule_id=molecule_id,
                )


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _parse_smiles_line(line, delimiter):
    if delimiter == "whitespace":
        return _split_whitespace(line)
    if delimiter == "space":
        return _split_single_delimiter(line, " ", "space-delimited")
    if delimiter == "tab":
        return _split_single_delimiter(line, "\t", "tab-delimited")
    if delimiter == "comma":
        return _split_single_delimiter(line, ",", "comma-delimited")
    if delimiter == "to-eol":
        return _split_to_eol(line)
    raise ValueError(f"unsupported SMILES file delimiter: {delimiter}")


def _split_whitespace(line):
    parts = line.split(None, 2)
    if not parts:
        raise ValueError("row must contain a SMILES field")
    molecule_id = parts[1] if len(parts) > 1 else None
    return parts[0], molecule_id


def _split_single_delimiter(line, delimiter, label):
    if delimiter not in line:
        raise ValueError(f"row must contain at least two {label} fields")
    smiles, molecule_id, *_extra = line.split(delimiter, 2)
    smiles = smiles.strip()
    molecule_id = molecule_id.strip()
    if not smiles:
        raise ValueError(f"row must contain at least two {label} fields")
    return smiles, molecule_id or None


def _split_to_eol(line):
    for index, character in enumerate(line):
        if character.isspace():
            smiles = line[:index]
            molecule_id = line[index + 1 :]
            if smiles:
                molecule_id = molecule_id.strip(" ")
                return smiles, molecule_id or None
            break
    raise ValueError("row must contain a whitespace to delimit the to-eol fields")
