"""Internal helpers for lightweight notebook representations."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from html import escape
from itertools import islice
from typing import Any, Protocol, cast


DEFAULT_PREVIEW_ROWS = 5


class _ToDicts(Protocol):
    def to_dicts(self) -> object:
        """Return row dictionaries for preview rendering."""


class _ToDict(Protocol):
    def to_dict(self) -> object:
        """Return one row dictionary for preview rendering."""


def text_summary(name: str, values: Mapping[str, object]) -> str:
    """Return a compact plain-text summary."""
    parts = [f"{key}={value!r}" for key, value in values.items()]
    return f"{name}({', '.join(parts)})"


def text_collection_summary(name: str, count: int) -> str:
    """Return a compact plain-text collection summary."""
    return f"{name}({int(count)} rows)"


def html_summary_card(
    title: object,
    values: Mapping[str, object],
    *,
    actions: Sequence[object] = (),
) -> str:
    """Return an escaped HTML summary card."""
    rows = "\n".join(
        "<tr>"
        f"<th>{escape(str(key))}</th>"
        f"<td>{escape(str(value))}</td>"
        "</tr>"
        for key, value in values.items()
    )
    actions_html = ""
    if actions:
        action_items = "".join(
            f"<li><code>{escape(str(action))}</code></li>" for action in actions
        )
        actions_html = (
            '<div class="oemmpa-actions">'
            "<strong>Next actions</strong>"
            f"<ul>{action_items}</ul>"
            "</div>"
        )
    return (
        '<div class="oemmpa-card">'
        f"<h3>{escape(str(title))}</h3>"
        f"<table>{rows}</table>"
        f"{actions_html}"
        "</div>"
    )


def html_preview_table(
    rows: object,
    *,
    max_rows: int = DEFAULT_PREVIEW_ROWS,
) -> str:
    """Return an escaped bounded HTML preview table."""
    row_limit = max(0, int(max_rows))
    try:
        preview_rows, total_count, has_more = _preview_rows_from(rows, row_limit)
    except Exception:
        # Display hooks must not affect computation if row serialization fails.
        return '<div class="oemmpa-preview-unavailable">Preview unavailable</div>'
    return _html_preview_from_rows(preview_rows, total_count, has_more)


def _html_preview_from_rows(
    preview_rows: Sequence[Mapping[str, object]],
    total_count: int | None,
    has_more: bool,
) -> str:
    if total_count == 0 or (total_count is None and not preview_rows and not has_more):
        return '<div class="oemmpa-preview-empty">No rows</div>'

    columns = _columns_for(preview_rows)
    header = "".join(f"<th>{escape(str(column))}</th>" for column in columns)
    body = []
    for row in preview_rows:
        cells = "".join(
            f"<td>{escape(_format_cell(row.get(column)))}</td>" for column in columns
        )
        body.append(f"<tr>{cells}</tr>")

    omitted = None if total_count is None else total_count - len(preview_rows)
    omitted_html = _omitted_html(omitted, has_more=has_more)

    return (
        '<div class="oemmpa-preview">'
        f"<table><thead><tr>{header}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
        f"{omitted_html}"
        "</div>"
    )


def html_collection_preview(
    title: object,
    rows: object,
    *,
    max_rows: int = DEFAULT_PREVIEW_ROWS,
) -> str:
    """Return an escaped titled preview for collection-like objects."""
    count = _safe_count(rows)
    row_limit = max(0, int(max_rows))
    try:
        preview_rows, total_count, has_more = _preview_rows_from(rows, row_limit)
    except Exception:
        return _html_collection_shell(
            title,
            count,
            '<div class="oemmpa-preview-unavailable">Preview unavailable</div>',
        )
    if count is None and total_count is not None:
        count = total_count
    return _html_collection_shell(
        title,
        count,
        _html_preview_from_rows(preview_rows, total_count, has_more),
    )


def _html_collection_shell(title: object, count: int | None, content: str) -> str:
    if count is None:
        label = escape(str(title))
    else:
        label = f"{escape(str(title))} ({count} rows)"
    return (
        '<div class="oemmpa-collection">'
        f"<h4>{label}</h4>"
        f"{content}"
        "</div>"
    )


def _preview_rows_from(
    rows: object,
    row_limit: int,
) -> tuple[list[dict[str, object]], int | None, bool]:
    if rows is None:
        return [], 0, False
    if isinstance(rows, Mapping):
        return [_string_keyed_dict(rows)][:row_limit], 1, row_limit == 0
    if _is_row_iterable(rows):
        total_count = _safe_count(rows)
        limit = row_limit if total_count is not None else row_limit + 1
        raw_preview = list(islice(cast(Iterable[object], rows), limit))
        has_more = total_count is None and len(raw_preview) > row_limit
        preview = raw_preview[:row_limit]
        return [_row_to_dict(row) for row in preview], total_count, has_more
    if hasattr(rows, "to_dicts"):
        return _preview_rows_from(cast(_ToDicts, rows).to_dicts(), row_limit)
    return [], 0, False


def _row_to_dict(row: object) -> dict[str, object]:
    if isinstance(row, Mapping):
        return _string_keyed_dict(row)
    if hasattr(row, "to_dict"):
        value = cast(_ToDict, row).to_dict()
        if isinstance(value, Mapping):
            return _string_keyed_dict(value)
    return {"value": row}


def _is_row_iterable(rows: object) -> bool:
    return isinstance(rows, Iterable) and not isinstance(rows, str | bytes | bytearray)


def _string_keyed_dict(row: Mapping[Any, Any]) -> dict[str, object]:
    return {str(key): value for key, value in row.items()}


def _safe_count(rows: object) -> int | None:
    if rows is None:
        return 0
    if isinstance(rows, Mapping):
        return 1
    try:
        return len(rows)  # type: ignore[arg-type]
    except Exception:
        return None


def _columns_for(rows: Sequence[Mapping[str, object]]) -> list[str]:
    columns: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                columns.append(key)
                seen.add(key)
    return columns


def _format_cell(value: object) -> str:
    if value is None:
        return ""
    return str(value)


def _omitted_html(omitted: int | None, *, has_more: bool = False) -> str:
    if omitted == 1:
        return '<div class="oemmpa-preview-more">1 more row</div>'
    if omitted is not None and omitted > 1:
        return f'<div class="oemmpa-preview-more">{omitted} more rows</div>'
    if has_more:
        return '<div class="oemmpa-preview-more">more rows available</div>'
    return ""
