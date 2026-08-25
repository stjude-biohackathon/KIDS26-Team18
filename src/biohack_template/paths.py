"""Per-analyst output directory helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from lib.paths import safe_name as _safe_name


@dataclass(frozen=True)
class AnalystPaths:
    analyst: str
    output_root: Path
    figures: Path
    tables: Path
    intermediate: Path


def build_analyst_output_dirs(project_root: Path, analyst: str) -> AnalystPaths:
    """Create isolated output folders for one notebook / team member."""
    safe_analyst = _safe_name(analyst)
    output_root = Path(project_root).resolve() / "result" / safe_analyst
    figures = output_root / "figures"
    tables = output_root / "tables"
    intermediate = output_root / "intermediate"

    for directory in (figures, tables, intermediate):
        directory.mkdir(parents=True, exist_ok=True)

    return AnalystPaths(
        analyst=safe_analyst,
        output_root=output_root,
        figures=figures,
        tables=tables,
        intermediate=intermediate,
    )
