"""Load demo marker panel configuration."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

from lib.config import load_yaml

_DEMO_PANELS_PATH = Path(__file__).resolve().parent / "demo_marker_panels.yaml"


@lru_cache(maxsize=1)
def load_demo_marker_panels(path: str | Path | None = None) -> dict[str, Any]:
    cfg = load_yaml(path or _DEMO_PANELS_PATH)
    panels = cfg.get("marker_panels", {})
    if not panels:
        raise ValueError("demo_marker_panels.yaml must define marker_panels")
    return cfg


def demo_panel_genes(cfg: dict[str, Any] | None = None) -> list[str]:
    cfg = cfg or load_demo_marker_panels()
    genes: list[str] = []
    for panel_genes in cfg["marker_panels"].values():
        genes.extend(panel_genes)
    background = cfg.get("background_genes", [])
    genes.extend(background)
    return list(dict.fromkeys(genes))
