"""Save/load helpers for expensive notebook computations."""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any, Callable

import anndata as ad

from lib.qc_common import ensure_output_dir


def save_result(result: dict[str, Any], path: Path) -> None:
    ensure_output_dir(str(path))
    with Path(path).open("wb") as handle:
        pickle.dump(result, handle)


def load_result(
    path: Path,
    adata: ad.AnnData | None = None,
    compute_fn: Callable[[ad.AnnData], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    path = Path(path)
    if path.exists():
        with path.open("rb") as handle:
            return pickle.load(handle)
    if adata is None or compute_fn is None:
        raise FileNotFoundError(f"Missing cached result at {path}")
    result = compute_fn(adata)
    save_result(result, path)
    return result
