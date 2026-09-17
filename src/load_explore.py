"""Bootstrap helpers for ``00_load_and_explore`` notebooks."""

from __future__ import annotations

import sys
from pathlib import Path


def find_repo_root(start: Path | None = None) -> Path:
    """Walk parents from ``start`` until the repo root (contains ``src/``) is found."""
    path = (start or Path.cwd()).resolve()
    for candidate in (path, *path.parents):
        if (candidate / "src").is_dir() and (candidate / "data").is_dir():
            return candidate
    raise FileNotFoundError(
        f"Could not find repo root (expected src/ and data/) starting from {path}"
    )


def repo_root_from_raw_dir(raw_dir: Path) -> Path:
    """Return repo root given a dataset ``raw_data`` directory."""
    return find_repo_root(raw_dir)


def setup_notebook_paths(
    *,
    repo_root: Path | None = None,
    raw_dir: Path | None = None,
) -> Path:
    """Insert ``src/`` on ``sys.path`` so notebooks can import shared helpers.

    Pass either ``repo_root`` or ``raw_dir`` (``.../data/Dataset_XX/raw_data``).
    Returns the resolved repo root.
    """
    if repo_root is None:
        if raw_dir is None:
            raise ValueError("Pass repo_root or raw_dir")
        repo_root = repo_root_from_raw_dir(raw_dir)
    else:
        repo_root = Path(repo_root).resolve()

    src = repo_root / "src"
    if not src.is_dir():
        raise FileNotFoundError(f"Missing src directory: {src}")

    src_str = str(src)
    if src_str not in sys.path:
        sys.path.insert(0, src_str)
    return repo_root
