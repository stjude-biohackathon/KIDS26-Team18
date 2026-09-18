"""Label normalization helpers for cross-method comparison."""

from __future__ import annotations

import re


def _strip_cell_suffix(s: str) -> str:
    if s.endswith(" cells"):
        prefix = s[:-6].rstrip()
        if len(prefix) == 1:
            return f"{prefix} cell"
        return prefix
    if s.endswith(" cell"):
        prefix = s[:-5].rstrip()
        words = prefix.split()
        if len(words) >= 2 or (len(words) == 1 and len(words[0]) > 1):
            return prefix
    return s


def _singularize_word(word: str) -> str:
    if len(word) <= 3:
        return word
    if word.endswith("phages"):
        return word[:-1]
    if word.endswith("cytes"):
        return word[:-1]
    if word.endswith("ies") and len(word) > 4:
        return word[:-3] + "y"
    if word.endswith("s") and not word.endswith("ss"):
        return word[:-1]
    return word


def normalize_cell_label(label: str) -> str:
    """
    Canonical form for comparing labels across methods.

    Handles case, trailing 'cell'/'cells', and simple trailing plurals
    (e.g. Smooth Muscle Cells vs Smooth Muscle Cell, Macrophages vs macrophage).
    """
    s = str(label).strip().lower()
    if s in ("nan", "none", ""):
        return ""
    if s in ("unassigned",):
        return "unassigned"
    if s in ("low quality", "low_quality"):
        return "low_quality"
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s)
    s = _strip_cell_suffix(s)
    words = s.split()
    if words:
        words[-1] = _singularize_word(words[-1])
        s = " ".join(words)
    return s


def format_label_display(label: str, *, canonical: str | None = None) -> str:
    """Consistent display casing for plots and tables."""
    base = canonical if canonical is not None else normalize_cell_label(label)
    if base == "unassigned":
        return "Unassigned"
    if base == "low_quality":
        return "Low quality"
    if not base:
        return "Unassigned"
    return base.title()
