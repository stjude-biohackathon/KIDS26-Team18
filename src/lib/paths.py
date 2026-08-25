"""Path helpers."""

from __future__ import annotations

import re


def safe_name(value: str) -> str:
    value = str(value)
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")
    return value if value else "sample"
