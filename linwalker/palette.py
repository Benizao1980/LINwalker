"""Colour palettes used across LINwalker."""

from __future__ import annotations

from typing import Dict

DEFAULT_PALETTE: Dict[str, str] = {
    "chicken": "#E6B800",
    "ruminant": "#4F7F6A",
    "pig": "#D081A3",
    "wild bird": "#6A4C93",
    "human": "#1A1A1A",
    "other": "#777777",
}

RESERVOIR_SOURCES = ["chicken", "ruminant", "pig", "wild bird"]


def get_palette(custom: Dict[str, str] | None = None) -> Dict[str, str]:
    """Return a palette dict, optionally updated with custom colours."""
    pal = dict(DEFAULT_PALETTE)
    if custom:
        pal.update(custom)
    return pal
