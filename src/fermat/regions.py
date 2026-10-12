"""Searched regions (results/search/regions.json) -> per-m k-limits.

Each region certifies: for every m in its range and every odd k <= kmax(m),
p = k*2^m+1 was tested against every F_n with nmin <= n <= m-2.
The union of such regions is again described by m -> max kmax(m).
"""
from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REGIONS_FILE = os.path.join(ROOT, "results", "search", "regions.json")


def load_regions(path: str = REGIONS_FILE) -> list[dict]:
    with open(path) as f:
        return [r for r in json.load(f) if r.get("completed")]


def region_kmax(r: dict, m: int) -> int:
    if m < r["nmin"] + 2:
        return 0
    if r["type"] == "pbound":           # all p < 2^B
        if m > r["B"] - 1:
            return 0
        return ((1 << r["B"]) - 2) >> m
    if r["type"] == "rect":             # odd k in [1, kmax], mmin <= m <= mmax
        if r["mmin"] <= m <= r["mmax"]:
            return r["kmax"]
        return 0
    raise ValueError(r["type"])


def union_mkmax(regions: list[dict], nmin: int, m_hi: int) -> dict[int, int]:
    """m -> kmax for regions whose nmin <= the requested nmin (so they cover n >= nmin)."""
    out = {}
    for m in range(nmin + 2, m_hi + 1):
        km = max((region_kmax(r, m) for r in regions if r["nmin"] <= nmin), default=0)
        if km >= 1:
            out[m] = km
    return out
