"""AlphaMissense gene-specific calibrated P/LP.

Gene-specific thresholds follow Chen/Pejaver 2026; the calibration table is a
user-supplied TSV. Bare `Moderate` normalizes to `PP3_Moderate`. `PP3_Supporting`
never promotes to Moderate.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Optional
from .config import AlphaMissenseParams


_STRENGTHS = ("PP3_Supporting", "PP3_Moderate", "PP3_Strong", "PP3_VeryStrong")
_STRENGTH_RANK = {s: i for i, s in enumerate(_STRENGTHS)}
_ALIASES = {
    "supporting": "PP3_Supporting",
    "moderate": "PP3_Moderate",
    "strong": "PP3_Strong",
    "verystrong": "PP3_VeryStrong",
    "very_strong": "PP3_VeryStrong",
    "pp3_supporting": "PP3_Supporting",
    "pp3_moderate": "PP3_Moderate",
    "pp3_strong": "PP3_Strong",
    "pp3_verystrong": "PP3_VeryStrong",
}


def normalize_strength(s: str) -> str:
    """Normalize a strength label. Bare 'Moderate' -> 'PP3_Moderate'. Unknown raises."""
    if s is None:
        raise ValueError("min_evidence_strength is required")
    key = str(s).strip().lower().replace(" ", "_")
    if key not in _ALIASES:
        raise ValueError(f"Unknown AlphaMissense evidence strength: {s!r}")
    return _ALIASES[key]


@dataclass(frozen=True)
class GeneCalibration:
    """Per-gene thresholds; higher key first for lookup convenience."""
    pp3_supporting: Optional[float] = None
    pp3_moderate: Optional[float] = None
    pp3_strong: Optional[float] = None
    pp3_verystrong: Optional[float] = None

    def threshold_for(self, min_strength: str) -> Optional[float]:
        norm = normalize_strength(min_strength)
        return {
            "PP3_Supporting": self.pp3_supporting,
            "PP3_Moderate": self.pp3_moderate,
            "PP3_Strong": self.pp3_strong,
            "PP3_VeryStrong": self.pp3_verystrong,
        }[norm]


def load_calibration_tsv(path: str) -> dict[str, GeneCalibration]:
    """Load a gene-specific calibration TSV.

    Expected columns (case-insensitive): gene, and any of
    pp3_supporting, pp3_moderate, pp3_strong, pp3_verystrong. Missing cells => None.
    """
    if not path:
        return {}
    table: dict[str, GeneCalibration] = {}
    with open(path, "r", encoding="utf-8") as fh:
        header: list[str] = []
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not header:
                header = [c.strip().lower() for c in fields]
                if "gene" not in header:
                    raise ValueError(f"AM calibration table missing 'gene' column: {path}")
                continue
            row = dict(zip(header, fields))
            gene = row["gene"].strip()
            if not gene:
                continue

            def _f(name: str) -> Optional[float]:
                v = row.get(name, "").strip()
                if v in ("", "NA", "NaN", "."):
                    return None
                return float(v)

            table[gene] = GeneCalibration(
                pp3_supporting=_f("pp3_supporting"),
                pp3_moderate=_f("pp3_moderate"),
                pp3_strong=_f("pp3_strong"),
                pp3_verystrong=_f("pp3_verystrong"),
            )
    return table


def is_plp(
    gene: str,
    am_score: Optional[float],
    calibration: dict[str, GeneCalibration],
    params: AlphaMissenseParams,
) -> bool:
    """Return True if `am_score` meets/exceeds the gene-specific threshold at
    the configured minimum evidence strength.

    - Bare 'Moderate' has already been normalized in params consumers (see
      `normalize_strength`); this function also normalizes defensively.
    - PP3_Supporting is a genuine evidence tier — it does NOT promote to Moderate.
    - If gene missing from calibration, fall back to `params.default_threshold`
      only at PP3_Supporting; higher requested tiers on an uncalibrated gene
      return False (no evidence for a stronger call).
    """
    if am_score is None:
        return False
    try:
        strength = normalize_strength(params.min_evidence_strength)
    except ValueError:
        raise
    gc = calibration.get(gene)
    if gc is None:
        if strength == "PP3_Supporting":
            return am_score >= params.default_threshold
        return False
    thr = gc.threshold_for(strength)
    if thr is None:
        return False
    return am_score >= thr
