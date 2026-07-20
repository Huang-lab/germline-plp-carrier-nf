"""Genotype-level QC, paralog artifact filters, and cohort-vs-gnomAD checks."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from .config import QCThresholds


# Known paralog / pseudogene pairs prone to mismapping in short-read WES.
# Variants in the first gene must be excluded if flagged as paralog artifacts
# by any of the reserved INFO/FILTER markers; the actual exclusion policy is
# codified in `is_paralog_artifact`.
PARALOG_PAIRS: dict[str, str] = {
    "PMS2": "PMS2CL",
    "CHEK2": "CHEK2P2",
    "NBN": "NBNP1",
}


def genotype_passes(dp: Optional[int], gq: Optional[int], ab: Optional[float],
                    is_het: bool, is_hom_alt: bool, th: QCThresholds) -> bool:
    if dp is None or dp < th.min_dp:
        return False
    if gq is None or gq < th.min_gq:
        return False
    if is_het:
        if ab is None or ab < th.het_ab_min or ab > th.het_ab_max:
            return False
    if is_hom_alt:
        if ab is None or ab < th.hom_ab_min:
            return False
    return True


def is_paralog_artifact(gene: str, mapq_flag: bool, region_flag: bool) -> bool:
    """A variant is a paralog artifact when it is in a paralog-prone gene AND
    the aligner-derived MAPQ or region flag is set."""
    if gene not in PARALOG_PAIRS:
        return False
    return bool(mapq_flag or region_flag)


@dataclass(frozen=True)
class FrequencyCheck:
    passed: bool
    cohort_af: float
    gnomad_af: Optional[float]
    ratio: Optional[float]


def cohort_vs_gnomad(
    cohort_ac: int,
    cohort_an: int,
    gnomad_popmax_af: Optional[float],
    max_ratio: float = 10.0,
) -> FrequencyCheck:
    """Flag variants inflated ≥ max_ratio over gnomAD popmax; if gnomAD AF is
    absent, only pass when cohort AF itself is <= 1e-3 (guard against novel-artifact spikes)."""
    if cohort_an <= 0:
        return FrequencyCheck(False, 0.0, gnomad_popmax_af, None)
    caf = cohort_ac / cohort_an
    if gnomad_popmax_af is None:
        return FrequencyCheck(caf <= 1e-3, caf, None, None)
    if gnomad_popmax_af <= 0:
        return FrequencyCheck(caf <= 1e-3, caf, gnomad_popmax_af, None)
    ratio = caf / gnomad_popmax_af
    return FrequencyCheck(ratio <= max_ratio, caf, gnomad_popmax_af, ratio)
