"""ClinVar P/LP classification.

Inputs are the raw CLNSIG / CLNREVSTAT strings as VEP `--custom` emits them
(underscore-joined, comma-separated multi-values). Star count is derived from
the CLNREVSTAT text per NCBI's review-status ladder.
"""
from __future__ import annotations
from .config import ClinVarParams


# CLNREVSTAT text → star rating (ClinVar's public ladder).
_REVSTAT_STARS: dict[str, int] = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,  # legacy spelling
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_interpretation_for_the_single_variant": 0,
}


def parse_stars(clnrevstat: str | None) -> int:
    if not clnrevstat:
        return 0
    # VEP encodes commas inside CSQ subfields as '&' (its list separator), so
    # 'criteria_provided,_multiple_submitters,_no_conflicts' arrives as
    # 'criteria_provided&_multiple_submitters&_no_conflicts'. Normalize back.
    key = clnrevstat.strip().lower().replace(" ", "_").replace("&", ",")
    return _REVSTAT_STARS.get(key, 0)


def _clnsig_terms(clnsig: str | None) -> list[str]:
    if not clnsig:
        return []
    # CLNSIG may combine values with ',', '/', or (VEP-encoded) '&'. Split on all.
    parts: list[str] = []
    normalized = clnsig.replace("&", ",").replace("/", ",")
    for sub in normalized.split(","):
        sub = sub.strip()
        if sub:
            parts.append(sub)
    return parts


def is_plp(clnsig: str | None, clnrevstat: str | None, params: ClinVarParams) -> bool:
    stars = parse_stars(clnrevstat)
    if stars < params.min_stars:
        return False
    terms = _clnsig_terms(clnsig)
    if not terms:
        return False
    plp_set = {t.lower() for t in params.plp_terms}
    conflict_markers = {"benign", "likely_benign"}
    hit = False
    for t in terms:
        tl = t.lower()
        if tl in plp_set:
            hit = True
        if tl in conflict_markers:
            return False
    return hit
