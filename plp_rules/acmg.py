"""ACMG post-processing over InterVar output.

InterVar emits a `InterVar` column like:
  'InterVar: Pathogenic PVS1=1 PS=[...] PM=[...] PP=[...] BA1=0 BS=[...] BP=[...]'
We treat P/LP as the labeled classification, but expose a hook for future
custom rule overrides (e.g. downgrade if BS1 fires in cohort).
"""
from __future__ import annotations
from dataclasses import dataclass
import re


_LABEL_RE = re.compile(r"InterVar:\s*([A-Za-z ]+?)(?:\s+PVS1=|\s+PS=|\s*$)")


@dataclass(frozen=True)
class ACMGCall:
    label: str  # 'Pathogenic' | 'Likely pathogenic' | 'Uncertain significance' | 'Likely benign' | 'Benign' | ''
    is_plp: bool


def parse_intervar_label(intervar_field: str | None) -> ACMGCall:
    if not intervar_field:
        return ACMGCall(label="", is_plp=False)
    m = _LABEL_RE.search(intervar_field)
    label = (m.group(1).strip() if m else "").rstrip(".")
    plp = label.lower() in {"pathogenic", "likely pathogenic"}
    return ACMGCall(label=label, is_plp=plp)


def classify_row(intervar_field: str | None) -> bool:
    """Return whether the InterVar row is ACMG P/LP."""
    return parse_intervar_label(intervar_field).is_plp
