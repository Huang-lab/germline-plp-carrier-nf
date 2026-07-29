"""ACMG/AMP classification over InterVar output.

InterVar (WGLab) implements the 2015 ACMG/AMP rules and reports, per variant, a
final classification plus the full triggered-criteria vector in a single field:

    InterVar: Pathogenic PVS1=1 PS=[0, 0, 0, 0, 0] PM=[1, 0, 0, 0, 0, 0, 0]
              PP=[0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]

This module parses that field comprehensively:
  * the final InterVar label,
  * every triggered criterion (PVS1, PSn, PMn, PPn, BA1, BSn, BPn),
  * a compact human-readable criteria string,
  * counts of pathogenic vs benign evidence,
and derives the P/LP boolean. Array lengths are read from the data (not
hard-coded), so the parser is robust to InterVar version differences.

The default P/LP call follows InterVar's own label. `reclassify()` exposes a
hook for cohort-specific post-processing (e.g. demoting on a strong benign
criterion) without changing the default behavior.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import re
from typing import Optional


# Final label: everything after 'InterVar:' up to the first evidence token.
_LABEL_RE = re.compile(r"InterVar:\s*(.+?)(?:\s+PVS1=|\s+PS=|\s*$)", re.IGNORECASE)
# Scalar criteria (PVS1, BA1) and bracketed arrays (PS/PM/PP/BS/BP).
_SCALAR_RE = re.compile(r"\b(PVS1|BA1)\s*=\s*(\d+)")
_ARRAY_RE = re.compile(r"\b(PS|PM|PP|BS|BP)\s*=\s*\[([^\]]*)\]")

_PLP_LABELS = {"pathogenic", "likely pathogenic"}
_BENIGN_LABELS = {"benign", "likely benign"}
# Evidence-code prefixes that are pathogenic vs benign.
_PATHOGENIC_PREFIXES = ("PVS", "PS", "PM", "PP")
_BENIGN_PREFIXES = ("BA", "BS", "BP")


@dataclass(frozen=True)
class ACMGEvidence:
    label: str                          # e.g. 'Pathogenic'
    criteria: tuple[str, ...] = ()      # triggered codes, e.g. ('PVS1','PM2','PP3')
    is_plp: bool = False

    @property
    def pathogenic_criteria(self) -> tuple[str, ...]:
        return tuple(c for c in self.criteria if c.startswith(_PATHOGENIC_PREFIXES))

    @property
    def benign_criteria(self) -> tuple[str, ...]:
        return tuple(c for c in self.criteria if c.startswith(_BENIGN_PREFIXES))

    @property
    def criteria_str(self) -> str:
        return ";".join(self.criteria)


def _triggered(prefix: str, values: list[int]) -> list[str]:
    """['PM1','PM4', ...] for the 1-valued positions (1-indexed)."""
    return [f"{prefix}{i + 1}" for i, v in enumerate(values) if v == 1]


def _ints(csv_ish: str) -> list[int]:
    out: list[int] = []
    for tok in csv_ish.split(","):
        tok = tok.strip()
        if tok == "":
            continue
        try:
            out.append(int(tok))
        except ValueError:
            out.append(0)
    return out


def parse_intervar_evidence(intervar_field: str | None) -> ACMGEvidence:
    """Parse the full InterVar 'InterVar: <label> PVS1=.. PS=[..] ...' field."""
    if not intervar_field:
        return ACMGEvidence(label="", criteria=(), is_plp=False)

    m = _LABEL_RE.search(intervar_field)
    label = (m.group(1).strip().rstrip(".") if m else "")

    criteria: list[str] = []
    # Scalars: PVS1, BA1
    for name, val in _SCALAR_RE.findall(intervar_field):
        if val == "1":
            criteria.append(name)
    # Arrays: PS/PM/PP/BS/BP
    for prefix, body in _ARRAY_RE.findall(intervar_field):
        criteria.extend(_triggered(prefix, _ints(body)))

    # Stable, readable ordering.
    order = {p: i for i, p in enumerate(["PVS", "PS", "PM", "PP", "BA", "BS", "BP"])}
    def _sort_key(code: str):
        pre = next((p for p in order if code.startswith(p)), "ZZ")
        num = re.sub(r"\D", "", code) or "0"
        return (order.get(pre, 99), int(num))
    criteria.sort(key=_sort_key)

    is_plp = label.lower() in _PLP_LABELS
    return ACMGEvidence(label=label, criteria=tuple(criteria), is_plp=is_plp)


def reclassify(ev: ACMGEvidence, *, demote_on_benign_standalone: bool = True) -> ACMGEvidence:
    """Optional post-processing hook. Default keeps InterVar's label.

    If `demote_on_benign_standalone` and BA1 (stand-alone benign) fired, force
    non-P/LP regardless of the label — a conservative safety net. Returns a new
    ACMGEvidence; does not mutate the input.
    """
    if demote_on_benign_standalone and "BA1" in ev.criteria and ev.is_plp:
        return ACMGEvidence(label=ev.label, criteria=ev.criteria, is_plp=False)
    return ev


# ---- Backward-compatible thin wrappers -------------------------------------

@dataclass(frozen=True)
class ACMGCall:
    label: str
    is_plp: bool


def parse_intervar_label(intervar_field: str | None) -> ACMGCall:
    ev = parse_intervar_evidence(intervar_field)
    return ACMGCall(label=ev.label, is_plp=ev.is_plp)


def classify_row(intervar_field: str | None) -> bool:
    """Return whether the InterVar row is ACMG P/LP (InterVar's own label)."""
    return parse_intervar_evidence(intervar_field).is_plp
