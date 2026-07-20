"""VEP CSQ header parsing + chunk validation helpers."""
from __future__ import annotations
from dataclasses import dataclass
import re
from typing import Iterable
from .config import ValidateParams


_FORMAT_RE = re.compile(r"Format:\s*([A-Za-z0-9_|\-\.]+)")


@dataclass(frozen=True)
class CSQSchema:
    fields: tuple[str, ...]

    def index_of(self, name: str) -> int:
        return self.fields.index(name)


def parse_csq_format(header_description: str) -> CSQSchema:
    """Parse a VEP CSQ INFO header Description to extract the pipe-separated field list.

    Accepts either the full `##INFO=<...,Description="...">` line or the raw
    Description text. Raises ValueError if `Format: ...` is absent.
    """
    m = _FORMAT_RE.search(header_description)
    if not m:
        raise ValueError("CSQ header Description missing 'Format:' spec")
    fields = tuple(f.strip() for f in m.group(1).split("|") if f.strip())
    if not fields:
        raise ValueError("CSQ Format list is empty")
    return CSQSchema(fields=fields)


def check_required_fields(schema: CSQSchema, required: Iterable[str]) -> list[str]:
    missing = [f for f in required if f not in schema.fields]
    return missing


def symbol_populated_fraction(csq_values: Iterable[str], schema: CSQSchema) -> float:
    """Given an iterable of raw CSQ INFO values (one per variant record), return
    the fraction whose first annotation block has a non-empty SYMBOL subfield."""
    idx = schema.index_of("SYMBOL")
    total = 0
    with_sym = 0
    for raw in csq_values:
        if raw is None or raw == "" or raw == ".":
            total += 1
            continue
        first = raw.split(",", 1)[0]
        parts = first.split("|")
        total += 1
        if idx < len(parts) and parts[idx].strip():
            with_sym += 1
    if total == 0:
        return 0.0
    return with_sym / total


def validate_chunk(
    header_description: str,
    csq_values: Iterable[str],
    params: ValidateParams,
) -> tuple[bool, list[str]]:
    errors: list[str] = []
    try:
        schema = parse_csq_format(header_description)
    except ValueError as e:
        return False, [str(e)]
    missing = check_required_fields(schema, params.required_csq_fields)
    if missing:
        errors.append(f"CSQ missing required subfields: {missing}")
    frac = symbol_populated_fraction(csq_values, schema)
    if frac < params.min_symbol_fraction:
        errors.append(
            f"SYMBOL populated in only {frac:.2%} of records; require ≥ {params.min_symbol_fraction:.0%}"
        )
    return (not errors), errors
