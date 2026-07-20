from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Any, Mapping
import json


@dataclass(frozen=True)
class QCThresholds:
    min_dp: int = 10
    min_gq: int = 20
    het_ab_min: float = 0.20
    het_ab_max: float = 0.80
    hom_ab_min: float = 0.90


@dataclass(frozen=True)
class ClinVarParams:
    min_stars: int = 2  # 0=no assertion, 1=single submitter, 2=multi/expert, 3=reviewed by expert panel, 4=practice guideline
    plp_terms: tuple[str, ...] = ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")


@dataclass(frozen=True)
class AlphaMissenseParams:
    calibration_table: str = ""            # path (Chen/Pejaver 2026 gene-specific table)
    min_evidence_strength: str = "PP3_Moderate"
    default_threshold: float = 0.564       # AM default if a gene is absent from calibration


@dataclass(frozen=True)
class ValidateParams:
    required_csq_fields: tuple[str, ...] = (
        "SYMBOL", "Consequence", "am_pathogenicity", "LoF",
        "ClinVar_CLNSIG", "gnomAD_AF", "gnomAD_AF_grpmax",
    )
    min_symbol_fraction: float = 0.5
    variant_count_ratio_min: float = 0.5  # QC precedes VEP; lenient
    variant_count_ratio_max: float = 1.10


@dataclass(frozen=True)
class Params:
    qc: QCThresholds = field(default_factory=QCThresholds)
    clinvar: ClinVarParams = field(default_factory=ClinVarParams)
    am: AlphaMissenseParams = field(default_factory=AlphaMissenseParams)
    validate: ValidateParams = field(default_factory=ValidateParams)
    gnomad_popmax_field: str = "AF_grpmax"  # gnomAD v4
    cancer_gene_panel: tuple[str, ...] = ()  # reporting/QC only, NEVER a filter

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_params(source: str | Mapping[str, Any] | None) -> Params:
    """Load Params from a JSON string, a path, or a mapping. `None` returns defaults."""
    if source is None:
        return Params()
    if isinstance(source, Mapping):
        data = dict(source)
    else:
        s = str(source).strip()
        if s.startswith("{"):
            data = json.loads(s)
        else:
            with open(s, "r", encoding="utf-8") as fh:
                data = json.load(fh)
    return _from_dict(data)


def _from_dict(d: Mapping[str, Any]) -> Params:
    return Params(
        qc=QCThresholds(**d.get("qc", {})),
        clinvar=ClinVarParams(**{k: (tuple(v) if isinstance(v, list) else v)
                                 for k, v in d.get("clinvar", {}).items()}),
        am=AlphaMissenseParams(**d.get("am", {})),
        validate=ValidateParams(**{k: (tuple(v) if isinstance(v, list) else v)
                                   for k, v in d.get("validate", {}).items()}),
        gnomad_popmax_field=d.get("gnomad_popmax_field", "AF_grpmax"),
        cancer_gene_panel=tuple(d.get("cancer_gene_panel", ())),
    )
