"""Pure-Python P/LP classification rules.

Cohort-agnostic; all thresholds injected via `config.Params`.
"""
from .config import Params, load_params, parse_classifiers, required_csq_for, VALID_CLASSIFIERS
from . import clinvar, alphamissense, qc
from .clinvar import parse_condition

# NOTE: ACMG-AMP classification is produced by fastVEP (Huang-lab/fastVEP) via
# acmg_fastvep/, not by a plp_rules module. The former plp_rules.acmg
# (ANNOVAR/InterVar evidence parser) was removed with the ANNOVAR/InterVar path.

__all__ = [
    "Params", "load_params", "parse_classifiers", "required_csq_for",
    "VALID_CLASSIFIERS", "parse_condition", "clinvar", "alphamissense", "qc",
]
