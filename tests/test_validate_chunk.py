import pytest
from plp_rules.csq import parse_csq_format, check_required_fields, symbol_populated_fraction, validate_chunk
from plp_rules.config import ValidateParams


CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|MANE_SELECT|LoF|LoF_filter|'
    'am_pathogenicity|am_class|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|gnomAD_AF|gnomAD_AF_grpmax">'
)


def test_parse_csq_format_extracts_fields():
    s = parse_csq_format(CSQ_HEADER)
    assert "SYMBOL" in s.fields
    assert "am_pathogenicity" in s.fields
    assert s.index_of("Consequence") == 1


def test_parse_csq_format_missing_raises():
    with pytest.raises(ValueError):
        parse_csq_format("no format here")


def test_required_fields_all_present():
    s = parse_csq_format(CSQ_HEADER)
    missing = check_required_fields(s, ValidateParams().required_csq_fields)
    assert missing == []


def test_symbol_populated_fraction():
    s = parse_csq_format(CSQ_HEADER)
    # SYMBOL is index 3.
    csqs = [
        "A|missense_variant|MODERATE|BRCA1|ENSG1|Transcript|ENST1|protein_coding|1|||0.9|likely_pathogenic|Pathogenic|reviewed_by_expert_panel|0.0001|0.0002",
        "A|missense_variant|MODERATE||ENSG2|Transcript|ENST2|protein_coding|1|||0.2|benign|Benign|criteria_provided,_single_submitter|0.01|0.02",
        "A|synonymous_variant|LOW|TP53|ENSG3|Transcript|ENST3|protein_coding|1|||0.1|benign|||0.05|0.06",
    ]
    frac = symbol_populated_fraction(csqs, s)
    assert abs(frac - 2 / 3) < 1e-9


def test_validate_chunk_ok_and_missing_field():
    csqs = [
        "A|missense_variant|MODERATE|BRCA1|ENSG1|Transcript|ENST1|protein_coding|1|||0.9|likely_pathogenic|Pathogenic|reviewed_by_expert_panel|0.0001|0.0002",
    ]
    ok, errs = validate_chunk(CSQ_HEADER, csqs, ValidateParams(min_symbol_fraction=0.5))
    assert ok, errs

    truncated = CSQ_HEADER.replace("|am_pathogenicity", "")
    ok2, errs2 = validate_chunk(truncated, csqs, ValidateParams(min_symbol_fraction=0.0))
    assert not ok2
    assert any("am_pathogenicity" in e for e in errs2)


def test_symbol_fraction_threshold_enforced():
    csqs = [
        "A|missense_variant|MODERATE||ENSG1|Transcript|ENST1|protein_coding|1|||0.9|likely_pathogenic|Pathogenic|reviewed_by_expert_panel|0.0001|0.0002",
        "A|missense_variant|MODERATE||ENSG2|Transcript|ENST2|protein_coding|1|||0.2|benign|Benign|criteria_provided,_single_submitter|0.01|0.02",
    ]
    ok, errs = validate_chunk(CSQ_HEADER, csqs, ValidateParams(min_symbol_fraction=0.5))
    assert not ok
    assert any("SYMBOL" in e for e in errs)
