from plp_rules.clinvar import parse_stars, is_plp, parse_condition
from plp_rules.config import ClinVarParams


def test_parse_condition_basic():
    # VEP encodes commas as '&'; conditions are '|'-separated.
    assert parse_condition("Breast-ovarian_cancer&_familial_1") == "Breast-ovarian_cancer,_familial_1"
    assert parse_condition("Cardiomyopathy|Long_QT_syndrome") == "Cardiomyopathy; Long_QT_syndrome"
    assert parse_condition("") == ""


def test_parse_condition_drops_placeholders():
    # 'not_provided'/'not_specified' dropped unless they're all there is.
    assert parse_condition("not_provided|Retinitis_pigmentosa") == "Retinitis_pigmentosa"
    assert parse_condition("not_specified") == "not_specified"


def test_parse_stars_ladder():
    assert parse_stars("practice_guideline") == 4
    assert parse_stars("reviewed_by_expert_panel") == 3
    assert parse_stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert parse_stars("criteria_provided,_single_submitter") == 1
    assert parse_stars("no_assertion_criteria_provided") == 0
    assert parse_stars(None) == 0
    assert parse_stars("") == 0


def test_is_plp_requires_min_stars():
    p2 = ClinVarParams(min_stars=2)
    assert is_plp("Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", p2) is True
    assert is_plp("Pathogenic", "criteria_provided,_single_submitter", p2) is False


def test_is_plp_pathogenic_and_likely_pathogenic():
    p1 = ClinVarParams(min_stars=1)
    assert is_plp("Likely_pathogenic", "criteria_provided,_single_submitter", p1) is True
    assert is_plp("Pathogenic/Likely_pathogenic", "criteria_provided,_single_submitter", p1) is True
    assert is_plp("Uncertain_significance", "criteria_provided,_single_submitter", p1) is False


def test_parse_stars_vep_ampersand_encoding():
    # VEP encodes commas as '&' inside CSQ subfields.
    assert parse_stars("criteria_provided&_multiple_submitters&_no_conflicts") == 2
    assert parse_stars("criteria_provided&_single_submitter") == 1


def test_is_plp_with_vep_encoded_revstat():
    p2 = ClinVarParams(min_stars=2)
    # Real-world VEP output form: Pathogenic, 2-star, ampersand-encoded revstat.
    assert is_plp("Pathogenic", "criteria_provided&_multiple_submitters&_no_conflicts", p2) is True
    assert is_plp("Likely_pathogenic", "criteria_provided&_single_submitter", p2) is False  # 1 star < 2


def test_is_plp_conflict_downweight():
    p1 = ClinVarParams(min_stars=1)
    # Any Benign / Likely_benign co-annotation → not P/LP.
    assert is_plp("Pathogenic,Benign", "criteria_provided,_single_submitter", p1) is False
