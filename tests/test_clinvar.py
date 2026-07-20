from plp_rules.clinvar import parse_stars, is_plp
from plp_rules.config import ClinVarParams


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


def test_is_plp_conflict_downweight():
    p1 = ClinVarParams(min_stars=1)
    # Any Benign / Likely_benign co-annotation → not P/LP.
    assert is_plp("Pathogenic,Benign", "criteria_provided,_single_submitter", p1) is False
