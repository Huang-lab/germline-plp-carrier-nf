import pytest
from plp_rules.config import parse_classifiers, required_csq_for, VALID_CLASSIFIERS


def test_parse_classifiers_default():
    assert parse_classifiers(None) == VALID_CLASSIFIERS
    assert parse_classifiers("") == VALID_CLASSIFIERS
    assert parse_classifiers([]) == VALID_CLASSIFIERS


def test_parse_classifiers_string():
    assert parse_classifiers("clinvar") == ("clinvar",)
    assert parse_classifiers("clinvar,am") == ("clinvar", "am")
    assert parse_classifiers("acmg, clinvar") == ("clinvar", "acmg")  # stable order


def test_parse_classifiers_list():
    assert parse_classifiers(["am", "clinvar"]) == ("clinvar", "am")


def test_parse_classifiers_bad():
    with pytest.raises(ValueError):
        parse_classifiers("clinvar,foo")


def test_required_csq_for_clinvar_only():
    req = required_csq_for(("clinvar",))
    assert "ClinVar_CLNSIG" in req
    assert "ClinVar_CLNREVSTAT" in req
    assert "am_pathogenicity" not in req
    assert "gnomAD_AF" not in req


def test_required_csq_for_am_only():
    req = required_csq_for(("am",))
    assert "am_pathogenicity" in req
    assert "ClinVar_CLNSIG" not in req


def test_required_csq_union_ordering_stable():
    req = required_csq_for(("clinvar", "acmg", "am"))
    assert req.index("SYMBOL") == 0
    assert "gnomAD_AF" in req
    assert "am_pathogenicity" in req
    assert "ClinVar_CLNSIG" in req
