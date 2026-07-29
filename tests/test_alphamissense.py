import pytest
from plp_rules.alphamissense import (
    normalize_strength,
    is_plp,
    load_calibration_tsv,
    GeneCalibration,
)
from plp_rules.config import AlphaMissenseParams


def test_bare_moderate_normalizes_to_pp3_moderate():
    assert normalize_strength("Moderate") == "PP3_Moderate"
    assert normalize_strength("moderate") == "PP3_Moderate"


def test_pp3_supporting_stays_supporting_no_promotion():
    assert normalize_strength("PP3_Supporting") == "PP3_Supporting"
    assert normalize_strength("Supporting") == "PP3_Supporting"


def test_unknown_strength_raises():
    with pytest.raises(ValueError):
        normalize_strength("moderate-ish")
    with pytest.raises(ValueError):
        normalize_strength("")


def test_is_plp_gene_specific_moderate():
    calib = {
        "BRCA1": GeneCalibration(pp3_supporting=0.5, pp3_moderate=0.85, pp3_strong=0.95),
    }
    p_mod = AlphaMissenseParams(min_evidence_strength="Moderate")
    assert is_plp("BRCA1", 0.86, calib, p_mod) is True
    assert is_plp("BRCA1", 0.80, calib, p_mod) is False


def test_is_plp_supporting_does_not_promote_to_moderate():
    calib = {
        "BRCA1": GeneCalibration(pp3_supporting=0.5, pp3_moderate=0.85),
    }
    # At Supporting level, a score of 0.6 counts even though it wouldn't at Moderate.
    p_sup = AlphaMissenseParams(min_evidence_strength="PP3_Supporting")
    p_mod = AlphaMissenseParams(min_evidence_strength="Moderate")
    assert is_plp("BRCA1", 0.6, calib, p_sup) is True
    assert is_plp("BRCA1", 0.6, calib, p_mod) is False


def test_is_plp_missing_gene_falls_back_only_at_supporting():
    p_sup = AlphaMissenseParams(min_evidence_strength="PP3_Supporting", default_threshold=0.564)
    p_mod = AlphaMissenseParams(min_evidence_strength="Moderate", default_threshold=0.564)
    assert is_plp("UNCAL_GENE", 0.9, {}, p_sup) is True
    assert is_plp("UNCAL_GENE", 0.9, {}, p_mod) is False


def test_is_plp_none_score():
    assert is_plp("BRCA1", None, {}, AlphaMissenseParams()) is False


def test_load_calibration_tsv(tmp_path):
    p = tmp_path / "am_calib.tsv"
    p.write_text(
        "gene\tpp3_supporting\tpp3_moderate\tpp3_strong\tpp3_verystrong\n"
        "BRCA1\t0.50\t0.85\t0.95\tNA\n"
        "TP53\t0.55\t0.80\t.\t\n",
        encoding="utf-8",
    )
    t = load_calibration_tsv(str(p))
    assert t["BRCA1"].pp3_moderate == 0.85
    assert t["BRCA1"].pp3_verystrong is None
    assert t["TP53"].pp3_strong is None
