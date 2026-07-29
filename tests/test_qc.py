from plp_rules.qc import genotype_passes, is_paralog_artifact, cohort_vs_gnomad
from plp_rules.config import QCThresholds


def test_genotype_passes_het():
    th = QCThresholds()
    assert genotype_passes(20, 40, 0.5, is_het=True, is_hom_alt=False, th=th) is True
    assert genotype_passes(20, 40, 0.10, is_het=True, is_hom_alt=False, th=th) is False
    assert genotype_passes(5, 40, 0.5, is_het=True, is_hom_alt=False, th=th) is False
    assert genotype_passes(20, 5, 0.5, is_het=True, is_hom_alt=False, th=th) is False


def test_genotype_passes_hom_alt():
    th = QCThresholds()
    assert genotype_passes(30, 60, 0.98, is_het=False, is_hom_alt=True, th=th) is True
    assert genotype_passes(30, 60, 0.80, is_het=False, is_hom_alt=True, th=th) is False


def test_paralog_artifact_only_when_flag_set():
    assert is_paralog_artifact("PMS2", mapq_flag=True, region_flag=False) is True
    assert is_paralog_artifact("PMS2", mapq_flag=False, region_flag=False) is False
    assert is_paralog_artifact("BRCA1", mapq_flag=True, region_flag=True) is False


def test_cohort_vs_gnomad():
    r = cohort_vs_gnomad(cohort_ac=1, cohort_an=1000, gnomad_popmax_af=1e-4, max_ratio=10.0)
    assert r.passed is True
    r2 = cohort_vs_gnomad(cohort_ac=100, cohort_an=1000, gnomad_popmax_af=1e-4, max_ratio=10.0)
    assert r2.passed is False
    r3 = cohort_vs_gnomad(cohort_ac=1, cohort_an=100000, gnomad_popmax_af=None)
    assert r3.passed is True
    r4 = cohort_vs_gnomad(cohort_ac=100, cohort_an=1000, gnomad_popmax_af=None)
    assert r4.passed is False
