from plp_rules.acmg import (
    parse_intervar_label,
    classify_row,
    parse_intervar_evidence,
    reclassify,
)

PATHOGENIC = ("InterVar: Pathogenic PVS1=1 PS=[0,0,0,0,0] PM=[1,0,0,0,0,0,0] "
              "PP=[1,0,0,0,0] BA1=0 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]")
LIKELY = ("InterVar: Likely pathogenic PVS1=0 PS=[0,0,0,0,0] PM=[0,1,0,0,0,0,0] "
          "PP=[1,0,0,0,0] BA1=0 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]")
BENIGN_BA1 = ("InterVar: Benign PVS1=0 PS=[0,0,0,0,0] PM=[0,0,0,0,0,0,0] "
              "PP=[0,0,0,0,0] BA1=1 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]")


def test_label_backcompat():
    assert parse_intervar_label(PATHOGENIC).label == "Pathogenic"
    assert classify_row(PATHOGENIC) is True
    assert classify_row(LIKELY) is True
    assert classify_row(BENIGN_BA1) is False
    assert classify_row("") is False
    assert classify_row(None) is False


def test_evidence_criteria_parsed():
    ev = parse_intervar_evidence(PATHOGENIC)
    assert ev.label == "Pathogenic"
    assert ev.is_plp is True
    assert ev.criteria == ("PVS1", "PM1", "PP1")
    assert ev.pathogenic_criteria == ("PVS1", "PM1", "PP1")
    assert ev.benign_criteria == ()
    assert ev.criteria_str == "PVS1;PM1;PP1"


def test_evidence_likely_pathogenic():
    ev = parse_intervar_evidence(LIKELY)
    assert ev.criteria == ("PM2", "PP1")
    assert ev.is_plp is True


def test_evidence_ba1_benign():
    ev = parse_intervar_evidence(BENIGN_BA1)
    assert ev.criteria == ("BA1",)
    assert ev.benign_criteria == ("BA1",)
    assert ev.is_plp is False


def test_reclassify_demotes_on_ba1():
    # A (hypothetical) P/LP label that also carries BA1 → demoted when requested.
    weird = ("InterVar: Pathogenic PVS1=1 PS=[0,0,0,0,0] PM=[0,0,0,0,0,0,0] "
             "PP=[0,0,0,0,0] BA1=1 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]")
    ev = parse_intervar_evidence(weird)
    assert ev.is_plp is True
    demoted = reclassify(ev, demote_on_benign_standalone=True)
    assert demoted.is_plp is False
    assert "BA1" in demoted.criteria


def test_variable_array_lengths():
    # Robust to InterVar version differences in array sizes.
    s = "InterVar: Likely pathogenic PVS1=0 PS=[0,0,1] PM=[0,0] PP=[1] BA1=0 BS=[0] BP=[0]"
    ev = parse_intervar_evidence(s)
    assert "PS3" in ev.criteria and "PP1" in ev.criteria
