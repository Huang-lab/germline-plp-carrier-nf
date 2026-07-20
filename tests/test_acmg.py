from plp_rules.acmg import parse_intervar_label, classify_row


def test_parse_intervar_pathogenic():
    s = "InterVar: Pathogenic PVS1=1 PS=[0,0,0,0,0] PM=[1,0,0,0,0,0,0] PP=[0,0,0,0,0] BA1=0 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]"
    c = parse_intervar_label(s)
    assert c.label == "Pathogenic"
    assert c.is_plp is True


def test_parse_intervar_likely_pathogenic():
    s = "InterVar: Likely pathogenic PVS1=0 PS=[0,0,0,0,0] PM=[0,1,0,0,0,0,0] PP=[1,0,0,0,0]"
    assert classify_row(s) is True


def test_parse_intervar_vus_or_benign():
    for s in [
        "InterVar: Uncertain significance PVS1=0 PS=[0,0,0,0,0]",
        "InterVar: Likely benign PVS1=0",
        "InterVar: Benign PVS1=0",
        "",
        None,
    ]:
        assert classify_row(s) is False
