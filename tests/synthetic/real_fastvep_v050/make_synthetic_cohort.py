#!/usr/bin/env python3
"""Synthetic genome, GFF3 and cohort VCF for tests/test_real_fastvep_output.py.

No real data: seeded random sequence, designed genes, 10 invented samples.

usage: make_synthetic_cohort.py OUTDIR   writes synth.fa synth.gff3 cohort.vcf truth.tsv and sa_src/ (see README.md)

Contig "1" (60,000 bp) and contig "X" (20,000 bp), Ensembl style names, + strand genes only.
Every gene follows one template (o = offset):
  exon1 o+1..o+200   5'UTR o+1..o+50   CDS o+51..o+200  (codon 1 ATG, codon 2 CAG, codon 3 GGA, 47 random sense codons)
  intron o+201..o+500 (donor GT at o+201..202, acceptor AG at o+499..500)
  exon2 o+501..o+700 CDS o+501..o+650 (last codon TAA at o+648..o+650), 3'UTR o+651..o+700
C>T at o+54 turns codon 2 CAG into TAG: stop_gained, exon 1 of 2, so NMD-triggering.

Genes (contig 1 unless noted):
  GENEA  o=1000   MANE Select + canonical
  GENEB  1101..1400 single exon, 5'UTR 1101..1250 (overlaps GENEA's CDS tail), CDS 1251..1400, MANE Select
  GENEC  o=5000   NO MANE; T1 canonical, T2 non-canonical copy of T1
  GENED  o=8000   MANE Plus Clinical only (T1, canonical)
  GENEF  o=12000  MANE Select T1 (canonical) + non-MANE T2 with an extra in-frame exon o+301..o+399
  GENEG  o=16000  MANE Select (used for the multi-allelic record)
  GENEH  o=20000  NO MANE, isolated (stop_gained there must still be reported by a pick-based call)
  GENEI  o=25000  NO MANE, 560 bp downstream of GENEJ (MANE): the masked-LoF case
  GENEJ  o=26400  MANE Select (flanks GENEI)  [placed 3.1 kb upstream-adjacent is not needed: 700 bp gap]
  GENEX  contig X, o=1000, MANE Select (non-PAR style hemizygous test)
"""
import random
import sys
from pathlib import Path

SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)
    contigs = {"1": [rng.choice("ACGT") for _ in range(60000)], "X": [rng.choice("ACGT") for _ in range(20000)]}

    def put(ctg, pos1, s):
        for i, ch in enumerate(s):
            contigs[ctg][pos1 - 1 + i] = ch

    def get(ctg, pos1, n=1):
        return "".join(contigs[ctg][pos1 - 1:pos1 - 1 + n])

    gff = ["##gff-version 3"]
    n_gene = [0]

    def add_gene(ctg, name, tx_specs):
        """tx_specs: list of dict(tid, tags, exons=[(s,e)], cds=[(s,e)], utr5=[(s,e)], utr3=[(s,e)])"""
        n_gene[0] += 1
        gid = f"ENSG{n_gene[0]:011d}"
        gs = min(e[0] for t in tx_specs for e in t["exons"])
        ge = max(e[1] for t in tx_specs for e in t["exons"])
        gff.append(f"{ctg}\tensembl\tgene\t{gs}\t{ge}\t.\t+\t.\tID=gene:{gid};Name={name};biotype=protein_coding;gene_id={gid};version=1")
        for t in tx_specs:
            tid = t["tid"]
            ts = min(e[0] for e in t["exons"])
            te = max(e[1] for e in t["exons"])
            tag = f";tag={t['tags']}" if t["tags"] else ""
            gff.append(f"{ctg}\tensembl\tmRNA\t{ts}\t{te}\t.\t+\t.\tID=transcript:{tid};Parent=gene:{gid};Name={name}-{tid[-3:]};biotype=protein_coding{tag};transcript_id={tid};version=1")
            for r, (s, e) in enumerate(t["exons"], 1):
                gff.append(f"{ctg}\tensembl\texon\t{s}\t{e}\t.\t+\t.\tParent=transcript:{tid};Name={tid}.e{r};exon_id={tid}.e{r};rank={r};version=1")
            for (s, e) in t["utr5"]:
                gff.append(f"{ctg}\tensembl\tfive_prime_UTR\t{s}\t{e}\t.\t+\t.\tParent=transcript:{tid}")
            for (s, e) in t["cds"]:
                gff.append(f"{ctg}\tensembl\tCDS\t{s}\t{e}\t.\t+\t0\tID=CDS:{tid}P;Parent=transcript:{tid};protein_id={tid}P;version=1")
            for (s, e) in t["utr3"]:
                gff.append(f"{ctg}\tensembl\tthree_prime_UTR\t{s}\t{e}\t.\t+\t.\tParent=transcript:{tid}")
        return gid

    tcount = [0]

    def tid():
        tcount[0] += 1
        return f"ENST{tcount[0]:011d}"

    def template(ctg, name, o, tags_list, extra_exon=False):
        cds1 = ["ATG", "CAG", "GGA"] + [rng.choice(SENSE) for _ in range(47)]
        cds2 = [rng.choice(SENSE) for _ in range(49)] + ["TAA"]
        put(ctg, o + 51, "".join(cds1))
        put(ctg, o + 501, "".join(cds2))
        put(ctg, o + 201, "GT")
        put(ctg, o + 203, "AAGT")
        put(ctg, o + 499, "AG")
        put(ctg, o + 480, "TTTCTTTCTTTC")
        if extra_exon:
            ex = ["GAA", "CAG"] + [rng.choice(SENSE) for _ in range(31)]
            put(ctg, o + 301, "".join(ex))
            put(ctg, o + 299, "AG")
            put(ctg, o + 400, "GT")
        txs = []
        for k, tags in enumerate(tags_list):
            t_extra = extra_exon and k == 1
            if t_extra:
                txs.append(dict(tid=tid(), tags=tags, exons=[(o + 1, o + 200), (o + 301, o + 399), (o + 501, o + 700)],
                                cds=[(o + 51, o + 200), (o + 301, o + 399), (o + 501, o + 650)],
                                utr5=[(o + 1, o + 50)], utr3=[(o + 651, o + 700)]))
            else:
                txs.append(dict(tid=tid(), tags=tags, exons=[(o + 1, o + 200), (o + 501, o + 700)],
                                cds=[(o + 51, o + 200), (o + 501, o + 650)],
                                utr5=[(o + 1, o + 50)], utr3=[(o + 651, o + 700)]))
        add_gene(ctg, name, txs)

    # contig 1
    template("1", "GENEA", 1000, ["Ensembl_canonical,MANE_Select"])
    # GENEB: single exon, 5'UTR overlaps GENEA CDS tail
    cdsb = ["ATG"] + [rng.choice(SENSE) for _ in range(48)] + ["TAA"]
    put("1", 1251, "".join(cdsb))
    add_gene("1", "GENEB", [dict(tid=tid(), tags="Ensembl_canonical,MANE_Select", exons=[(1101, 1400)],
                                 cds=[(1251, 1400)], utr5=[(1101, 1250)], utr3=[])])
    template("1", "GENEC", 5000, ["Ensembl_canonical", ""])
    template("1", "GENED", 8000, ["Ensembl_canonical,MANE_Plus_Clinical"])
    template("1", "GENEF", 12000, ["Ensembl_canonical,MANE_Select", ""], extra_exon=True)
    template("1", "GENEG", 16000, ["Ensembl_canonical,MANE_Select"])
    template("1", "GENEH", 20000, ["Ensembl_canonical"])
    template("1", "GENEI", 25000, ["Ensembl_canonical"])
    template("1", "GENEJ", 26400, ["Ensembl_canonical,MANE_Select"])
    template("X", "GENEX", 1000, ["Ensembl_canonical,MANE_Select"])

    with open(out / "synth.fa", "w") as fh:
        for c, s in contigs.items():
            fh.write(f">{c}\n")
            seq = "".join(s)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    (out / "synth.gff3").write_text("\n".join(gff) + "\n")

    # cohort VCF: 10 samples, GT:DP:GQ
    samples = [f"S{i}" for i in range(1, 11)]
    contig_hdr = "##contig=<ID=1,length=60000>\n##contig=<ID=X,length=20000>\n"
    hdr = ("##fileformat=VCFv4.2\n" + contig_hdr +
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
           '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
           '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

    def rec(ctg, pos, ref, alt, gts, dp="30", gq="60"):
        assert len(gts) == 10, (ctg, pos, gts)
        cols = "\t".join(f"{g}:{dp}:{gq}" for g in gts)
        return (ctg, pos, f"{ctg}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:DP:GQ\t{cols}\n")

    het10 = ["0/1", "1/1", "1", "./1", "0|1", "./.", "0/0", "1/.", "1|1", "0"]   # the genotype zoo
    simple = ["0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]
    recs = []
    truth = []

    def stop(ctg, o, label, gts, extra=0, alt="T"):
        pos = o + 54 + extra
        assert get(ctg, pos) == "C", (label, get(ctg, pos))
        recs.append(rec(ctg, pos, "C", alt, gts))
        truth.append((label, ctg, pos, "C", alt))

    stop("1", 1000, "A_stop_zoo", het10)                      # GENEA MANE stop_gained, every genotype shape
    stop("1", 5000, "C_stop_noMANE_2tx", simple)              # GENEC no MANE gene, two transcripts, stop in both
    stop("1", 8000, "D_stop_PlusClinical_only", simple)       # GENED MANE Plus Clinical only
    stop("1", 20000, "H_stop_isolated_noMANE", simple)        # GENEH isolated, no MANE
    stop("1", 25000, "I_stop_noMANE_next_to_MANE", simple)    # GENEI no MANE, GENEJ (MANE) 700 bp away
    stop("X", 1000, "X_stop_chrX", het10)                     # GENEX on X
    # GENEF: stop_gained only in the non-MANE transcript's extra exon (o+304, codon 2 of the extra exon)
    pos = 12000 + 304
    assert get("1", pos) == "C", get("1", pos)
    recs.append(rec("1", pos, "C", "T", simple))
    truth.append(("F_stop_in_nonMANE_exon_only", "1", pos, "C", "T"))
    # multi-allelic record in GENEG: allele 1 = stop_gained (T), allele 2 = missense (G)
    pos = 16054
    assert get("1", pos) == "C"
    recs.append(rec("1", pos, "C", "T,G", ["0/2", "0/1", "1/2", "2/2", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]))
    truth.append(("G_multiallelic_T_stop_G_missense", "1", pos, "C", "T,G"))
    # variant in GENEA coding tail / GENEB 5'UTR overlap: a stop in GENEA at 1100+ ? use a missense-able codon: pos 1160
    # GENEA CDS 1051..1200; codon boundaries start at 1051 + 3k. 1159 = 1051+108 -> codon 37 start. Make it a stop: set to TAT? use C>T on a CAG-like base.
    # simply place a SNV, whatever the codon, and let fastVEP say what it is (truth: two genes, GENEA coding, GENEB 5'UTR)
    pos = 1160
    recs.append(rec("1", pos, get("1", pos), {"A": "G", "C": "T", "G": "A", "T": "C"}[get("1", pos)], simple))
    truth.append(("AB_overlap_snv", "1", pos, get("1", pos), "?"))
    # frameshift insertion in GENEA exon 1 (1 base), padded representation
    pos = 1100
    recs.append(rec("1", pos, get("1", pos), get("1", pos) + "A", simple))
    truth.append(("A_frameshift_ins", "1", pos, get("1", pos), get("1", pos) + "A"))
    recs.sort(key=lambda r: ({"1": 1, "X": 2}[r[0]], r[1]))
    (out / "cohort.vcf").write_text(hdr + "".join(r[2] for r in recs))
    with open(out / "truth.tsv", "w") as fh:
        fh.write("label\tchr\tpos\tref\talt\n")
        for t in truth:
            fh.write("\t".join(map(str, t)) + "\n")
    print("wrote", out)


def write_sa_src(d):
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    genes = [f"GENE{c}" for c in "ABCDEFGHIJ"] + ["GENEX"]
    with open(d / "constraint.tsv", "w") as fh:
        fh.write("gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z\n")
        for i, g in enumerate(genes):
            fh.write(f"{g}\tENST9{i:010d}\t0\t40.0\t0.00\t0.10\t1.0000\t3.0\t0.1\n")
    with open(d / "genemap2.txt", "w") as fh:
        fh.write("# Generated\n# Copyright OMIM\n")
        for i, g in enumerate(genes):
            inh = "X-linked recessive" if g == "GENEX" else "Autosomal dominant"
            fh.write(f"1\tp36\t1:1-2\tGene\t\t{600000+i}\tprotein\t\t{g}\t\t\t\tSynthetic cancer syndrome {i}, {610000+i} (3), {inh}\n")
    with open(d / "gnomad.vcf", "w") as fh:
        fh.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=60000>\n##contig=<ID=X,length=20000>\n')
        fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="af">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="an">\n')
        fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">\n##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="hc">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t30010\t.\tA\tG\t.\tPASS\tAF=0.2;AN=100000;AC=20000;nhomalt=2000\n")


if __name__ == "__main__":
    main(sys.argv[1])
    write_sa_src(Path(sys.argv[1]) / "sa_src")
