process CARRIER_MATRIX {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path clinvar_tsv
    path acmg_tsv
    path am_tsv
    path keep_samples

    output:
    path("${chunk_id}.carrier_matrix.tsv"), emit: long_tsv

    script:
    def keep_arg = keep_samples ? "--keep ${keep_samples}" : ""
    """
    set -euo pipefail
    # Subset annotated VCF to union of qualifying variants (once), then single bcftools query.
    python3 - <<'PY'
import csv
keys = set()
for p, f in [("${clinvar_tsv}", "is_clinvar_PLP"),
             ("${acmg_tsv}", "is_acmg_PLP"),
             ("${am_tsv}", "is_AM_PLP")]:
    with open(p) as fh:
        r = csv.DictReader(fh, delimiter="\\t")
        for row in r:
            if int(row.get(f, "0") or 0):
                keys.add((row["chr"], row["pos"], row["ref"], row["alt"]))
with open("qualifying.bed", "w") as out:
    for c, p, ref, alt in sorted(keys):
        start = int(p) - 1
        end = int(p) + max(len(ref), len(alt)) - 1
        out.write(f"{c}\\t{start}\\t{end}\\n")
PY

    if command -v bcftools >/dev/null 2>&1; then
        bcftools query -R qualifying.bed \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' ${vcf} \\
          | awk -F'\\t' 'BEGIN{OFS="\\t"} {
                for(i=5; i<=NF; i++) {
                    split(\$i, kv, "=");
                    print \$1,\$2,\$3,\$4,kv[1],kv[2]
                }
            }' > gt.tsv
    else
        # Test-profile fallback: derive per-sample GT from the fixture VCF directly.
        python3 ${projectDir}/bin/fixture_gt_extract.py --vcf ${vcf} --bed qualifying.bed --out gt.tsv
    fi

    ${projectDir}/bin/build_carrier_matrix.py \\
        --gt gt.tsv \\
        --clinvar ${clinvar_tsv} \\
        --acmg ${acmg_tsv} \\
        --am ${am_tsv} \\
        ${keep_arg} \\
        --out-long ${chunk_id}.carrier_matrix.tsv
    """
}
