process CARRIER_MATRIX {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path clinvar_tsv
    path acmg_tsv
    path am_tsv
    path keep_samples
    val classifiers_csv

    output:
    path("${chunk_id}.carrier_matrix.tsv"), emit: long_tsv

    script:
    def has_cv = clinvar_tsv?.size() > 0 ? "--clinvar ${clinvar_tsv}" : ""
    def has_ac = acmg_tsv?.size()    > 0 ? "--acmg ${acmg_tsv}"       : ""
    def has_am = am_tsv?.size()      > 0 ? "--am ${am_tsv}"           : ""
    def keep_arg = (keep_samples && keep_samples.size() > 0) ? "--keep ${keep_samples}" : ""
    """
    set -euo pipefail
    # Union of qualifying variants (across whichever classifiers ran) → BED.
    python3 - <<'PY'
import csv, os
keys = set()
for p, f in [("${clinvar_tsv}", "is_clinvar_PLP"),
             ("${acmg_tsv}",    "is_acmg_PLP"),
             ("${am_tsv}",      "is_AM_PLP")]:
    if not p or not os.path.exists(p) or os.path.getsize(p) == 0:
        continue
    with open(p) as fh:
        for row in csv.DictReader(fh, delimiter="\\t"):
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
        python3 ${projectDir}/bin/fixture_gt_extract.py --vcf ${vcf} --bed qualifying.bed --out gt.tsv
    fi

    ${projectDir}/bin/build_carrier_matrix.py \\
        --gt gt.tsv \\
        ${has_cv} ${has_ac} ${has_am} \\
        ${keep_arg} \\
        --out-long ${chunk_id}.carrier_matrix.tsv
    """
}
