process PER_GENE_QC {
    tag "${chunk_id}"

    input:
    tuple val(chunk_id), path(vcf)
    path clinvar_tsv
    path acmg_tsv
    path am_tsv
    val classifiers_csv

    output:
    path("${chunk_id}.qc_per_gene.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
import csv, os
from collections import Counter
inputs = []
for path, flag in [("${clinvar_tsv}", "is_clinvar_PLP"),
                   ("${acmg_tsv}",    "is_acmg_PLP"),
                   ("${am_tsv}",      "is_AM_PLP")]:
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        inputs.append((path, flag))
counts = Counter(); plps = Counter()
for path, flag in inputs:
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\\t")
        for row in r:
            g = row.get("gene", "") or "."
            counts[g] += 1
            if int(row.get(flag, "0") or 0):
                plps[(g, flag)] += 1
with open("${chunk_id}.qc_per_gene.tsv", "w") as out:
    out.write("gene\\tn_annotated\\tn_clinvar_PLP\\tn_acmg_PLP\\tn_AM_PLP\\n")
    for g in sorted(counts):
        out.write(f"{g}\\t{counts[g]}\\t{plps[(g,'is_clinvar_PLP')]}\\t{plps[(g,'is_acmg_PLP')]}\\t{plps[(g,'is_AM_PLP')]}\\n")
PY
    """
}
