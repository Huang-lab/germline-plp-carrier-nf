process CARRIER_GT {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf), path(clinvar_tsv), path(acmg_tsv), path(am_tsv)

    output:
    tuple val(chunk_id), path("${chunk_id}.gt.tsv"), emit: gt

    script:
    """
    set -euo pipefail
    # Qualifying positions (CHROM<TAB>POS, 1-based) from whichever classifier
    # TSVs ran. Column located by header name, so column order is irrelevant.
    : > qualifying.txt
    for spec in "${clinvar_tsv}:is_clinvar_PLP" "${acmg_tsv}:is_acmg_PLP" "${am_tsv}:is_AM_PLP"; do
        f="\${spec%%:*}"; flag="\${spec##*:}"
        [ -s "\$f" ] || continue
        awk -F'\\t' -v flag="\$flag" '
            NR==1 { for (i=1;i<=NF;i++) c[\$i]=i; next }
            (c[flag] && \$(c[flag])==1) { print \$(c["chr"])"\\t"\$(c["pos"]) }
        ' "\$f" >> qualifying.txt
    done
    sort -u qualifying.txt > qualifying.pos.txt || true

    if [ ! -s qualifying.pos.txt ]; then
        # No qualifying variants in this chunk → empty GT table (header-only matrix downstream).
        : > ${chunk_id}.gt.tsv
    elif command -v bcftools >/dev/null 2>&1; then
        # C-speed: stream the VCF once, keep only qualifying positions (-T), emit per-sample GT.
        bcftools query -T qualifying.pos.txt \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' ${vcf} \\
          | awk -F'\\t' 'BEGIN{OFS="\\t"} {
                for (i=5; i<=NF; i++) { split(\$i, kv, "="); print \$1,\$2,\$3,\$4,kv[1],kv[2] }
            }' > ${chunk_id}.gt.tsv
    else
        # Test-profile fallback (no bcftools): pure-Python extraction.
        python3 ${projectDir}/bin/fixture_gt_extract.py \\
            --vcf ${vcf} --positions qualifying.pos.txt --out ${chunk_id}.gt.tsv
    fi
    """
}
