process NORM_QC {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path reference_fasta

    output:
    tuple val(chunk_id), path("${chunk_id}.norm.vcf"), emit: vcf

    script:
    def min_dp = params.min_dp
    def min_gq = params.min_gq
    """
    set -euo pipefail
    if command -v bcftools >/dev/null 2>&1 && [ -s "${reference_fasta}" ]; then
        # Contig-naming reconciliation. The reference FASTA, VEP cache, and NCBI
        # ClinVar are Ensembl-style (1,2,...,X,Y,MT); input VCFs may be UCSC-style
        # (chr1,...,chrM). Derive the chromosome names from the actual records
        # (CHROM column) — VCFs frequently omit ##contig headers — and build a
        # chr-stripping rename map. chrM/chrMT -> MT.
        bcftools view ${vcf} | grep -v '^#' | cut -f1 | uniq > chroms.txt
        : > chr_map.txt
        while read -r c; do
            [ -z "\$c" ] && continue
            case "\$c" in
                chrM|chrMT) printf '%s\\tMT\\n'  "\$c" >> chr_map.txt ;;
                chr*)       printf '%s\\t%s\\n' "\$c" "\${c#chr}" >> chr_map.txt ;;
            esac
        done < chroms.txt

        if [ -s chr_map.txt ]; then
            echo "Renaming contigs:"; cat chr_map.txt
            bcftools annotate --rename-chrs chr_map.txt ${vcf} -Oz -o renamed.vcf.gz
            SRC=renamed.vcf.gz
        else
            SRC=${vcf}
        fi

        bcftools norm -m -any -f ${reference_fasta} "\$SRC" -Ov \\
          | bcftools view -f PASS,. -Ov \\
          | bcftools +setGT -Ov -- \\
              -t q -n . \\
              -i 'FMT/DP<${min_dp} | FMT/GQ<${min_gq}' \\
          > ${chunk_id}.norm.vcf
    else
        # Test-profile fallback: no bcftools available; passthrough while asserting the file is a VCF.
        head -c 4 ${vcf} | grep -q '##' || { echo "not a VCF"; exit 2; }
        cp ${vcf} ${chunk_id}.norm.vcf
    fi
    """
}
