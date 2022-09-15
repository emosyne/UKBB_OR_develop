process GENERATESNPLISTS {
    // debug true
    container 'emosyne/r_docker:1.7'
    stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(SNVs), path(GWAS), path(pheno)
    each path(chain_38_19)

    output:
    tuple val(meta), path("GWAS_*_hg19.bed"), path("ENH_*_hg19.bed"), path(pheno), path("ENH_*_hg19.csv"), emit: processed_ENH_SNP_lists_hg19
    path("versions.yml")         , emit: versions

    script:
    """
    generateSNPlist.R ${SNVs} ${GWAS} ${chain_38_19} ${meta}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
        R_GenomicRanges: \$(Rscript -e 'packageVersion("GenomicRanges")' | awk '{print \$2}' | tr -d "‘’")
        R_rtracklayer: \$(Rscript -e 'packageVersion("rtracklayer")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
