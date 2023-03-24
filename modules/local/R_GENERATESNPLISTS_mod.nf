process GENERATESNPLISTS {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(SNVs), path(GWAS), path(pheno)
    each path(chain_38_19)

    output:
    tuple val(meta), path("GWAS_*_hg19.bed"), path("ENH_*_hg19.bed"), path(pheno), path("ENH_*_hg19.csv"), emit: processed_ENH_SNP_lists_hg19

    script:
    """
    generateSNPlist.R ${SNVs} ${GWAS} ${chain_38_19} ${meta}
    
   
    """
}
