process PLINK2_ASSOC_GLM {
    tag "$meta"
    // debug true
    label 'process_high_memory'
    cache 'lenient'
    container 'emosyne/plink2:1.23'
    // errorStrategy 'ignore'

    input: 
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log), path(SNVs_hg19), path(pheno)
    each path(UKBB_covariates)
    each path(Neural_significant_enh)
    

    output:
    tuple val(meta), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid"), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid"), path("*_ORs_PLINK2_logistic_firth_fallback_covar_dominant.PHENO1.glm.logistic.hybrid"), emit: associations
    path("*.snplist")
    path("*.frq")
    

    
    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bedfilepath.baseName} \\
        --extract bed1 ${Neural_significant_enh} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm recessive firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_recessive
     plink \\
            --bfile ${bedfilepath.baseName} \\
            --extract ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_recessive.snplist \\
            --freq \\
            --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_recessive

    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bedfilepath.baseName} \\
        --extract bed1 ${Neural_significant_enh} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_standard
    plink \\
            --bfile ${bedfilepath.baseName} \\
            --extract ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_standard.snplist \\
            --freq \\
            --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_standard
    
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bedfilepath.baseName} \\
        --extract bed1 ${Neural_significant_enh} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm dominant firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_dominant
     plink \\
            --bfile ${bedfilepath.baseName} \\
            --extract ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_dominant.snplist \\
            --freq \\
            --out ${meta}_ORs_PLINK2_logistic_firth_fallback_covar_dominant

    """
}
