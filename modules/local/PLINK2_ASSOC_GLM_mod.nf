process PLINK2_ASSOC_GLM {
    tag "$condition"
    // debug true
    label 'process_high_memory'
    cache 'lenient'
    container 'emosyne/plink2:1.23'
    // errorStrategy 'ignore'

    input: 
    // [SCZ, SCZ_ALLCHR_SCZ_QC.bed, SCZ_ALLCHR_SCZ_QC.bim, SCZ_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, Neural_significant_enh.bed]
    tuple val(condition), path (bed_QC), path (bim_QC), path (fam_QC), path(GWAS_QC), path(enhancers_bed)
    each path(UKBB_covariates)
    

    output:
    tuple val(condition), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid"), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid"), path("*_ORs_PLINK2_logistic_firth_fallback_covar_dominant.PHENO1.glm.logistic.hybrid"),path("*fallback_covar_recessive.frq"), emit: associations
    path("*.snplist")
    
    

    
    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bed_QC.baseName} \\
        --extract bed1 ${enhancers_bed} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm recessive firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_recessive
     plink \\
            --bfile ${bed_QC.baseName} \\
            --extract ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_recessive.snplist \\
            --freq \\
            --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_recessive

    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bed_QC.baseName} \\
        --extract bed1 ${enhancers_bed} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_standard
    plink \\
            --bfile ${bed_QC.baseName} \\
            --extract ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_standard.snplist \\
            --freq \\
            --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_standard
    
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bed_QC.baseName} \\
        --extract bed1 ${enhancers_bed} \\
        --chr 1-22 \\
        --write-snplist \\
        --glm dominant firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar ${UKBB_covariates} \\
        --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_dominant
     plink \\
            --bfile ${bed_QC.baseName} \\
            --extract ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_dominant.snplist \\
            --freq \\
            --out ${condition}_ORs_PLINK2_logistic_firth_fallback_covar_dominant

    """
}
