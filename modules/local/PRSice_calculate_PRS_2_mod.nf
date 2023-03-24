process PRSice_calculate_PRS_substituted {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/r_docker:1.97'
    cache 'lenient'


    input:
    tuple val(meta), path(original_dedup_GWAS), path(substituted_GWAS), path(tissue_EPeQTL_associations), path(tissue_facet_associations), path(all_TS_EPs_associations), path(merged_GWAS), path (QCbed), path (QCbim), path (QCfam), path (QClog), path(caco_pheno), path(EUR_phase3_autosomes_hg19_bed), path(EUR_phase3_autosomes_hg19_bim), path(EUR_phase3_autosomes_hg19_fam)
    each path(UKBB_covariates)
    each path(UKBB_rs34380086_cases)

    output:
    tuple val(meta), path("*_substituted_vars_less0.05_PRS.summary"), path("*_substituted_vars_less0.05_PRS.prsice"), path("*_substituted_vars_less0.05_PRS.best"), emit: PRS_text_results
    tuple val(meta), path("*_substituted_vars_less0.05_PRS_rs34380086.summary"), path("*_substituted_vars_less0.05_PRS_rs34380086.prsice"), path("*_substituted_vars_less0.05_PRS_rs34380086.best")
    tuple val(meta), path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log
    // tuple val(meta), path ("*.*")
    // path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    def mem_mb = (task.memory * 0.95).toMega()
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno
    
    #modified substituted_GWAS PRS
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${substituted_GWAS} \\
        --pheno ${caco_pheno} \\
        --ld ${EUR_phase3_autosomes_hg19_bed.baseName} \\
        --clump-kb 3M   \\
        --clump-p 1     \\
        --clump-r2 0.1  \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
        --keep-ambig \\
        --stat OR --or \\
        --target ${QCbed.baseName} \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno --cov-factor sex \\
        --quantile 10 --quant-ref 1 \\
        --out ${QCbed.baseName}_substituted_vars_less0.05_PRS \\
        --thread $task.cpus \\
        --memory $mem_mb 
    
    #modified substituted_GWAS PRS with --quant-extract rs34380086_cases
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${substituted_GWAS} \\
        --pheno ${caco_pheno} \\
        --ld ${EUR_phase3_autosomes_hg19_bed.baseName} \\
        --clump-kb 3M   \\
        --clump-p 1     \\
        --clump-r2 0.1  \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
        --keep-ambig \\
        --stat OR --or \\
        --target ${QCbed.baseName} \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno --cov-factor sex \\
        --quantile 10 --quant-ref 1 \\
        --quant-extract ${UKBB_rs34380086_cases} \\
        --out ${QCbed.baseName}_substituted_vars_less0.05_PRS_rs34380086 \\
        --thread $task.cpus \\
        --memory $mem_mb 
    

    
    
    """
}



    