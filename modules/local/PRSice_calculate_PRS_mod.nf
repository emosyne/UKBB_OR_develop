process PRSice_calculate_PRS_original {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/prsice_gwama_exec:1.0'
    cache 'lenient'


    input:
    // tuple val(meta), path(original_dedup_GWAS), path (substituted_GWAS), path(tissue_EPeQTL_associations), path(tissue_facet_associations), path(all_TS_EPs_associations), path(merged_GWAS), path(all_TS_EPs_ZEROP_associations)
    tuple val(meta), path(original_dedup_GWAS), path(substituted_GWAS), path(tissue_EPeQTL_associations), path(tissue_facet_associations), path(all_TS_EPs_associations), path(merged_GWAS), path(all_TS_EPs_ZEROP_associations), path (QCbed), path (QCbim), path (QCfam), path (QClog), path(caco_pheno), path(EUR_phase3_autosomes_hg19_bed), path(EUR_phase3_autosomes_hg19_bim), path(EUR_phase3_autosomes_hg19_fam)
    each path(UKBB_covariates)
    each path(UKBB_rs34380086_cases)

    output:
    tuple val(meta), path("*_original_PRS.summary"), path("*_original_PRS.prsice"), path("*_original_PRS.best"), emit: PRS_text_results
    tuple val(meta), path("*_original_PRS_rs34380086.summary"), path("*_original_PRS_rs34380086.prsice"), path("*_original_PRS_rs34380086.best")
    tuple val(meta), path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log
    // tuple val(meta), path ("*.*")
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    def mem_mb = (task.memory * 0.95).toMega()
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno
    
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${original_dedup_GWAS} \\
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
        --out ${QCbed.baseName}_original_PRS \\
        --thread $task.cpus \\
        --memory $mem_mb 
    
    #original PRS with --quant-extract rs34380086_cases
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${original_dedup_GWAS} \\
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
        --out ${QCbed.baseName}_original_PRS_rs34380086 \\
        --thread $task.cpus \\
        --memory $mem_mb 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    
    """
}



    
    
    // #modified substituted_GWAS PRS
    // PRSice.R \\
    //     --prsice /usr/local/bin/PRSice_linux \\
    //     --base ${substituted_GWAS} \\
    //     --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
    //     --keep-ambig --no-clump  \\
    //     --stat OR --or \\
    //     --target ${QCbed.baseName} \\
    //     --binary-target T --prevalence 0.01 \\
    //     --cov covariates.pheno \\
    //     --quantile 10 --quant-ref 1 \\
    //     --out ${QCbed.baseName}_modified_PRS_included_vars_lessp0.01 \\
    //     --thread $task.cpus \\
    //     --memory $mem_mb 
    
    // #modified substituted_GWAS PRS with --quant-extract rs34380086_cases
    // PRSice.R \\
    //     --prsice /usr/local/bin/PRSice_linux \\
    //     --base ${substituted_GWAS} \\
    //     --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
    //     --keep-ambig --no-clump  \\
    //     --stat OR --or \\
    //     --target ${QCbed.baseName} \\
    //     --binary-target T --prevalence 0.01 \\
    //     --cov covariates.pheno \\
    //     --quantile 10 --quant-ref 1 \\
    //     --quant-extract ${UKBB_rs34380086_cases} \\
    //     --out ${QCbed.baseName}_modified_PRS_included_vars_lessp0.01_rs34380086 \\
    //     --thread $task.cpus \\
    //     --memory $mem_mb 
    
    // #partitioned E-P_eQTL GWAS
    // PRSice.R \\
    //     --prsice /usr/local/bin/PRSice_linux \\
    //     --base ${all_TS_EPs_associations} \\
    //     --pheno ${caco_pheno} \\
    //     --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
    //     --keep-ambig \\
    //     --stat OR --or \\
    //     --target ${QCbed.baseName} \\
    //     --binary-target T --prevalence 0.01 \\
    //     --cov covariates.pheno \\
    //     --quantile 10 --quant-ref 1 \\
    //     --out ${QCbed.baseName}_all_TS_EPs_partition_PRS \\
    //     --thread $task.cpus \\
    //     --memory $mem_mb 