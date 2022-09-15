process PRSice_calculate_4_partition {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/prsice_gwama_exec:1.0'
    cache 'lenient'


    input:
    //out tuple val(meta),  path("*_original_NoEPs_associations.tsv"), path("*_NonOverlapOriginal_TS_EPs_associations.tsv"), emit: split_GWASes
    tuple val(meta), path(original_NoEPs_associations_clumped), path(NonOverlapOriginal_TS_EPs_associations_clumped), path (QCbed), path (QCbim), path (QCfam), path (QClog), path(caco_pheno), path(EUR_phase3_autosomes_hg19_bed), path(EUR_phase3_autosomes_hg19_bim), path(EUR_phase3_autosomes_hg19_fam)
    each path(UKBB_covariates)
    each path(UKBB_rs34380086_cases)

    output:
    tuple val(meta), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.summary"), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.prsice"), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.best"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.summary"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.prsice"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.best"), emit: PRS_text_split_results
    tuple val(meta), path("*.png"),  path("*.log") //figures, quantiles text and log
    // tuple val(meta), path ("*.*")
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    def mem_mb = task.memory.toMega()
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno

    #NonOverlapOriginal_TS_EPs_associations
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${NonOverlapOriginal_TS_EPs_associations_clumped} \\
        --pheno ${caco_pheno} \\
        --no-clump  \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
        --keep-ambig \\
        --stat OR --or \\
        --quantile 10 --quant-ref 1 \\
        --target ${QCbed.baseName} \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno --cov-factor sex \\
        --model rec \\
        --out ${QCbed.baseName}_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS \\
        --thread $task.cpus \\
        --memory $mem_mb 
    #original_NoEPs_associations
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${original_NoEPs_associations_clumped} \\
        --pheno ${caco_pheno} \\
        --no-clump  \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P  \\
        --keep-ambig \\
        --stat OR --or \\
        --quantile 10 --quant-ref 1 \\
        --target ${QCbed.baseName} \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno --cov-factor sex \\
        --model add \\
        --out ${QCbed.baseName}_ADD_original_NoEPs_associations_partition_PLINKclump_PRS \\
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



    // --no-full --bar-levels 0.05 --fastscore \\