process PLINK2_ASSOC_GLM {
    tag "$meta"
    // debug true
    label 'process_high_memory'
    cache 'lenient'
    container 'emosyne/plink2:1.2'
    errorStrategy 'ignore'

    input: 
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log), path(SNVs_hg19), path(pheno)
    each path (UKBB_covariates)
    

    output:
    tuple val(meta), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid"), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid"), emit: associations
    path("*")
    path "versions2.yml"           , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    //def prefix = task.ext.prefix 
    // if( "$PLINKbgenfiles" == "${prefix}.bgen" ) error "Input and output names are the same, use \"task.ext.prefix\" in modules.config to disambiguate!"
    def mem_mb = task.memory.toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bedfilepath.baseName} \\
        --extract bed1 ${SNVs_hg19} \\
        $args \\
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
        --extract bed1 ${SNVs_hg19} \\
        $args \\
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
      

    cat <<-END_VERSIONS > versions2.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
//         --pheno ${pheno} \\