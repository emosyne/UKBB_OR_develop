process PRSice_calculate_PRS_split_partitions {
    // debug true
    tag "${ENH_list}_${CTthreshold}_${EPWAS_model}"
    label 'process_high_resource_short'
    // clusterOptions "--partition=shared_52c_384g" //only for LISA
    container 'emosyne/prsice_gwama_exec:1.0'
    cache "lenient"
    // maxForks 5
    // errorStrategy 'ignore'


    input:
    // [SCZ_ALLCHR_SCZ_QC.bed, SCZ_ALLCHR_SCZ_QC.bim, SCZ_ALLCHR_SCZ_QC.fam, Neural_significant_enh, 
        // Neural_significant_enh_REC_SCZ_X_1_clumped_EPWAS.tsv.gz, Neural_significant_enh_REC_SCZ_clumped_residual_GWAS_compartment.tsv.gz, 1, SCZ, REC, 
        // SCZ_clumped_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/input/biobank/non_missing_10PCs_Jun22.covariate.gz, 
        // SCZ, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam, 
        // 0.5]
    tuple path(cohort_bed_QC),  path(cohort_bim_QC), path(cohort_fam_QC), val(ENH_list), \
        path(clumped_EPWAS), path(clumped_residual_GWAS_compartment), val(multiplier), val(condition),  val(EPWAS_model), \
        path(clumped_GWAS_QC_nodups), path(UKBB_covariates), \
        val(condition), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)\
        val(CTthreshold)
    
    output:
    tuple val("${ENH_list}_${CTthreshold}_${EPWAS_model}"), path("*_clumped_EPWAS_*.summary"), path("*_clumped_EPWAS_*.prsice"), path("*_clumped_EPWAS_*.best"), \
             emit: clumped_EPWAS_PRS
    tuple val("${ENH_list}_${CTthreshold}_${EPWAS_model}"), path("*_clumped_residual_GWAS_compartment.summary"), path("*_clumped_residual_GWAS_compartment.prsice"), path("*_clumped_residual_GWAS_compartment.best"), \
             emit: clumped_residual_GWAS_compartment_PRS
    tuple val("${ENH_list}_${CTthreshold}_${EPWAS_model}"), path("*_original_GWAS.summary"), path("*_original_GWAS.prsice"), path("*_original_GWAS.best"), val(ENH_list),  val(CTthreshold), val(condition),  path(cohort_fam_QC), val(EPWAS_model),        \
             emit: clumped_original_GWAS_PRS
    tuple  path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log

    script:
    def mem_Gb = (task.memory * 0.95).toGiga()
    def max_cpus = Math.round(task.cpus * 4/5)
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno
    echo memory: ${mem_Gb}Gb
    echo cpus: $max_cpus

    echo clumped_EPWAS - ORIGINAL OR
    #CHR    POS     SNP     A1      A2      P       OR      OR_by_measure1  OR_by_measure2  measure1       measure2
    # ORIGINAL OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_EPWAS} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${condition}_${ENH_list}_${CTthreshold}_${EPWAS_model}_clumped_EPWAS_originalOR

    echo clumped_EPWAS  - OR by measure 1
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_EPWAS} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure1 --or \\
        --out ${condition}_${ENH_list}_${CTthreshold}_${EPWAS_model}_mult_${multiplier}_clumped_EPWAS_OR_by_measure1
    
    echo clumped_EPWAS  - OR by measure 2
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_EPWAS} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure2 --or \\
        --out ${condition}_${ENH_list}_${CTthreshold}_${EPWAS_model}_mult_${multiplier}_clumped_EPWAS_OR_by_measure2

    echo clumped_residual_GWAS_compartment - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_residual_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${condition}_${ENH_list}_${CTthreshold}_${EPWAS_model}_clumped_residual_GWAS_compartment
    
    echo ORIGINAL GWAS LOO
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_GWAS_QC_nodups} \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat BETA --beta \\
        --target ${cohort_bed_QC.simpleName} \\
        --ld ${LD_ref_bed.baseName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --out ${condition}_${ENH_list}_original_GWAS
    """

}
        