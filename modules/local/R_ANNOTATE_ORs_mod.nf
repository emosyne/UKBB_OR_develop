process R_ANNOTATE_ORs {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_high_resource_short'
    tag "$condition"
    cache 'lenient'

    input:
    //[SCZ, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, 
        // SCZ_ORs_PLINK2_logistic_firth_fallback_covar_dominant.PHENO1.glm.logistic.hybrid, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.frq, 
        // /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz]
    tuple val(condition), path(PLINK_ORs_recessive),  path(PLINK_ORs_additive),  \
        path(PLINK_ORs_dominant),   path(SNP_frq), path(full_GWAS)


    output:
    path("*_EPWAS_SNPs.bed"),                   emit: EPWAS_SNPs
    tuple path("*.tsv.gz"),                     emit: EPWAS_files_HCM_format
    path "figs/*_UKBB.pdf"
    

    script:
    """
    R_annotate_ORs.R ${condition} ${PLINK_ORs_recessive} ${PLINK_ORs_dominant} ${PLINK_ORs_additive} ${full_GWAS} ${SNP_frq} 
    
   
    """
}
    