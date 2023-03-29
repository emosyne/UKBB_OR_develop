process R_ANNOTATE_ORs {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_high_resource_short'
    tag "$condition"
    cache 'lenient'

    input:
    //[SCZ, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, 
        // ENH_SCZ_hg19.csv, /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz]
    tuple val(condition), path(PLINK_ORs_recessive),  path(PLINK_ORs_additive),  \
        path(PLINK_ORs_dominant),  path(full_GWAS), path(SNP_frq)


    output:
    tuple val(condition), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    path "figs/*_UKBB.png"                      ,       emit: ORfigure
    

    script:
    """
    R_annotate_ORs.R ${condition} ${PLINK_ORs_recessive} ${PLINK_ORs_dominant} ${PLINK_ORs_additive} ${full_GWAS} ${SNP_frq} 
    
   
    """
}
    