process R_ANNOTATE_ORs {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_high_resource_short'
    tag "$meta"
    cache 'lenient'

    input:
    //[SCZ, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, ENH_SCZ_hg19.csv, /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz]
    tuple val(meta), path(associations_recessive),  path(associations_standardGLM),  path(processed_SNPlists_hg19), path(full_GWAS)
    each path(hg38ToHg19_chain)
    each path(GW_LD_blocks)

    output:
    tuple val(meta), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    path "figs/*_UKBB.png"                      ,       emit: ORfigure
    

    script:
    """
    R_annotate_ORs.R ${meta} ${associations_recessive} ${associations_standardGLM} ${processed_SNPlists_hg19} ${full_GWAS} ${hg38ToHg19_chain} ${GW_LD_blocks}
    
   
    """
}
    