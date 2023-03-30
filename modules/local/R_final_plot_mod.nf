process R_final_plot {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high_short'
    tag "$ENHlist_thresh_model"
    cache "lenient"
    // errorStrategy 'ignore'

    input: 
    // [Neural_significant_enh_0.05_REC, 
        // SCZ_Neural_significant_enh_0.05_REC_clumped_EPWAS_originalOR.summary, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure1.summary, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure2.summary, 
        // SCZ_Neural_significant_enh_0.05_REC_clumped_EPWAS_originalOR.prsice, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure1.prsice, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure2.prsice, 
        // SCZ_Neural_significant_enh_0.05_REC_clumped_EPWAS_originalOR.best, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure1.best, SCZ_Neural_significant_enh_0.05_REC_mult_1_clumped_EPWAS_OR_by_measure2.best, 
        // SCZ_Neural_significant_enh_0.05_REC_clumped_residual_GWAS_compartment.summary, SCZ_Neural_significant_enh_0.05_REC_clumped_residual_GWAS_compartment.prsice, SCZ_Neural_significant_enh_0.05_REC_clumped_residual_GWAS_compartment.best,
        // SCZ_Neural_significant_enh_original_GWAS.summary, SCZ_Neural_significant_enh_original_GWAS.prsice, SCZ_Neural_significant_enh_original_GWAS.best, 
        // Neural_significant_enh, 0.05, SCZ, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, REC, enh_ES, enh_TS_tpm]      
    tuple val(ENHlist_thresh_model), \
        path(EPWAS_originalOR_summary), path(EPWAS_OR_by_measure1_summary), path(EPWAS_OR_by_measure2_summary),  \
        path(EPWAS_originalOR_prsice), path(EPWAS_OR_by_measure1_prsice), path(EPWAS_OR_by_measure2_prsice),  \
        path(EPWAS_originalOR_best), path(EPWAS_OR_by_measure1_best), path(EPWAS_OR_by_measure2_best),  \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(original_GWAS_summary), path(original_GWAS_prsice), path (original_GWAS_best),\
        val(ENH_list), val(CTthreshold), val(condition), path(cohort_fam), val(model),\
        val(modif_name_1),val(modif_name_2)

    output:
    // path("*.txt")
    path("*/*.pdf")
    path("*/*.log")

    script:
    """
    
    R_final_plot.R $task.cpus ${ENH_list} ${cohort_fam} \
        ${EPWAS_originalOR_summary} ${EPWAS_originalOR_best}\
        ${EPWAS_OR_by_measure1_summary} ${EPWAS_OR_by_measure1_best}\
        ${EPWAS_OR_by_measure2_summary} ${EPWAS_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${EPWAS_originalOR_prsice} ${EPWAS_OR_by_measure1_prsice} ${EPWAS_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice}   \
        ${original_GWAS_summary} ${original_GWAS_prsice} ${original_GWAS_best}\
        ${modif_name_1} ${modif_name_2} ${CTthreshold} ${condition} ${ENHlist_thresh_model} 
    """
}
    