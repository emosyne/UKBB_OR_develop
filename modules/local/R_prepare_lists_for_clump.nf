process R_prepare_lists_for_clump {
    // debug true
    // errorStrategy 'terminate'
    container 'emosyne/r_docker:1.97'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${ENH_list}_${EPWAS_model}"
    cache "lenient"
    

    input:
    // [SCZ, SCZ_ALLCHR_SCZ_QC.bed, SCZ_ALLCHR_SCZ_QC.bim, SCZ_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, Neural_significant_enh, Neural_significant_enh.bed, 
        // REC, UKBB_ENH_associations_REC.tsv.gz]
    tuple val(condition), path(bed_QC),  path(bim_QC), path(fam_QC), path (QCed_GWAS),  val(ENH_list), path(ENH_bed), \
        val(EPWAS_model), path(ENH_EPwas)

    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path("*_noclump_EPWAS.tsv.gz"), path("*_noclump_residual_GWAS_compartment.tsv.gz"),  \
        val(condition), val(EPWAS_model),  emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${QCed_GWAS}  ${condition} ${EPWAS_model} ${ENH_EPwas}
    
   
    """
}
    
    