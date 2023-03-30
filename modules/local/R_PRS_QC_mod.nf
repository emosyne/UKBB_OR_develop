process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$condition"
    cache 'lenient'


    input: 
    //    val(condition), path (bedfilepath), path (bim), path (fam), path ("*.prune.in"), path ("*.het"),  emit: pruned_variants_het
    tuple val(condition), path(bed), path(bim), path(fam), path (prune), path (het), path(GWAS_QC)
    
    

    output:
    tuple val(condition), path(bed), path(bim), path(fam), path(GWAS_QC), path ("*_het_valid_out_vs_GWAS.sample"), path("*_a1_cohort_bim_vs_GWAS*"), path("*_mismatching_SNPs_vs_GWAS*"),  emit: QC_het_a1_mismatch
    


    script:
    """
    R_PRS_QC2.R ${het} ${bim} ${GWAS_QC} ${condition} 
    
   
    """
}