process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$condition"
    cache 'lenient'


    input: 
    //    val(condition), path (bedfilepath), path (bim), path (fam), path ("*.prune.in"), path ("*.het"),  emit: pruned_variants_het
    tuple val(condition), path (bedfilepath), path (bim), path (fam), path (prune), path (het), path(GWAS)
    
    

    output:
    tuple val(condition), path ("*_het_valid_out.sample"), path("*_a1_bim"), path("*_mismatching_SNPs"),  emit: QC_het_a1_mismatch
    


    script:
    """
    R_PRS_QC.R ${condition} ${het} ${bim} ${GWAS} 
    
   
    """
}