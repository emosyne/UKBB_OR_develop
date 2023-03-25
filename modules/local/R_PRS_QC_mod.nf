process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'


    input: 
    tuple val(meta), path (prune), path (het), path (bedfilepath), path (bim), path (fam), path (log), path(GWAS)
    
    

    output:
    tuple val(meta), path ("*_het_valid_out.sample"), path("*_a1_bim"), path("*_mismatching_SNPs"),  emit: QC_het_a1_mismatch
    path "versions.yml"           , emit: versions
    


    script:
    """
    R_PRS_QC.R ${meta} ${het} ${bim} ${GWAS} 
    
   
    """
}