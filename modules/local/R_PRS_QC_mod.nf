process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$condition"
    cache 'lenient'


    input: 
    // [SCZ, SCZ_ALLCHR.prune.in, SCZ_ALLCHR.het, SCZ_GWAS_QC_nodups.tsv.gz]
    tuple val(condition), path (prune), path (het), path (bedfilepath), path (bim), path (fam), path (log), path(GWAS)
    
    

    output:
    tuple val(condition), path ("*_het_valid_out.sample"), path("*_a1_bim"), path("*_mismatching_SNPs"),  emit: QC_het_a1_mismatch
    


    script:
    """
    R_PRS_QC.R ${condition} ${het} ${bim} ${GWAS} 
    
   
    """
}