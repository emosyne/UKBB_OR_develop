process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "${ENH_list}_${EPWAS_model}"
    cache "lenient"
    // errorStrategy 'ignore'
    

    input:
    //  [SCZ_ALLCHR_SCZ_QC.bed, SCZ_ALLCHR_SCZ_QC.bim, SCZ_ALLCHR_SCZ_QC.fam, Neural_significant_enh, 
        // SCZ_REC_Neural_significant_enh_noclump_EPWAS.tsv.gz, SCZ_REC_Neural_significant_enh_noclump_residual_GWAS_compartment.tsv.gz, SCZ_Neural_significant_enh_clumped_SNPs.clumped, SCZ, REC]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(condition), val(EPWAS_model)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path("*_clumped_EPWAS.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), val("1"), \
        val(condition), val(EPWAS_model),      emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${ENH_list}_${EPWAS_model}_${condition}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_EPWAS} ${EP_ES_gene_brain_exp} "1"
    
    
   
    """
}
    