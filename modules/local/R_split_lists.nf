process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "${ENH_list}_${EPWAS_model}"
    cache "lenient"
    // errorStrategy 'ignore'
    

    input:
    //     tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment), path("*clumped_SNPs.clumped"), \
            // val(condition), val(EPWAS_model),   emit: clumped_SNPs_and_noclump_lists
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
    