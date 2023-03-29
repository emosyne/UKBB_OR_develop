process R_extract_GWAS_SNPs_into_bed {
    // debug true
    maxForks 2
    container 'emosyne/r_docker:1.97'
    // container 'emosyne/simpler:1.1'
    label 'process_high_resource_short'
    // tag "$meta"
    cache "lenient"

    
    input:
    path(collected_bed_files_for_enhancers)
    tuple path(GWAS_QC_nodups), path(GWAS_QC_nodups_clumped_SNPs), val(condition)
    

    output:
    path("*clumped_GWAS_SNPs_plus_those_in_bed_files.bed"), emit: clumped_GWAS_SNPs_plus_those_in_bed_files
    path("*clumped_GWAS_QC_nodups.tsv.gz"),      emit: clumped_GWAS
    

    script:
    """
    echo ${collected_bed_files_for_enhancers} | tr ' ' '\\n' > collected_bed_files_for_enhancers.txt

    R_extract_GWAS_SNPs_into_bed.R ${GWAS_QC_nodups} collected_bed_files_for_enhancers.txt ${GWAS_QC_nodups_clumped_SNPs} ${condition}
    
    """
}
