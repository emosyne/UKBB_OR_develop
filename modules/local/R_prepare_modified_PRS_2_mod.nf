process R_PREPARE_MODIF_PRS_2_LISTS {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(clumped_SNPs), path(original_dedup_GWAS), path (substituted_GWAS), path(tissue_EPeQTL_associations), path(tissue_facet_associations), path(all_TS_EPs_associations), path(merged_GWAS), path(all_TS_EPs_ZEROP_associations)
    // each path(GW_LD_blocks)

    output:
    tuple val(meta),  path("*_originalGWAS_SNPs_clumped_NoEPs_overlap.tsv"), path("*_TS_EPs_GWAS_REC_clumped_NonOverlapOriginal.tsv"), emit: split_GWASes
    // path("versions.yml")                        ,       emit: versions

    script:
    """
    R_prepare_modified_PRS2.R ${meta} ${clumped_SNPs} ${merged_GWAS}
    
    
   
    """
}
    