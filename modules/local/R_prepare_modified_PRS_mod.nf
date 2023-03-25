process R_PREPARE_MODIF_PRS {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(full_GWAS_hg19), path (annotated_ORs)
    // each path(GW_LD_blocks)

    output:
    tuple val(meta),  path("*_original_dedup_GWAS.tsv"), path("*_substituted_GWAS.tsv"), path("*_only_tissue_EPeQTL_associations.tsv"), path("*_only_tissue_facet_associations.tsv"), path("*_all_TS_EPs_associations.tsv"), path("*_merged_GWAS.tsv"), path("*_all_TS_EPs_associations_Pdivide250.tsv"), emit: orig_and_modified_GWASes
    // path("*.csv")
    // path("versions.yml")                        ,       emit: versions

    script:
    """
    R_prepare_modified_PRS.R ${meta} ${full_GWAS_hg19} ${annotated_ORs}
    

   
    """
}
    