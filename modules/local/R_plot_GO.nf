process R_plot_GO {
    // debug true
    container 'emosyne/r_docker:1.97'
    // stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(annotated_ORs), path(full_GWAS_hg19), path(ENH_PROM_hg38)

    output:
    path("*")
    

    script:
    """
    R_plot_GO.R ${meta} ${annotated_ORs} ${full_GWAS_hg19} ${ENH_PROM_hg38}
    

   
    """
}
    