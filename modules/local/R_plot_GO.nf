process R_plot_GO {
    // debug true
    container 'emosyne/r_docker:1.97'
    stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(annotated_ORs), path(full_GWAS_hg19), path(ENH_PROM_hg38)

    output:
    path("*")
    path("versions.yml")                        ,       emit: versions

    script:
    """
    R_plot_GO.R ${meta} ${annotated_ORs} ${full_GWAS_hg19} ${ENH_PROM_hg38}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
    