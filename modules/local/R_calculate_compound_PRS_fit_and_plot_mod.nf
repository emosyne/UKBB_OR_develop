process R_calculate_compound_PRS_fit_and_plot {
    // debug true
    container 'emosyne/r_docker:1.97'
    stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(TS_EPs_summary), path(TS_EPs_prsice), path (TS_EPs_best), path(original_NoEPsOverlap_summary), path(original_NoEPsOverlap_prsice), path(original_NoEPsOverlap_best), path(original_PRS_summary), path(original_PRS_prsice), path (original_PRS_best), path(caco_pheno)
    // each path(GW_LD_blocks)

    output:
    path("*")
    path("versions.yml")                        ,       emit: versions

    script:
    """
    R_calculate_compound_PRS_fit_and_plot.R ${meta} ${caco_pheno} ${TS_EPs_summary} ${TS_EPs_prsice} ${TS_EPs_best} ${original_NoEPsOverlap_summary} ${original_NoEPsOverlap_prsice} ${original_NoEPsOverlap_best} ${original_PRS_summary} ${original_PRS_best}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
    