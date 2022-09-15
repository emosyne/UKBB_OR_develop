process R_PRS_PPV_plotting {
    // debug true
    container 'emosyne/r_docker:1.7'
    stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(original_PRS_summary), path(original_PRS_prsice), path (original_PRS_best), path(TS_EPs_summary), path(TS_EPs_prsice), path (TS_EPs_best), path(original_NoEPsOverlap_summary), path(original_NoEPsOverlap_prsice), path(original_NoEPsOverlap_best), path(TS_partition_PRS_summary), path(TS_partition_PRS_prsice), path (TS_partition_PRS_best)
    // each path(GW_LD_blocks)

    output:
    // tuple val(meta), path("*_modified_GWAS.csv"),       emit: modified_GWAS
    path("*.png")
    path("versions.yml")         ,                      emit: versions

    script:
    """
    R_PRS_PPV_plotting.R ${meta} ${original_PRS_prsice} ${original_NoEPsOverlap_prsice} ${TS_partition_PRS_prsice} ${TS_EPs_prsice} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
    