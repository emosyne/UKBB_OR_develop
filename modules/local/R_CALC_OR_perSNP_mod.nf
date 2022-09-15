process R_CALC_OR_perSNP {
    debug true
    container 'emosyne/r_docker:1.7'
    stageInMode 'copy'
    label 'process_high_memory'
    tag "$meta"
    maxRetries 3

    input:
    tuple val(meta), path (bim), path (bedfilepath),  path (fam), path (log), path(SNVs_hg19), path(pheno)

    output:
    path("*_odds_ratios_UKBB.csv")     ,     emit: odds_ratios_UKBB


    script:
    """
    R_CALC_OR_perSNP.R ${bedfilepath} ${pheno}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}