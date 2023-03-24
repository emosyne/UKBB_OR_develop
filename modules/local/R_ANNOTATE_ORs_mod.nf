process R_ANNOTATE_ORs {
    // debug true
    container 'emosyne/r_docker:1.97'
    stageInMode 'copy'
    label 'process_high_memory'
    tag "$meta"
    cache 'lenient'

    input:
    tuple val(meta), path(associations_recessive),  path(associations_standardGLM),  path(processed_SNPlists_hg19), path(full_GWAS)
    each path(hg38ToHg19_chain)
    each path(GW_LD_blocks)

    output:
    tuple val(meta), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    path "figs/*_UKBB.png"                      ,       emit: ORfigure
    path("versions.yml")                        ,       emit: versions

    script:
    """
    R_annotate_ORs.R ${meta} ${associations_recessive} ${associations_standardGLM} ${processed_SNPlists_hg19} ${full_GWAS} ${hg38ToHg19_chain} ${GW_LD_blocks}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
        R_GenomicRanges: \$(Rscript -e 'packageVersion("GenomicRanges")' | awk '{print \$2}' | tr -d "‘’")
        R_rtracklayer: \$(Rscript -e 'packageVersion("rtracklayer")' | awk '{print \$2}' | tr -d "‘’")
        ggnewscale: \$(Rscript -e 'packageVersion("ggnewscale")' | awk '{print \$2}' | tr -d "‘’")
        biomaRt: \$(Rscript -e 'packageVersion("biomaRt")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
    