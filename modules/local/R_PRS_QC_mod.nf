process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.7'
    stageInMode 'copy'
    label 'process_low'
    tag "$meta"
    cache 'lenient'


    input: 
    tuple val(meta), path (prune), path (het), path (bedfilepath), path (bim), path (fam), path (log), path(GWAS)
    
    

    output:
    tuple val(meta), path ("*_het_valid_out.sample"), path("*_a1_bim"), path("*_mismatching_SNPs"),  emit: QC_het_a1_mismatch
    path "versions.yml"           , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    R_PRS_QC.R ${meta} ${het} ${bim} ${GWAS} 
    
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