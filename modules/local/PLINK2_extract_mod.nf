process PLINK2_EXTRACT {
    // debug true
    label 'process_high_memory'
    // maxForks 4
    tag "$meta"
    // errorStrategy = 'finish' 
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input: 
    tuple val(meta), path(GWAS_ENH_MERGE_BED_hg19), path(ENH_bed), path(pheno), path(ENH_csv), path(chr_bgenfile), path(chr_samplefile)
    

    output:
    tuple val(meta), path("*.bim"), path("*.bed"), path ("*.fam"), path("*.log"), emit: SNPextracted_by_chromosome
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bgen ${chr_bgenfile} ref-first \\
        --sample ${chr_samplefile} \\
        --pheno $pheno \\
        --rm-dup force-first --make-bed \\
        --extract bed1 ${GWAS_ENH_MERGE_BED_hg19} \\
        --out ${GWAS_ENH_MERGE_BED_hg19.simpleName}_${chr_bgenfile.simpleName}


    """
}
