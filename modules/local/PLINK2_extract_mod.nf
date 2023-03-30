process PLINK2_EXTRACT {
    // debug true
    label 'process_high_memory'
    // maxForks 4
    // tag "$meta"
    // errorStrategy = 'finish' 
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input: 
    // [SCZ_clumped_GWAS_SNPs_plus_those_in_bed_files.bed, ukb22828_c13_b0_v3.bgen, ukb22828_c13_b0_v3.sample]
    tuple path(extracted_GWAS_SNPs_bed), path(chr_bgenfile), path(chr_samplefile)
    

    output:
    tuple val("SCZ"), path("*.bim"), path("*.bed"), path ("*.fam"), path("*.log"), emit: SNPextracted_by_chromosome
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bgen ${chr_bgenfile} ref-first \\
        --sample ${chr_samplefile} \\
        --rm-dup force-first --make-bed \\
        --extract bed1 ${extracted_GWAS_SNPs_bed} \\
        --out extracted_GWAS_SNPs_${chr_bgenfile.simpleName}


    """
}
