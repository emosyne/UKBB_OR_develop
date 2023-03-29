process PLINK_PRODUCE_QC_DATASET {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'


    input:
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log), path (het_valid), path(a1_bim), path(mismatch)
    

    output:
    tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam") , emit: all_chromosomes_QC




    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink \\
        --bfile ${bedfilepath.baseName} \\
        --make-bed \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ${bedfilepath.baseName}_QC

    
    """
}
// --a1-allele ${a1_bim} \\
// --keep ${het_valid} \\
// --exclude ${mismatch} \\