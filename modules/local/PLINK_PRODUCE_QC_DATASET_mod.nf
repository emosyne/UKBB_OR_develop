process PLINK_PRODUCE_QC_DATASET {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'


    input:
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log), path (het_valid), path(a1_bim), path(mismatch)
    

    output:
    tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    path "versions1.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    """
    plink \\
        --bfile ${bedfilepath.baseName} \\
        --make-bed \\
        --threads $task.cpus \\
        --out ${bedfilepath.baseName}_QC


    cat <<-END_VERSIONS > versions1.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    
    """
}
// --a1-allele ${a1_bim} \\
// --keep ${het_valid} \\
// --exclude ${mismatch} \\