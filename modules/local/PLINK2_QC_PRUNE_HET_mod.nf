process PLINK2_QC_PRUNE_HET {
    tag "$condition"
    // debug true
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input: 
    tuple val(condition), path (bedfilepath), path (bim), path (fam), path (log)
    
    

    output:
    
    tuple val(condition), path (bedfilepath), path (bim), path (fam), path ("*.prune.in"), path ("*.het"),  emit: pruned_variants_het
    // path("*.log")
    

    
    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink \\
      --bfile ${bedfilepath.baseName} \\
      --indep-pairwise 200 50 0.25 \\
      --out ${bedfilepath.baseName} \\
      --threads $task.cpus

    plink \\
      --bfile ${bedfilepath.baseName} \\
      --extract ${bedfilepath.baseName}.prune.in \\
      --het \\
      --out ${bedfilepath.baseName} \\
      --threads $task.cpus
      
    

    """
}

