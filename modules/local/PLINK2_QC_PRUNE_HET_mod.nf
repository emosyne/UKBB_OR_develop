process PLINK2_QC_PRUNE_HET {
    tag "$meta"
    // debug true
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input: 
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log)
    
    

    output:
    tuple val(meta), path ("*.prune.in"), path ("*.het"), emit: pruned_variants_het
    // path("*.log")
    path "versions2.yml"           , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    //def prefix = task.ext.prefix 
    // if( "$PLINKbgenfiles" == "${prefix}.bgen" ) error "Input and output names are the same, use \"task.ext.prefix\" in modules.config to disambiguate!"
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
      
    

    cat <<-END_VERSIONS > versions2.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

