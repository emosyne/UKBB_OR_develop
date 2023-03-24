process PLINK2_SUBSET {
    // debug true
    label 'process_high_memory'
    tag "$meta"
    cache 'lenient'
    container 'emosyne/plink2:1.23'

    input: 
    tuple val(meta), path (bedfilepath), path (bim), path (fam), path (log), path(SNVs_hg19)

    output:
    tuple val(meta), path("*_subsetted_SNVs_hg19.bim"), path("*_subsetted_SNVs_hg19.bed"), path("*_subsetted_SNVs_hg19.fam"), path("*_subsetted_SNVs_hg19.log"), emit: subsetted_genotypes
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix 
    // if( "$PLINKbgenfiles" == "${prefix}.bgen" ) error "Input and output names are the same, use \"task.ext.prefix\" in modules.config to disambiguate!"
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile ${bedfilepath.baseName} \\
        --extract bed1 ${SNVs_hg19} \\
        --out ${bedfilepath.simpleName}_subsetted_SNVs_hg19 \\
        --make-bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
