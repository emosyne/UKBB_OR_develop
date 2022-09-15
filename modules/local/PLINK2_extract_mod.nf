process PLINK2_EXTRACT {
    // debug true
    label 'process_high_memory'
    // maxForks 4
    tag "$meta"
    errorStrategy = 'finish' 
    container 'emosyne/plink2:1.2'
    cache 'lenient'

    input: 
    tuple val(meta), path(GWAS_ENH_MERGE_BED_hg19), path(ENH_bed), path(pheno), path(ENH_csv), path(chr_bgenfile), path(chr_samplefile)
    

    output:
    tuple val(meta), path("*.bim"), path("*.bed"), path ("*.fam"), path("*.log"), emit: SNPextracted_by_chromosome
    path "versions.yml"                , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix 
    // if( "$PLINKbgenfiles" == "${prefix}.bgen" ) error "Input and output names are the same, use \"task.ext.prefix\" in modules.config to disambiguate!"
    def mem_mb = task.memory.toMega()
    """ 
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bgen ${chr_bgenfile} ref-first \\
        --sample ${chr_samplefile} \\
        --pheno $pheno \\
        --rm-dup force-first --make-bed \\
        --extract bed1 ${GWAS_ENH_MERGE_BED_hg19} \\
        --out ${prefix}_${GWAS_ENH_MERGE_BED_hg19.simpleName}_${chr_bgenfile.simpleName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
