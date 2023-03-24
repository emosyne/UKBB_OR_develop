process PLINK_MERGE {
    // debug true
    tag "$meta"
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input:
    tuple val(meta), val(bedfiles)
    each path (PLINKethinicityRelatedness)
    

    output:
    tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_extracted
    path "chr_file_list.txt"
    path "versions1.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    // def mem_mb = task.memory.toMega()

    // if( "$bed" == "${prefix}.bed" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # IMPORT LIST OF FILES INTO chr_file_list
    echo ${bedfiles} | \\
        tr -d ',[]' | \\
        tr ' ' '\\n' | \\
        sed 's/.bed//g' > \\
        chr_file_list.txt
    # take first chrom file
    first_chr=\$(head -n 1 chr_file_list.txt)
    echo \$first_chr
    # REMOVE FIRST CHR FROM FILE
    sed -i '1d' chr_file_list.txt 
    
    plink --bfile \$first_chr \\
        --merge-list chr_file_list.txt \\
        --make-bed \\
        --remove $PLINKethinicityRelatedness \\
        --maf 0.01 \\
        --hwe 1e-6 \\
        --geno 0.01 \\
        --mac 100 \\
        --threads $task.cpus \\
        --out ${prefix}_${meta} 


    cat <<-END_VERSIONS > versions1.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    
    """
}
