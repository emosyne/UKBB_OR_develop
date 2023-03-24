process PLINK_MERGE {
    // debug true
    tag "$meta"
    label 'vlarge'
    container 'emosyne/plink2:1.23'
    cache 'lenient'

    input:
    tuple val(meta), val(bedfiles)
    each path (PLINKethinicityRelatedness)
    

    output:
    tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_extracted
    path "chr_file_list.txt"



    script:
    def mem_mb = (task.memory * 0.95).toMega()
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
        --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ${meta}_mergedfile

    
    """
}
