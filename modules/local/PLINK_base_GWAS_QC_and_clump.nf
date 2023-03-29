process PLINK_base_GWAS_QC_and_clump {
    // tag "$cohort"
    // debug true
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple val(condition), path(full_GWAS), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)
    

    output:
    tuple val(condition), path ("*GWAS_QC_nodups.tsv.gz"),                 emit: GWAS_QC_noClump
    tuple val(condition), path ("*GWAS_QC_nodups_clump.clumped"),          emit: clumped_SNPs
    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    #remove SNPs with MAF < 0.01 (in this case main allele freq <0.99) and remove duplicated SNPs
    zcat ${full_GWAS} | awk 'NR==1 || (\$6 < 0.99) {print}' | awk '!seen[\$1]++' | sed -E 's/,/\t/g' | gzip > ${condition}_GWAS_QC_nodups.tsv.gz

    plink  \\
       --clump ${condition}_GWAS_QC_nodups.tsv.gz \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --maf 0.01 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out ${condition}_GWAS_QC_nodups_clump  \\
       --threads $task.cpus \\
       --memory $mem_mb

    """
}

