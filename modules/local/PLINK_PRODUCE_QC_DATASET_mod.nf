process PLINK_PRODUCE_QC_DATASET {
    // debug true
    tag "$condition"
    label 'process_high_memory'
    container 'emosyne/plink2:1.23'
    cache 'lenient'


    input:
    tuple val(condition), path(bed), path(bim), path(fam), path(HCM_GWAS_QC), path (het_valid), path(a1_bim), path(mismatch)
    

    output:
    tuple val(condition), path ("*_QC.bed"), path ("*_QC.bim"), path ("*_QC.fam"), path(HCM_GWAS_QC),            emit: target_QC
    path ("*.log")



    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
        plink \\
        --bfile ${bed.simpleName} \\
        --make-bed \\
        --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \\
        --a1-allele ${a1_bim} \\
        --keep ${het_valid} \\
        --exclude ${mismatch} \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ${bed.simpleName}_${condition}_QC

    
    """
}
// --a1-allele ${a1_bim} \\
// --keep ${het_valid} \\
// --exclude ${mismatch} \\