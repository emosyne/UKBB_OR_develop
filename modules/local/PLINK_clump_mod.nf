process PLINK_clump {
    // debug true
    tag "${ENH_list}_${EPWAS_model}"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input:
    // [SCZ, SCZ_ALLCHR_SCZ_QC.bed, SCZ_ALLCHR_SCZ_QC.bim, SCZ_ALLCHR_SCZ_QC.fam, 
        // Neural_significant_enh, SCZ_ADD_Neural_significant_enh_noclump_EPWAS.tsv.gz, SCZ_ADD_Neural_significant_enh_noclump_residual_GWAS_compartment.tsv.gz, 
        // ADD, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
    tuple val(condition), path(bed_QC),  path(bim_QC), path(fam_QC), \
        val(ENH_list), path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment), \
        val(EPWAS_model), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment), path("*clumped_SNPs.clumped"), \
        val(condition), val(EPWAS_model),   emit: clumped_SNPs_and_noclump_lists
    path("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink  \\
       --clump ${noclump_EPWAS},${noclump_residual_GWAS_compartment} \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out ${condition}_${ENH_list}_clumped_SNPs  \\
       --threads $task.cpus \\
       --memory $mem_mb

    
    """
}
