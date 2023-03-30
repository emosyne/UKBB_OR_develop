include { PLINK2_EXTRACT }                  from '../modules/local/PLINK2_extract_mod.nf'
include { PLINK_MERGE }                     from '../modules/local/plink_mergechromfiles_mod.nf'
include { PLINK2_ASSOC_GLM }                from '../modules/local/PLINK2_ASSOC_GLM_mod.nf'
include { R_ANNOTATE_ORs }                  from '../modules/local/R_ANNOTATE_ORs_mod.nf'
include { PLINK2_QC_PRUNE_HET }             from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                        from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }        from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { PLINK_clump }                     from '../modules/local/PLINK_clump_mod.nf'
include { PLINK_base_GWAS_QC_and_clump }    from '../modules/local/PLINK_base_GWAS_QC_and_clump.nf'
include { R_extract_GWAS_SNPs_into_bed }    from '../modules/local/R_extract_GWAS_SNPs_into_bed.nf'
include { R_prepare_lists_for_clump }       from '../modules/local/R_prepare_lists_for_clump.nf'
include { R_split_lists }                   from '../modules/local/R_split_lists.nf'
include { PRSice_calculate_PRS_split_partitions }   from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include { R_final_plot }                    from '../modules/local/R_final_plot_mod.nf'





full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz", checkIfExists: true) 
    .map{it-> ["SCZ", it]}

dx_UKBB_pheno =    Channel.fromPath("./input/biobank/SCZ.pheno", checkIfExists: true).map{it-> ["SCZ", it]}



// // chain file
// hg38ToHg19_chain = Channel
//     .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)
// //LD blocks 1000 genomes
// GW_LD_blocks = Channel
//     .fromPath( "./input/LD/EUR_phase3_autosomes_hg19.blocks.det.gz", checkIfExists: true)

// #### UKBB PLINK input files ####
genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/*c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }



UKBBethinicityRelatedness = Channel.fromPath( './input/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
UKBB_covariates = Channel.fromPath( './input/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)
UKBB_rs34380086_cases = Channel.fromPath( './input/biobank/rs34380086_cases.pheno', checkIfExists: true)


//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect().map{it-> ["SCZ", it]}

// ################################ INTERNAL VALIDATION INPUTS ################################

full_GWAS_HCMformat = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz", checkIfExists: true) 
    .map{it-> ["SCZ", it]}

enhancer_lists_bed_files = 
    Channel.from("Neural_significant_enh")
            .map { ENH_list -> ["${ENH_list}", 
                file("./input/validation/enhancer_files/${ENH_list}.bed", checkIfExists: true)]
            } 

workflow UKBB_OR_develop {
    
     // ################################ EPWAS development ################################

    PLINK_base_GWAS_QC_and_clump (
        full_GWAS_HCMformat
            .join(LD_reference)
            .map{it.flatten()}
    )
    
    

    R_extract_GWAS_SNPs_into_bed ( 
        // THIS MODULE IMPORTS 
        // GWAS (hg19), and selects all SNPs in input bed files and all GWAS clumped SNPs and outputs a bed file
        enhancer_lists_bed_files.map{it -> it[1]}.collect(),
        PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump
            .join(PLINK_base_GWAS_QC_and_clump.out.clumped_SNPs)
        )
    
    chromosomes_by_condition_plus_SNPs = 
        R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
            .combine(genotype_chr_files) //The combine operator combines (cartesian product) the items emitted by two channels
            
    // chromosomes_by_condition_plus_SNPs.view()

    PLINK2_EXTRACT ( 
        // extract genotypes at bed file locations
        chromosomes_by_condition_plus_SNPs

        //out tuple val(meta), path("*.bim"), path("*.bed"), path ("*.fam"), path("*.log"), emit: SNPextracted_by_chromosome
        )
    
    // PLINK2_EXTRACT.out.SNPextracted_by_chromosome.view()
    // [SCZ, GWAS_SCZ_SNV_merge_hg19_ukb22828_c17_b0_v3.bim, GWAS_SCZ_SNV_merge_hg19_ukb22828_c17_b0_v3.bed, GWAS_SCZ_SNV_merge_hg19_ukb22828_c17_b0_v3.fam, GWAS_SCZ_SNV_merge_hg19_ukb22828_c17_b0_v3.log]

    PLINK2_EXTRACT.out.SNPextracted_by_chromosome
        .branch{
            SCZ: it =~ /SCZ/
            HCM: it =~ /HCM/
        }
        .set{ SNPextracted_by_chromosome_byMeta }
        
    
    //split channel by meta, collect all bed files per chr, and generate tuple
    bedfiles = Channel
        .empty()
        .mix(
            SNPextracted_by_chromosome_byMeta.SCZ
                .map{it-> it[2]}
                .collect()
                .map{it-> ["SCZ", it]}
                // .view()
                )
        .mix(
            SNPextracted_by_chromosome_byMeta.HCM
                .map{it-> it[2]}
                .collect()
                .map{it-> ["HCM", it]}
                // .view()
        )
    // bedfiles.view()
    
    

    PLINK_MERGE(
        bedfiles, 
        UKBBethinicityRelatedness,
        dx_UKBB_pheno.map{it-> [it[1]]}

        //out tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_extracted
        )


    // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // produce prune.in and het files
    PLINK2_QC_PRUNE_HET (
        PLINK_MERGE.out.all_chromosomes_extracted
    )
    
    // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //         .join(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump)
    //         .view()

    // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // produce QC_het_a1_mismatch, 
    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
             .join(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump)
    )
    // R_PRS_QC.out.QC_het_a1_mismatch.view()
    //[SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_ALLCHR.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_ALLCHR.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_ALLCHR.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_het_valid_out_vs_GWAS.sample, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_a1_cohort_bim_vs_GWAS, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/da/7f6f852e60e9d996747fda3f1d5aa1/SCZ_mismatching_SNPs_vs_GWAS]

    // TARGET QC 3:  
    // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // also remove mismatching SNPs
    // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
    PLINK_PRODUCE_QC_DATASET ( 
        R_PRS_QC.out.QC_het_a1_mismatch
    )

    // PLINK_PRODUCE_QC_DATASET.out.target_QC
    //         .combine(enhancer_lists_bed_files.map{it -> it[1]})
    //         .view()
    // [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/input/validation/enhancer_files/Neural_significant_enh.bed]

    PLINK2_ASSOC_GLM(
        PLINK_PRODUCE_QC_DATASET.out.target_QC
            .combine(enhancer_lists_bed_files.map{it -> it[1]}),// ENH hg19 bed file
        UKBB_covariates
        )
    
    // PLINK2_ASSOC_GLM.out.associations // ORs
    //     .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
    //     .view() 
    // /[SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/0b/5bf167f5ec07a5519081d3a4d1847a/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/0b/5bf167f5ec07a5519081d3a4d1847a/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/0b/5bf167f5ec07a5519081d3a4d1847a/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_dominant.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/0b/5bf167f5ec07a5519081d3a4d1847a/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.frq, /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz]
    
    R_ANNOTATE_ORs(
        // annotate ORs from previous step with GWAS results and other info,
        //produce OR plots
        PLINK2_ASSOC_GLM.out.associations // ORs
            .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
        
        // out: tuple val(meta), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    )
    
// ################################ INTERNAL VALIDATION ################################


    enhancer_EPWAS_files =
        R_ANNOTATE_ORs.out.EPWAS_HCM_format_REC
            .concat(R_ANNOTATE_ORs.out.EPWAS_HCM_format_ADD)
            .concat(R_ANNOTATE_ORs.out.EPWAS_HCM_format_DOM)
    
    // enhancer_EPWAS_files.view()
    
    
    PLINK_PRODUCE_QC_DATASET.out.target_QC
        .combine(enhancer_lists_bed_files)
        .combine(enhancer_EPWAS_files)
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    
    // cohort_GWAS_enh_list.view()
    // [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_GWAS_QC_nodups.tsv.gz, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/input/validation/enhancer_files/Neural_significant_enh.bed, REC, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/29/55fbfe539be2822470efedd01e1707/UKBB_ENH_associations_REC.tsv.gz]
    // [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_GWAS_QC_nodups.tsv.gz, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/input/validation/enhancer_files/Neural_significant_enh.bed, ADD, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/29/55fbfe539be2822470efedd01e1707/UKBB_ENH_associations_ADD.tsv.gz]
    // [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/d7/d302e50863cd085dbe4f45b217ae25/SCZ_GWAS_QC_nodups.tsv.gz, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/input/validation/enhancer_files/Neural_significant_enh.bed, DOM, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/29/55fbfe539be2822470efedd01e1707/UKBB_ENH_associations_DOM.tsv.gz]

    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
    R_prepare_lists_for_clump.out.lists_before_clump
        .combine(LD_reference)
        .view()

    
//     PLINK_clump (
//         //CLUMPING of enhancer-based SNP compartments together 
//         R_prepare_lists_for_clump.out.lists_before_clump
//             .combine(LD_reference)
//     )
//     // PLINK_clump.out.clumped_SNPs_and_noclump_lists.view()
//     // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/SCZ_REC_Neural_significant_enh_noclump_EPWAS.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/SCZ_REC_Neural_significant_enh_noclump_residual_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/82/3cd69072a998abd4ba99f1bae51356/SCZ_Neural_significant_enh_clumped_SNPs.clumped, SCZ, REC]
    

//     R_split_lists (
//         // first annotate SNPs with ES of relevant E-P - for ENH SNP list
//         // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
//         // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
//         // output separate lists to calculate split PRSs and also merged one
//         PLINK_clump.out.clumped_SNPs_and_noclump_lists,//.map { [it, "1"].flatten() }, //######################## multiplier can be set here ########################
//         annotations
//     )

    
//     R_split_lists.out.partitioned 
//         .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
//         .combine(UKBB_covariates)
//         .combine(LD_reference)
//         .map { [it, "0.5"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
//         .set{combined_splitlists_bedfile_QCeddata_LDdata_05}
//     R_split_lists.out.partitioned 
//         .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
//         .combine(UKBB_covariates)
//         .combine(LD_reference)
//         .map { [it, "0.05"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
//         .set{combined_splitlists_bedfile_QCeddata_LDdata_005}
    
//     combined_splitlists_bedfile_QCeddata_LDdata = combined_splitlists_bedfile_QCeddata_LDdata_05.mix(combined_splitlists_bedfile_QCeddata_LDdata_005)
//     // combined_splitlists_bedfile_QCeddata_LDdata.view()
//     // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, Neural_significant_enh, Neural_significant_enh_REC_SCZ_X_1_clumped_EPWAS.tsv.gz, Neural_significant_enh_REC_SCZ_clumped_residual_GWAS_compartment.tsv.gz, 1, SCZ, REC, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b0/fcec9ad982f3a82c3b2bef76691d0a/SCZ_clumped_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam, 0.05]
//     // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, Neural_significant_enh, Neural_significant_enh_REC_SCZ_X_1_clumped_EPWAS.tsv.gz, Neural_significant_enh_REC_SCZ_clumped_residual_GWAS_compartment.tsv.gz, 1, SCZ, REC, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b0/fcec9ad982f3a82c3b2bef76691d0a/SCZ_clumped_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam, 0.5]

    
//     PRSice_calculate_PRS_split_partitions(
//         combined_splitlists_bedfile_QCeddata_LDdata
//     )
    
//     // ########################################### SET NAMES OF MULTIPLIERS ###########################################
//     PRS_results = 
//         PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
//             .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
//             .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
//             .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_GWAS_PRS)
//             .map { [it, "enh_ES", "enh_TS_tpm"].flatten() }


//     PRS_results.view()
    

// //     R_final_plot (
// //         PRS_results
// //     )


}



