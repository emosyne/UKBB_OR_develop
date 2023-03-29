include { GENERATESNPLISTS }                from '../modules/local/R_GENERATESNPLISTS_mod.nf'
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
// include {  }     from '../modules/local/'
// include {  }     from '../modules/local/'
// include {  }     from '../modules/local/'
// include {  }     from '../modules/local/'
// include {  }     from '../modules/local/'





enhancer_plus_GWAS_coords = Channel.from("SCZ") 
    .map { condition -> ["${condition}", 
           file("./input/textfiles/ENH_${condition}_hg38.csv.gz", checkIfExists: true), 
           file("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz", checkIfExists: true), 
           file("./input/biobank/${condition}.pheno", checkIfExists: true)] } 
        //    .view()



full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz", checkIfExists: true) 
    .map{it-> ["SCZ", it]}
full_GWAS_HCMformat = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz", checkIfExists: true) 
    .map{it-> ["SCZ", it]}


// chain file
hg38ToHg19_chain = Channel
    .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)
//LD blocks 1000 genomes
GW_LD_blocks = Channel
    .fromPath( "./input/LD/EUR_phase3_autosomes_hg19.blocks.det.gz", checkIfExists: true)

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



workflow UKBB_OR_develop {
    // ################################ EPWAS development ################################
    GENERATESNPLISTS( 
        // THIS MODULE IMPORTS E-PS LIST (hg38) AND GWAS results (hg19), converts them to hg 19 and merges them
        // outputting bed files
        enhancer_plus_GWAS_coords, 
        hg38ToHg19_chain

        //out tuple val(meta), path("GWAS_*_hg19.bed"), path("ENH_*_hg19.bed"), path(pheno), path("ENH_*_hg19.csv"), emit: processed_ENH_SNP_lists_hg19
        )
    // GENERATESNPLISTS.out.processed_ENH_SNP_lists_hg19.view()
    
    chromosomes_by_condition_plus_SNPs = 
        GENERATESNPLISTS.out.processed_ENH_SNP_lists_hg19
            .combine(genotype_chr_files) //The combine operator combines (cartesian product) the items emitted by two channels
            // .view()
        


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
        // merge all bed files into one:
        bedfiles, 
        UKBBethinicityRelatedness

        //out tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_extracted
        )
    
    // PLINK_MERGE.out.all_chromosomes_extracted.view() [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/8c/41dfb0b869376a5314f9b4a404da3e/SCZ_mergedfile.bed, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/8c/41dfb0b869376a5314f9b4a404da3e/SCZ_mergedfile.bim, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/8c/41dfb0b869376a5314f9b4a404da3e/SCZ_mergedfile.fam, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/8c/41dfb0b869376a5314f9b4a404da3e/SCZ_mergedfile.log]
    // PLINK_MERGE.out.chrfilelist.view()


    PLINK2_ASSOC_GLM(

        PLINK_MERGE.out.all_chromosomes_extracted
            //join will join all_chromosomes_extracted with the SNP list output from step 1 by condition (meta)
            .join(GENERATESNPLISTS.out.processed_ENH_SNP_lists_hg19.map{it->[it[0],it[2]]}, by: [0])//join ENH hg19 bed file
            .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0]), // also join pheno file
        UKBB_covariates,
        Channel.fromPath("./input/textfiles/Neural_significant_enh.bed", checkIfExists: true) 

        //out tuple val(meta), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid"), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid"), emit: associations
        )
        //  PLINK2_ASSOC_GLM.out.associations // ORs
        //     .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
        //     .view() 
        // [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/6f/d017d89d2b671cc0ff910a7ced8502/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/6f/d017d89d2b671cc0ff910a7ced8502/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/6f/d017d89d2b671cc0ff910a7ced8502/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_dominant.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/6f/d017d89d2b671cc0ff910a7ced8502/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.frq, /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz]
    
    R_ANNOTATE_ORs(
        // annotate ORs from previous step with GWAS results and other info,
        //produce OR plots
        PLINK2_ASSOC_GLM.out.associations // ORs
            .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
        
        // out: tuple val(meta), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    )
    
// ################################ INTERNAL VALIDATION ################################

    PLINK_base_GWAS_QC_and_clump (
        full_GWAS_HCMformat
            .combine(LD_reference)
            .map { it.flatten() }
    )
    
    

    // R_extract_GWAS_SNPs_into_bed ( 
    //     // THIS MODULE IMPORTS 
    //     // GWAS (hg19), and selects all SNPs in input bed files and all GWAS clumped SNPs and outputs a bed file
    //     enhancer_lists_bed_files.map{it -> it[1]}.mix(Channel.fromPath("./input/EPWAS/EP_WAS.bed")).collect(),
    //     PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump
    //         .combine(PLINK_base_GWAS_QC_and_clump.out.clumped_SNPs)
    //         .map { [it, condition].flatten() }
        
    //     )
//     // R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
//     //     .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
//     //     .view()
//     chromosomes_by_condition_plus_SNPs = 
//         // PGC_GWAS_plus_allEPlists_SNPs_hg19
//         R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
//             .combine(genotype_chr_files) //The combine operator combines (cartesian product) the items emitted by two channels
            
        
//     // chromosomes_by_condition_plus_SNPs.view()

//     // GENERATE UKBB UNIQUE FILE
//     PLINK2_EXTRACT ( 
//         // extract genotypes at bed file locations
//         chromosomes_by_condition_plus_SNPs

//         //out tuple val(meta), path("*.bim"), path("*.bed"), path ("*.fam"),  emit: SNPextracted_by_chromosome
//         )
    
    

//     PLINK_MERGE( // SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL
//         // merge all bed files into one:
//         PLINK2_EXTRACT.out.SNPextracted_by_chromosome.collect(),
//         UKBBethinicityRelatedness,
//         dx_UKBB_pheno
//         //out tuplepath ("*.bed"), path ("*.bim"), path ("*.fam"),  emit: all_chromosomes_extracted
//         )
//     // PLINK_MERGE.out.all_chromosomes_extracted.view()

//     // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
//     // produce prune.in and het files
//     PLINK2_QC_PRUNE_HET (
//         PLINK_MERGE.out.all_chromosomes_extracted
//     )
    
//     // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
//     //         .combine(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC)
//     //         .view()
//     // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.prune.in, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.het, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/0b/036c27bada52d3b859916ac5896889/GWAS_QC.gz]

//     // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
//     // produce QC_het_a1_mismatch, 
//     R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
//         PLINK2_QC_PRUNE_HET.out.pruned_variants_het
//             .combine(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump)
//             .map { [it, condition].flatten() }
//     )
//     // R_PRS_QC.out.QC_het_a1_mismatch.view()
//     //[/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_het_valid_out_vs_HCM_GWAS.sample, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_a1_cohort_bim_vs_HCM_GWAS, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_mismatching_SNPs_vs_HCM_GWAS, SCZ]

//     // TARGET QC 3:  
//     // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
//     // also remove mismatching SNPs
//     // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
//     PLINK_PRODUCE_QC_DATASET ( //   SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL     --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 \\

//         R_PRS_QC.out.QC_het_a1_mismatch
//     )

//     // PLINK_PRODUCE_QC_DATASET.out.target_QC.view()
//     //[GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz]
    
    
//     PLINK_PRODUCE_QC_DATASET.out.target_QC
//         .combine(enhancer_lists_bed_files)
//         .combine(enhancer_EPWAS_files)
//         .map { it.flatten() }
//         .set{cohort_GWAS_enh_list}
    
//     // cohort_GWAS_enh_list.view()
//     // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, SCZ, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/enh_bedfiles/Neural_significant_enh.bed, REC, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/EPWAS/UKBB_ENH_associations_REC.tsv.gz]
//     // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, SCZ, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/enh_bedfiles/Neural_significant_enh.bed, DOM, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/EPWAS/UKBB_ENH_associations_DOM.tsv.gz]
//     // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, SCZ, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/enh_bedfiles/Neural_significant_enh.bed, ADD, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/EPWAS/UKBB_ENH_associations_ADD.tsv.gz]
    
//     // BASE subsetting
//     R_prepare_lists_for_clump (
//         // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
//         // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
//         cohort_GWAS_enh_list
//     )
    
    
//     //    R_prepare_lists_for_clump.out.lists_before_clump
//     //         .combine(LD_reference)
//     //         .view()
//     // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8c/5832511f6ff56c817623ff8511e5a4/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8c/5832511f6ff56c817623ff8511e5a4/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8c/5832511f6ff56c817623ff8511e5a4/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8c/5832511f6ff56c817623ff8511e5a4/SCZ_REC_Neural_significant_enh_noclump_EPWAS.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8c/5832511f6ff56c817623ff8511e5a4/SCZ_REC_Neural_significant_enh_noclump_residual_GWAS_compartment.tsv.gz, SCZ, REC, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
//     // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/ce/148eaf9bc44e293f60bb7e97ca30b7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/ce/148eaf9bc44e293f60bb7e97ca30b7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/ce/148eaf9bc44e293f60bb7e97ca30b7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/ce/148eaf9bc44e293f60bb7e97ca30b7/SCZ_DOM_Neural_significant_enh_noclump_EPWAS.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/ce/148eaf9bc44e293f60bb7e97ca30b7/SCZ_DOM_Neural_significant_enh_noclump_residual_GWAS_compartment.tsv.gz, SCZ, DOM, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
    
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



