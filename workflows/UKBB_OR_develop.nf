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

dx_UKBB_pheno =    Channel.fromPath("../private_input_files/biobank/SCZ.pheno", checkIfExists: true).map{it-> ["SCZ", it]}




// #### UKBB PLINK input files ####
genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/*c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }



UKBBethinicityRelatedness = Channel.fromPath( '../private_input_files/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
UKBB_covariates = Channel.fromPath( '../private_input_files/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)
UKBB_rs34380086_cases = Channel.fromPath( '../private_input_files/biobank/rs34380086_cases.pheno', checkIfExists: true)


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
annotations = Channel.fromPath( "../private_input_files/ES_multipliers/2023-01-18_2023-02-17_NEURAL_ENH_EXP_significant_plus_100_noOverlap_HCMformat.csv.gz", checkIfExists: true)


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

    PLINK2_ASSOC_GLM(
        PLINK_PRODUCE_QC_DATASET.out.target_QC
            .combine(enhancer_lists_bed_files.map{it -> it[1]}),// ENH hg19 bed file
        UKBB_covariates
        )
    
    // PLINK2_ASSOC_GLM.out.associations // ORs
    //     .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
    //     .view() 
    
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

    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    // R_prepare_lists_for_clump.out.lists_before_clump
    //     .combine(LD_reference, by: [0])
    //     .map{it.flatten()}
    //     .view()
    
    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference, by: [0])
            .map{it.flatten()}
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists.view()
    

    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists,//.map { [it, "1"].flatten() }, //######################## multiplier can be set here ########################
        annotations
    )

    
    R_split_lists.out.partitioned 
        .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
        .combine(UKBB_covariates)
        .combine(LD_reference)
        .map { [it, "0.5"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_05}
    R_split_lists.out.partitioned 
        .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
        .combine(UKBB_covariates)
        .combine(LD_reference)
        .map { [it, "0.05"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_005}
    
    combined_splitlists_bedfile_QCeddata_LDdata = combined_splitlists_bedfile_QCeddata_LDdata_05.mix(combined_splitlists_bedfile_QCeddata_LDdata_005)
    // combined_splitlists_bedfile_QCeddata_LDdata.first().view()
    
    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### SET NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_EPWAS_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_GWAS_PRS)
            .map { [it, "enh_ES", "enh_TS_tpm"].flatten() }


    // PRS_results.view()

    R_final_plot (
        PRS_results
    )


}



