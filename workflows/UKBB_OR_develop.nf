include { GENERATESNPLISTS }                from '../modules/local/R_GENERATESNPLISTS_mod.nf'
include { PLINK2_EXTRACT }                  from '../modules/local/PLINK2_extract_mod.nf'
include { PLINK_MERGE }                     from '../modules/local/plink_mergechromfiles_mod.nf'
include { PLINK2_ASSOC_GLM }                from '../modules/local/PLINK2_ASSOC_GLM_mod.nf'
include { R_ANNOTATE_ORs }                  from '../modules/local/R_ANNOTATE_ORs_mod.nf'
include { PLINK2_QC_PRUNE_HET }             from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                        from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }        from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_PREPARE_MODIF_PRS }             from '../modules/local/R_prepare_modified_PRS_mod.nf'
include { PLINK_clump }                     from '../modules/local/PLINK_clump_mod.nf'
include { R_PREPARE_MODIF_PRS_2_LISTS }     from '../modules/local/R_prepare_modified_PRS_2_mod.nf'
include { PRSice_calculate_PRS_original }   from '../modules/local/PRSice_calculate_PRS_mod.nf'
include { PRSice_calculate_PRS_substituted }        from '../modules/local/PRSice_calculate_PRS_2_mod.nf'
include { PRSice_calculate_PRS_TS_partition }       from '../modules/local/PRSice_calculate_PRS_3_mod.nf'
include { PRSice_calculate_4_partition }    from '../modules/local/PRSice_calculate_PRS_4_partition.nf'
include { R_PRS_PPV_plotting }              from '../modules/local/R_PRS_PPV_plotting_mod.nf'
include { R_calculate_compound_PRS_fit_and_plot }   from '../modules/local/R_calculate_compound_PRS_fit_and_plot_mod.nf'
include { R_plot_GO }                       from '../modules/local/R_plot_GO.nf'




enhancer_plus_GWAS_coords = Channel.from("SCZ") //,"HCM"
    .map { condition -> ["${condition}", 
           file("./input/textfiles/ENH_${condition}_hg38.csv.gz", checkIfExists: true), 
           file("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz", checkIfExists: true), 
           file("./input/biobank/${condition}.pheno", checkIfExists: true)] } 
        //    .view()



//  SCHIZO and neural lists ##############
full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_overInfo.8_OR.tsv.gz", checkIfExists: true) 
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
// LD_reference = Channel.from("bed","bim","fam") 
//     .map { ext -> ["${ext}", 
//             file("/mnt/storage/emanuele/LDblocks/1000genomes/EUR_phase3_autosomes_hg19.${ext}")] }
//             .collect()
// LD_reference2 = Channel
//         .empty()
//         .mix(LD_reference.map{it-> ["SCZ", it[1],it[3],it[5]]})
//         .mix(LD_reference.map{it-> ["HCM", it[1],it[3],it[5]]})

//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect().map{it-> ["SCZ", it]}



workflow UKBB_OR_develop {
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
        UKBB_covariates 

        //out tuple val(meta), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid"), path ("*_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid"), emit: associations
        )
        // PLINK2_ASSOC_GLM.out.associations // ORs [SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/ed/61a5a66ddb71ab22a2c860a368b432/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/ed/61a5a66ddb71ab22a2c860a368b432/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid]
        //     .join(GENERATESNPLISTS.out.processed_ENH_SNP_lists_hg19.map{it->[it[0],it[4]]}, by: [0])//join ENH hg19 csv file
        //     .join(full_GWAS_hg19, by: [0]) //join full GWAS by condition
        //     .view() //[SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/ed/61a5a66ddb71ab22a2c860a368b432/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/ed/61a5a66ddb71ab22a2c860a368b432/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_standard.PHENO1.glm.logistic.hybrid, /rds/general/ephemeral/user/eosimo/ephemeral/UKBB_OR_develop/work/4e/fc02bef580148888c104f3ea7c06ed/ENH_SCZ_hg19.csv, /rds/general/user/eosimo/home/largedirs/scz_GWAS/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz]
    
    R_ANNOTATE_ORs(
        // annotate ORs from previous step with GWAS results and other info,
        //produce OR plots
        PLINK2_ASSOC_GLM.out.associations // ORs
            .join(GENERATESNPLISTS.out.processed_ENH_SNP_lists_hg19.map{it->[it[0],it[4]]}, by: [0])//join ENH hg19 csv file
            .join(full_GWAS_hg19, by: [0]), //join full GWAS by condition
     
        hg38ToHg19_chain,
        GW_LD_blocks
        
        // out: tuple val(meta), path("*_annotated_ORs.csv"),       emit: annotated_ORs
    )
    // R_ANNOTATE_ORs.out.annotated_ORs.view()
    
    // PLINK2_QC_PRUNE_HET (
    //     PLINK_MERGE.out.all_chromosomes_extracted
        
    //     //out tuple val(meta), path ("*.prune.in"), path ("*.het"), emit: pruned_variants_het
    // )
    // R_PRS_QC (
    //     PLINK2_QC_PRUNE_HET.out.pruned_variants_het         //het file
    //         .join(PLINK_MERGE.out.all_chromosomes_extracted, by: [0]) //join merged genotype files
    //         .join(full_GWAS_hg19, by: [0])                       //join full GWAS by condition
        
    //     //out tuple val(meta), path ("*_het_valid_out.sample"), path("*_a1_bim"), path("*_mismatching_SNPs"),  emit: QC_het_a1_mismatch
    // )

    // PLINK_PRODUCE_QC_DATASET (
    //     PLINK_MERGE.out.all_chromosomes_extracted
    //         .join(R_PRS_QC.out.QC_het_a1_mismatch, by: [0]) //join het, A1 bim, mismatching SNPs from previous step
    //         // .view()
    //     //out tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    // )

    // R_PREPARE_MODIF_PRS (
    //     full_GWAS_hg19                                      // full GWAS by condition
    //         .join(R_ANNOTATE_ORs.out.annotated_ORs, by: [0]),// *_annotated_ORs.csv
    //     // GW_LD_blocks
    //     // out tuple val(meta),  path("*_original_dedup_GWAS.tsv"), path("*_substituted_GWAS.tsv"), path("*_tissue_EPeQTL_associations.tsv"), path("*_tissue_facet_associations.tsv"), path("*_all_TS_EPs_associations.tsv"), path("*_merged_GWAS.tsv"), path("*_all_TS_EPs_ZEROP_associations.tsv"), emit: orig_and_modified_GWASes
    // )

    // PLINK_clump (
    //     //tuple val(meta),  path("*_original_dedup_GWAS.tsv"), path("*_substituted_GWAS.tsv"), path("*_tissue_EPeQTL_associations.tsv"), path("*_tissue_facet_associations.tsv"), path("*_all_TS_EPs_associations.tsv"), path("*_merged_GWAS.tsv"), path("*_all_TS_EPs_ZEROP_associations.tsv"), emit: orig_and_modified_GWASes
    //     R_PREPARE_MODIF_PRS.out.orig_and_modified_GWASes    // merged_GWAS for clumping
    //         .join(LD_reference, by: [0])                   //bed bim fam 1000 genomes ref files by DX
        
    // )

    // R_PREPARE_MODIF_PRS_2_LISTS (
    //     PLINK_clump.out.clumped_SNPs 
    //         .join(R_PREPARE_MODIF_PRS.out.orig_and_modified_GWASes, by: [0])
    //     //out tuple val(meta),  path("*_original_NoEPs_associations_clumped.tsv"), path("*_NonOverlapOriginal_TS_EPs_associations_clumped.tsv"), emit: split_GWASes
    // )

    
    // PRSice_calculate_PRS_original(
    //     R_PREPARE_MODIF_PRS.out.orig_and_modified_GWASes          // modified and original deduplicated GWAS
    //         //tuple tuple val(meta),  path("*_original_dedup_GWAS.tsv"), path("*_substituted_GWAS.tsv"), path("*_tissue_EPeQTL_associations.tsv"), path("*_tissue_facet_associations.tsv"), path("*_all_TS_EPs_associations.tsv"), path("*_merged_GWAS.tsv"), emit: orig_and_modified_GWASes
    //         .join(PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC, by: [0])         //QCed UKBB genotypes
    //         //tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    //         .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0])       // also join pheno file
    //         //tuple val(meta), pheno
    //         .join(LD_reference2, by: [0]),
    //         //bed bim fam 1000 genomes ref files by DX
    //     UKBB_covariates,
    //     UKBB_rs34380086_cases
    // )
    // // PRSice_calculate_PRS_substituted(
    // //     R_PREPARE_MODIF_PRS.out.orig_and_modified_GWASes          // modified and original deduplicated GWAS
    // //         //tuple val(meta),  path("*_original_dedup_GWAS.csv"), path("*_substituted_GWAS.csv"), path("*_tissue_EPeQTL_associations.csv"), path("*_tissue_facet_associations.csv"), path("*_all_TS_EPs_associations.csv"), path("*_original_NoEPs_associations.csv"), path("*_NonOverlapOriginal_TS_EPs_associations.csv"),   emit: orig_and_modified_GWASes
    // //         .join(PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC, by: [0])         //QCed UKBB genotypes
    // //         //tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    // //         .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0])       // also join pheno file
    // //         //tuple val(meta), pheno
    // //         .join(LD_reference2, by: [0]),
    // //         //bed bim fam 1000 genomes ref files by DX
    // //     UKBB_covariates,
    // //     UKBB_rs34380086_cases
    // // )
    // PRSice_calculate_PRS_TS_partition(
    //     R_PREPARE_MODIF_PRS.out.orig_and_modified_GWASes          // modified and original deduplicated GWAS
    //         //tuple val(meta),  path("*_original_dedup_GWAS.csv"), path("*_substituted_GWAS.csv"), path("*_tissue_EPeQTL_associations.csv"), path("*_tissue_facet_associations.csv"), path("*_all_TS_EPs_associations.csv"), path("*_original_NoEPs_associations.csv"), path("*_NonOverlapOriginal_TS_EPs_associations.csv"),   emit: orig_and_modified_GWASes
    //         .join(PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC, by: [0])         //QCed UKBB genotypes
    //         //tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    //         .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0])       // also join pheno file
    //         //tuple val(meta), pheno
    //         .join(LD_reference2, by: [0]),
    //         //bed bim fam 1000 genomes ref files by DX
    //     UKBB_covariates,
    //     UKBB_rs34380086_cases
    // )
    // PRSice_calculate_4_partition(
    //     R_PREPARE_MODIF_PRS_2_LISTS.out.split_GWASes
    //         //out tuple val(meta),  path("*_original_NoEPs_associations_clumped.tsv"), path("*_NonOverlapOriginal_TS_EPs_associations_clumped.tsv"), emit: split_GWASes
    //         .join(PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC, by: [0])         //QCed UKBB genotypes
    //         //tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_QC
    //         .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0])       // also join pheno file
    //         //tuple val(meta), pheno
    //         .join(LD_reference2, by: [0]),
    //         //bed bim fam 1000 genomes ref files by DX
    //     UKBB_covariates,
    //     UKBB_rs34380086_cases
    // )
    
    // R_PRS_PPV_plotting(
    //     //tuple val(meta), path("*_original_PRS.summary"), path("*_original_PRS.prsice"), path("*_original_PRS.best")
    //     PRSice_calculate_PRS_original.out.PRS_text_results
    //         //tuple val(meta), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.summary"), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.prsice"), path("*_REC_NonOverlapOriginal_TS_EPs_associations_partition_PLINKclump_PRS.best"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.summary"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.prsice"), path("*_ADD_original_NoEPs_associations_partition_PLINKclump_PRS.best"), emit: PRS_text_split_results
    //         .join(PRSice_calculate_4_partition.out.PRS_text_split_results, by: [0])
    //         .join(PRSice_calculate_PRS_TS_partition.out.PRS_text_results, by: [0])
    // )

    // R_calculate_compound_PRS_fit_and_plot (
    //     PRSice_calculate_4_partition.out.PRS_text_split_results
    //         .join(PRSice_calculate_PRS_original.out.PRS_text_results, by: [0])
    //         .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[3]]}, by: [0]) // join Diagnosis

    // )

    R_plot_GO (
        R_ANNOTATE_ORs.out.annotated_ORs
            .join(full_GWAS_hg19, by: [0])
            .join(enhancer_plus_GWAS_coords.map{it->[it[0],it[1]]}, by: [0]) //initial ENH_P list
    )




}



