process PLINK_clump {
    // debug true
    tag "$meta"
    label 'process_low'
    container 'emosyne/plink2:1.2'
    cache 'lenient'

    input:
    // tuple val(meta),  path("*_original_dedup_GWAS.tsv"), path("*_substituted_GWAS.tsv"), path("*_tissue_EPeQTL_associations.tsv"), path("*_tissue_facet_associations.tsv"), path("*_all_TS_EPs_associations.tsv"), path("*_merged_GWAS.tsv"), path("*_all_TS_EPs_ZEROP_associations.tsv"), emit: orig_and_modified_GWASes
    tuple val(meta), path(original_dedup_GWAS), path(substituted_GWAS), path(tissue_EPeQTL_associations), path(tissue_facet_associations), path(all_TS_EPs_associations), path(merged_GWAS), path(all_TS_EPs_ZEROP_associations), path(EUR_phase3_autosomes_hg19_bed), path(EUR_phase3_autosomes_hg19_bim), path(EUR_phase3_autosomes_hg19_fam)
    // each path (PLINKethinicityRelatedness)
    

    output:
    tuple val(meta), path("*.clumped"), emit: clumped_SNPs
    path("*")
    // tuple val(meta), path ("*.bed"), path ("*.bim"), path ("*.fam"), path ("*.log") , emit: all_chromosomes_extracted
    // path "chr_file_list.txt"
    path "versions1.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args 
    def prefix = task.ext.prefix 
    // def mem_mb = task.memory.toMega()

    """
    plink  \\
        --clump ${all_TS_EPs_ZEROP_associations},${original_dedup_GWAS} \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${EUR_phase3_autosomes_hg19_bed.baseName} \\
       --out clumped_merged_SNPs_${meta}  \\
       --threads $task.cpus


    cat <<-END_VERSIONS > versions1.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    
    """
}
