nextflow.enable.dsl=2

params.samplesheet = "./samplesheet.csv"

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> 
        tuple(
            row.sample, 
            file(row.fastq1), 
            file(row.fastq2), 
            row.library, 
            row.platform, 
            row.platform_unit
        ) 
    }
    .groupTuple(by: 0)
    .set { grouped_fastqs }

// 1. Mapping, sorting, deduplication, indel realignment (alt-aware or not)
process map_with_rg {
    tag "$sample_id"
    publishDir "results/mapped", mode: 'copy'
    input:
    tuple val(sample_id), val(fastq_tuples)
    output:
    tuple val(sample_id), path("${sample_id}_*.bam")
    script:
    """
    i=0
    for tuple in \$(seq 0 \$((\${#fastq_tuples[@]} / 4 - 1))); do
        fq1=\${fastq_tuples[\$((tuple*4+0))]}
        fq2=\${fastq_tuples[\$((tuple*4+1))]}
        lib=\${fastq_tuples[\$((tuple*4+2))]}
        pl=\${fastq_tuples[\$((tuple*4+3))]}
        pu=\${fastq_tuples[\$((tuple*4+4))]}
        rgid="${sample_id}.\$i"
        # Choose alt-aware or not based on genomeType param
        if [[ "\$GENOME_TYPE" == "has_alt_contigs" ]]; then
            mapSortDedupMarkIndels_alt_aware_phoenix.sh -S ${sample_id} -1 \$fq1 -2 \$fq2 -L \$lib -P \$pl -U \$pu -R \$REF_FASTA -o . -I \$rgid
        else
            mapSortDedupMarkIndels_no_alt_phoenix.sh -S ${sample_id} -1 \$fq1 -2 \$fq2 -L \$lib -P \$pl -U \$pu -R \$REF_FASTA -o . -I \$rgid
        fi
        mv ${sample_id}.*.marked.sort.bwa.*.bam ${sample_id}_\${i}.bam
        i=\$((i+1))
    done
    """
}

// 2. Merge BAMs for each sample
process merge_bams {
    tag "$sample_id"
    publishDir "results/merged", mode: 'copy'
    input:
    tuple val(sample_id), path(bams)
    output:
    tuple val(sample_id), path("${sample_id}.merged.bam")
    script:
    """
    sambamba merge -t 8 ${sample_id}.merged.bam ${bams.join(' ')}
    """
}

// 3. BQSR covariate calculation
process bqsr_covariates {
    tag "$sample_id"
    publishDir "results/bqsr", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.recal.table"), path(bam)
    script:
    """
    gatk BaseRecalibrator \
        -R \$REF_FASTA \
        -I $bam \
        --known-sites \$KNOWN_SITES \
        -O ${sample_id}.recal.table
    """
}

// 4. Apply BQSR (arrayed by interval if needed)
process apply_bqsr {
    tag "$sample_id"
    publishDir "results/bqsr", mode: 'copy'
    input:
    tuple val(sample_id), path(recal_table), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam")
    script:
    """
    gatk ApplyBQSR \
        -R \$REF_FASTA \
        -I $bam \
        --bqsr-recal-file $recal_table \
        -O ${sample_id}.bqsr.bam
    """
}

// 5. HaplotypeCaller (arrayed by interval if needed)
process haplotype_caller {
    tag "$sample_id"
    publishDir "results/gvcf", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz")
    script:
    """
    gatk HaplotypeCaller \
        -R \$REF_FASTA \
        -I $bam \
        -O ${sample_id}.g.vcf.gz \
        -ERC GVCF
    """
}

// 6. Gather gVCFs
process gather_gvcfs {
    tag "$sample_id"
    publishDir "results/gvcf", mode: 'copy'
    input:
    tuple val(sample_id), path(gvcf)
    output:
    tuple val(sample_id), path("${sample_id}.gathered.g.vcf.gz")
    script:
    """
    gatk GatherVcfs \
        -R \$REF_FASTA \
        -I $gvcf \
        -O ${sample_id}.gathered.g.vcf.gz
    """
}

// 7. Collect WGS metrics
process collect_wgs_metrics {
    tag "$sample_id"
    publishDir "results/metrics", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.wgs.metrics.txt")
    script:
    """
    picard CollectWgsMetrics \
        I=$bam \
        O=${sample_id}.wgs.metrics.txt \
        R=\$REF_FASTA
    """
}

// 8. Joint genotyping
process joint_genotyping {
    publishDir "results/joint", mode: 'copy'
    input:
    set val(sample_id), path(gvcf) from gather_gvcfs.out.collect()
    output:
    path("joint_genotyped.vcf.gz")
    script:
    """
    gatk GenotypeGVCFs \
        -R \$REF_FASTA \
        $(ls *.gathered.g.vcf.gz | sed 's/^/-V /') \
        -O joint_genotyped.vcf.gz
    """
}

// 9. Merge and VQSR
process merge_vqsr {
    publishDir "results/final", mode: 'copy'
    input:
    path vcf from joint_genotyping.out
    output:
    path("final.vcf.gz")
    script:
    """
    gatk VariantRecalibrator \
        -R \$REF_FASTA \
        -V $vcf \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \$HAPMAP \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 \$OMNI \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 \$G1000 \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$DBSNP \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
        -mode SNP \
        -O recalibrate_SNP.recal \
        --tranches-file recalibrate_SNP.tranches \
        --rscript-file recalibrate_SNP_plots.R
    gatk ApplyVQSR \
        -R \$REF_FASTA \
        -V $vcf \
        --recal-file recalibrate_SNP.recal \
        --tranches-file recalibrate_SNP.tranches \
        -mode SNP \
        -O final.vcf.gz
    """
}

workflow {
    grouped_fastqs | map_with_rg | merge_bams | bqsr_covariates | apply_bqsr
    apply_bqsr.out | haplotype_caller | gather_gvcfs | joint_genotyping | merge_vqsr
    apply_bqsr.out | collect_wgs_metrics
}