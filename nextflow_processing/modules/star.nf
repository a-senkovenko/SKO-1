process STAR {

    tag "$sample"

    publishDir "${params.outdir}/results_star", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample),
          path("${sample}_star/**")

    script:
    """
    mkdir -p ${sample}_star

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${params.star_index} \
        --readFilesIn ${fastq_1} ${fastq_2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample}_star/${sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --quantMode GeneCounts \
        --twopassMode None \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --limitBAMsortRAM 8000000000 \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimOutType WithinBAM SoftClip \
        --chimScoreMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.05 \
        --outFilterMultimapNmax 20 \
        --seedSearchStartLmax 30

    mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.sorted.bam

    rm -rf _STARtmp
    """
}
