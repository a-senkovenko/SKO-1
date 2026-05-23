process MIXCR {

    tag "$sample"

    publishDir "${params.outdir}/results_mixcr", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample),
          path("${sample}_mixcr/**")

    script:
    """
    mkdir -p ${sample}_mixcr

    mixcr analyze rna-seq \
        --threads ${task.cpus} \
        --species hsa \
        --assemble-longest-contigs \
        ${fastq_1} \
        ${fastq_2} \
        ${sample}_mixcr/${sample}
    """
}