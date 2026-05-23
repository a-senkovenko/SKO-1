process TRUST4  {

    publishDir "${params.outdir}/results_trust4/${sample}_trust4", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample),
          path("${sample}_report.tsv"),
          path("${sample}_airr_align.tsv"),
          path("${sample}_airr.tsv"),
          path("${sample}_annot.fa"),
          path("${sample}_cdr3.out")

    script:
    """
    export PATH=\$PATH:/home/bioinf2026/sko1/TRUST4
    run-trust4 \
      -f ${params.trust4_fa} \
      --ref ${params.trust4_ref} \
      -1 ${fastq_1} \
      -2 ${fastq_2} \
      -t ${task.cpus} \
      -o ${sample}
    """
}