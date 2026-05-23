process KALLISTO {

    tag "$sample"

    publishDir "${params.outdir}/results_kallisto/${sample}_kallisto", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    path(params.transcriptome_fasta)

    output:
    path("${sample}_kallisto"), emit: quant


    script:
    """
    # create index
    if [ ! -f "${params.kallisto_index}" ]; then
        echo "[kallisto] Building index from ${params.transcriptome_fasta}"
        kallisto index -i ${params.kallisto_index} ${params.transcriptome_fasta}
    fi

    echo "[kallisto] Quantifying ${sample}"
    kallisto quant \\
        -i ${params.kallisto_index} \\
        -o ${sample}_kallisto \\
        -t ${task.cpus} \\
        ${fastq_1} ${fastq_2}
    """
}