process FASTP {

    tag "$sample"

    publishDir "${params.outdir}/fastp_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample),
          path("${sample}_trimmed_1.fastq.gz"),
          path("${sample}_trimmed_2.fastq.gz"),
          emit: reads

    path("${sample}.fastp.html"), emit: html
    path("${sample}.fastp.json"), emit: json

    script:
    """
    fastp \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 30 \
        -i ${reads1} \
        -I ${reads2} \
        -o ${sample}_trimmed_1.fastq.gz \
        -O ${sample}_trimmed_2.fastq.gz \
        -h ${sample}.fastp.html \
        -j ${sample}.fastp.json
    """
}