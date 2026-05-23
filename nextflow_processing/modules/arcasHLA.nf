process ARCAS_HLA {

    tag "$sample"

    scratch false

    publishDir "${params.outdir}/arcas_hla", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample),
          path("${sample}_arcas/**")

    script:
    """
    mkdir -p ${sample}_arcas

    arcasHLA extract \
        ${bam} \
        -o tmp_${sample} \
        -t ${task.cpus} \
        --temp .

    arcasHLA genotype \
        tmp_${sample}/${sample}.sorted.extracted.1.fq.gz \
        tmp_${sample}/${sample}.sorted.extracted.2.fq.gz \
        -o ${sample}_arcas \
        -t ${task.cpus} \
        --temp .

    rm -rf tmp_${sample}
    """
}