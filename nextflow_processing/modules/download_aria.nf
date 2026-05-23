process DOWNLOAD_ARIA {

    tag "$srr"

    publishDir "${params.fastq_dir}", mode: 'copy'

    input:
    tuple val(sample), val(srr)

    output:
    tuple val(sample),
      path("${srr}_1.fastq.gz"),
      path("${srr}_2.fastq.gz")

    script:
    """
    # Run script from bin directory
    ${projectDir}/bin/download_fastq_from_csv.sh ${srr}
    """
}