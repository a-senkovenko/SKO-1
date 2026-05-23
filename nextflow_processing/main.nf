nextflow.enable.dsl=2

params.outdir = "results"

params.run_download = true
params.run_kallisto = false
params.run_fastp = false
params.run_star   = true
params.run_mixcr   = false
params.run_trust4  = false
params.run_arcashla = false

params.trust4_ref  = "/home/bioinf2026/sko1/TRUST4/human_IMGT+C.fa"
params.trust4_fa   = "/home/bioinf2026/sko1/TRUST4/hg38_bcrtcr.fa"

params.transcriptome_fasta = "/home/bioinf2026/sko1/reference/gencode.v38.transcripts.fa.gz"
params.kallisto_index     = "/home/bioinf2026/sko1/reference/human_transcripts.idx"

params.samples_csv = "sample_data/sample.csv"
params.fastq_dir   = "fastqs"
params.bam_dir = "results/results_star/"
params.star_index  = "/home/bioinf2026/sko1/reference_star"


include { DOWNLOAD_ARIA }       from './modules/download_aria.nf'
include { KALLISTO }            from './modules/kallisto.nf'
include { STAR }                from './modules/star.nf'
include { MIXCR }               from './modules/MIXCR.nf'
include { TRUST4 }              from './modules/TRUST4.nf'
include { FASTP }               from './modules/fastp.nf'
include { ARCAS_HLA }           from './modules/arcasHLA.nf'

workflow {

    /*
    STEP 1
    DOWNLOAD FASTQ (optional)
    */

    if(params.run_download) {

        log.info "Mode: DOWNLOAD"

        fastq_pairs_ch = Channel
            .fromPath(params.samples_csv)
            .splitCsv(header:false)
            .map { row -> tuple(row[0].trim(), row[0].trim()) }
            | DOWNLOAD_ARIA

    } else {

        log.info "Mode: LOCAL FASTQ"

        fastq_pairs_ch = Channel.fromFilePairs(
            "${params.fastq_dir}/*_{1,2}.fastq.gz",
            flat: true
        )
    }

    fastq_pairs_ch.view { "INPUT: $it" }

    /*
    STEP 2
    FASTP / STAR / MIXCR
    */

    if(params.run_fastp){
        fastp_results = FASTP(fastq_pairs_ch)
    }

    if(params.run_star){
        star_bams = STAR(fastq_pairs_ch)
    }

    if(params.run_mixcr){
        mixcr_results = MIXCR(fastq_pairs_ch)
    }

    /*
    STEP 3
    KALLISTO
    */

    if(params.run_kallisto){
        KALLISTO(fastq_pairs_ch, params.transcriptome_fasta)
    }

    /*
    STEP 4
    TRUST4
    */

    if(params.run_trust4){
        trust_fastq = TRUST4(fastq_pairs_ch)
    }

    /*
    STEP 5
    ARCAS_HLA
    */

    if(params.run_arcashla){

        if(!params.skip_star){
            ARCAS_HLA(star_bams)
        } else {
            //error "ARCAS_HLA requires STAR BAMs"
            star_bams = Channel.fromPath("${params.bam_dir}/**/*.sorted.bam")
                .map { bam ->
                    tuple(bam.baseName.replace(".sorted",""), bam)
            }
            ARCAS_HLA(star_bams)
        }
    }
}