process fastqc {
    label 'fastqc'
    tag "$sample_id"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """

}
