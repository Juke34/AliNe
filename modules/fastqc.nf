process fastqc {
    label 'fastqc'
    tag "$sample_id"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val outpath

    output:
    path ("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -q ${reads}
    """

}
