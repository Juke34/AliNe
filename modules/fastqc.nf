process fastqc {
    label 'fastqc'
    tag "$sample_id"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val outpath
    val suffix

    output:
    path ("*logs")

    script:

    // Suffix to separate different runs
    def add_suffix = suffix ? "_${suffix}_" : '_'

    """
    mkdir fastqc_${sample_id}${add_suffix}logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}${add_suffix}logs -q ${reads}
    """

}
