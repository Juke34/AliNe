process fastqc {
    label 'fastqc'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        val outpath
        val suffix

    output:
        path ("*logs")

    script:

        // Suffix to separate different runs
        def sample_id = meta.id
        def add_suffix = suffix ? "_${suffix}_" : '_'

        """
        mkdir fastqc_${sample_id}${add_suffix}logs
        fastqc -t ${task.cpus} -o fastqc_${sample_id}${add_suffix}logs -q ${reads}
        """
}

// To take in consideration the index coming along when aligned files are provided
process fastqc_ali {
    label 'fastqc'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(reads), path(index)
        val outpath
        val suffix

    output:
        path ("*logs")

    script:

        // Suffix to separate different runs
        def sample_id = meta.id
        def add_suffix = suffix ? "_${suffix}_" : '_'

        """
        mkdir fastqc_${sample_id}${add_suffix}logs
        fastqc -t ${task.cpus} -o fastqc_${sample_id}${add_suffix}logs -q ${reads}
        """
}