process fastqc {
    label 'fastqc'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        val outpath
        val suffix

    output:
        path ("*logs")

    script:

        // Suffix to separate different runs
        def add_suffix = suffix ? "_${suffix}_" : '_'

        // catch output file prefix 
        def fileName = "fastqc_${meta.uid}${add_suffix}logs"

        """
        mkdir fastqc_${fileName}
        fastqc -t ${task.cpus} -o fastqc_${fileName} -q ${reads}
        """
}

// To take in consideration the index coming along when aligned files are provided
process fastqc_ali {
    label 'fastqc'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(reads), path(index)
        val outpath
        val suffix

    output:
        path ("*logs")

    script:

        // Suffix to separate different runs
        def file_id = meta.uid
        def add_suffix = suffix ? "_${suffix}_" : '_'

        """
        mkdir fastqc_${file_id}${add_suffix}logs
        fastqc -t ${task.cpus} -o fastqc_${file_id}${add_suffix}logs -q ${reads}
        """
}