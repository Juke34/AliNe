process bwamem3_index {
    label 'bwamem3'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:
        """
        bwa-mem3 index -t ${task.cpus} -p ${genome_fasta.baseName} $genome_fasta
        """
}

process bwamem3 {
    label 'bwamem3'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path("*.sam"), emit: tuple_sample_sam
        path "*.log", emit: bwamem3_summary

    script:
        def bwamem3_options = meta.bwamem3_options ?: ""
        def fileName = meta.uid + meta.suffix + "_bwamem3"

        if (meta.paired){
            """
            bwa-mem3 mem ${bwamem3_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log
            """
        } else {
            """
            bwa-mem3 mem ${bwamem3_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log
            """
        }
}
