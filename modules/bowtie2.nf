process bowtie2_index {
    label 'bowtie2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.bt2')

    script:

        """
        bowtie2-build --threads ${task.cpus} $genome_fasta ${genome_fasta.baseName}
        """
}

process bowtie2 {
    label 'bowtie2'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bowtie2_summary

    script:
        // options for bowtie2
        def bowtie2_options = meta.bowtie2_options ?: ""
        
        // catch output file prefix 
        def filename = meta.uid + meta.suffix + "_bowtie2"

        if (meta.paired){
        """
            bowtie2 ${bowtie2_options} \\
                -p ${task.cpus} \\
                -x ${genome.baseName} \\
                -S ${filename}.sam \\
                -1 ${reads[0]} -2 ${reads[1]}  2> ${filename}_sorted.log
        """
        } else {
        """
            bowtie2 ${bowtie2_options} \\
                    -p ${task.cpus} \\
                    -x ${genome.baseName} \\
                    -S ${filename}.sam \\
                    -U ${reads} 2> ${filename}_sorted.log
        """
        }

}