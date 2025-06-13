process bowtie_index {
    label 'bowtie'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.ebwt')

    script:

        """
        bowtie-build --threads ${task.cpus} $genome_fasta ${genome_fasta.baseName}
        """
}

process bowtie {
    label 'bowtie'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bowtie.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*bowtie.log",  emit: bowtie_summary

    script:

        // options for bowtie
        def bowtie_options = meta.bowtie_options ?: ""

        // catch filename
        def outBaseName = AlineUtils.getCleanName(reads) + "_bowtie"

        if (meta.paired){
        """
            bowtie ${bowtie_options} \\
                -p ${task.cpus} \\
                -x ${genome.baseName} \\
                -S ${outBaseName}.sam \\
                -1 ${reads[0]} -2 ${reads[1]}  2> ${outBaseName}.log
        """
        } else {
        """
            bowtie ${bowtie_options} \\
                    -p ${task.cpus} \\
                    -x ${genome.baseName} \\
                    -S ${reads} > ${outBaseName}.sam 2> ${outBaseName}.log
        """
        }

}