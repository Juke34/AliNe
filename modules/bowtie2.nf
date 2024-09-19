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
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*bowtie2.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path genome
        path hisat2_index_files
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam
        path "*bowtie2.log",  emit: bowtie2_summary

    script:
    
    if (params.single_end){
    """
        bowtie2 ${params.bowtie2_options} \\
                -p ${task.cpus} \\
                -x ${genome.baseName} \\
                -S ${reads.baseName.replace('.fastq','')}_bowtie2.sam \\
                -U ${reads} 2> ${reads.baseName}_bowtie2.log
    """
    } else {
    """
        bowtie2 ${params.bowtie2_options} \\
            -p ${task.cpus} \\
            -x ${genome.baseName} \\
            -S ${reads[0].baseName.replace('.fastq','')}_bowtie2.sam \\
            -1 ${reads[0]} -2 ${reads[1]}  2> ${reads[0].baseName.replace('.fastq','')}_bowtie2.log
    """
    }

}