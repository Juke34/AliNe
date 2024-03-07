process hisat2_index {
    label 'hisat2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/Hisat2_indicies", mode: 'copy'

    input:
    path(genome_fasta)

    output:
    path('*.ht2')

    script:
    """
    hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
    """
}

process hisat2 {
    label 'hisat2'
    tag "$sample"
    publishDir "${params.outdir}/Hisat2_alignments", pattern: "*splicesite.txt", mode: 'copy'
    publishDir "${params.outdir}/Hisat2_alignments", pattern: "*hisat2-summary.txt", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path hisat2_index_files

    output:
        tuple val(sample), path ("${sample}.sam"), emit: tuple_sample_sam
        path "${sample}_splicesite.txt"
        path "*hisat2-summary.txt",  emit: hisat2_summary

    script:
        index_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
        if (params.single_end){
        """
        hisat2 ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
            --new-summary --summary-file ${sample}.hisat2-summary.txt \\
            -p ${task.cpus} -x $index_basename -U $reads > ${sample}.sam
        """
        } else {
        """
        hisat2 ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
            --new-summary --summary-file ${sample}.hisat2-summary.txt \\
            -p ${task.cpus} -x $index_basename -1 ${reads[0]} -2 ${reads[1]} > ${sample}.sam
        """
    }

}