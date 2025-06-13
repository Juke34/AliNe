/* ------------ HISTA2 -----------

*/


process hisat2_index {
    label 'hisat2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.ht2')

    script:

        """
        hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
        """
}

process hisat2 {
    label 'hisat2'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path hisat2_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*_splicesite.txt"
        path "*_summary.txt",  emit: hisat2_summary

    script:
        // options for hisat2
        def hisat2_options = meta.hisat2_options ?: ""

        // catch index basename
        index_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
        
        // catch filename
        def filename = AlineUtils.getCleanName(reads) + "_hisat2"
       
        // alignment
        if (meta.paired) {
            """
            hisat2 ${hisat2_options} --novel-splicesite-outfile ${filename}_splicesite.txt \\
                --new-summary --summary-file ${filename}_sorted_summary.txt \\
                -p ${task.cpus} -x $index_basename -1 ${reads[0]} -2 ${reads[1]} > ${filename}.sam
            """
            } else {
            """
            hisat2 ${hisat2_options} --novel-splicesite-outfile ${filename}_splicesite.txt \\
                --new-summary --summary-file ${filename}_sorted_summary.txt \\
                -p ${task.cpus} -x $index_basename -U $reads > ${filename}.sam
            """
        }

}