/*
Here are described all processes related to fastp
fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
See https://github.com/OpenGene/fastp
*/
 
 
process fastp {
    label 'fastp'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'
   
    input:
        tuple val(meta), path(fastq)
        val outpath

    output:
        tuple val(meta), path("*_trim.fastq.gz"), emit: trimmed
        path("${meta.id}_fastp_report.html"), emit: report
   
    script:

        // set input/output according to short_paired parameter
        def input = "-i ${fastq[0]}" 
        def fastqBase0 = AlineUtils.getCleanName(fastq[0])
        def output = "-o ${fastqBase0}_trim.fastq.gz" 
        if ( meta.paired ){
            def fastqBase1 = AlineUtils.getCleanName(fastq[1])
            input = "-i ${fastq[0]} -I ${fastq[1]}"
            output = "-o ${fastqBase0}_trim.fastq.gz -O ${fastqBase1}_trim.fastq.gz"
        }

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${meta.id}_fastp_report.html
        """
}