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
        def output = "-o ${fastq[0].baseName.replaceAll(/\.(fastq|fq)$/, '')}_trim.fastq.gz" 
        if ( meta.paired ){
            input = "-i ${fastq[0]} -I ${fastq[1]}"
            output = "-o ${fastq[0].baseName.replaceAll(/\.(fastq|fq)$/, '')}_trim.fastq.gz -O ${fastq[1].baseName.replace('.fastq','')}_trim.fastq.gz"
        }

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${meta.id}_fastp_report.html
        """
}