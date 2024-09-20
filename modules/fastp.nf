/*
Here are described all processes related to fastp
fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
See https://github.com/OpenGene/fastp
*/
 
 
process fastp {
    label 'fastp'
    publishDir "${params.outdir}/${outpath}", mode: 'copy'
   
    input:
        tuple val(id), path(fastq)
        val outpath

    output:
        tuple val(id), path("*_clean.fastq.gz"), emit: trimmed
        path("${id}_fastp_report.html"), emit: report
   
    script:

        // set input/output according to short_paired parameter
        def input = "-i ${fastq[0]}" 
        def output = "-o ${fastq[0].baseName.replace('.fastq','')}_clean.fastq.gz" 
        if (params.read_type == "short_paired"){
            input = "-i ${fastq[0]} -I ${fastq[1]}"
            output = "-o ${fastq[0].baseName.replace('.fastq','')}_clean.fastq.gz -O ${fastq[1].baseName.replace('.fastq','')}_clean.fastq.gz"
        }

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${id}_fastp_report.html
        """
}