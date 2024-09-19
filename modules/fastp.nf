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

        // set input according to single_end parameter
        def input = params.single_end ? "-i ${fastq}" : "-i ${fastq[0]} -I ${fastq[1]}" // if short reads check paired or not
        // set output according to single_end parameter
        def output = params.single_end ? "-o ${fastq[0].baseName.replace('.fastq','')}_clean.fastq.gz" : "-o ${fastq[0].baseName.replace('.fastq','')}_clean.fastq.gz -O ${fastq[1].baseName.replace('.fastq','')}_clean.fastq.gz" // if short reads check paired or not

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${id}_fastp_report.html
        """
}