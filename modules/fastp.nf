/*
Here are described all processes related to fastp
fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
See https://github.com/OpenGene/fastp
*/
 
 
process fastp {
    label 'fastp'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'
   
    input:
        tuple val(meta), path(fastq)
        val outpath

    output:
        tuple val(meta), path("*_fastp.fastq.gz"), emit: trimmed
        path("*_report.html"), emit: report
   
    script:

        // add suffix to meta for output files
        def suffix = meta.suffix ? "${meta.suffix}_fastp" : '_fastp'
        meta = meta + [suffix: suffix]

        // set input/output according to short_paired parameter
        def input = "-i ${fastq[0]}" 
        def output = "-o ${meta.uid}${suffix}.fastq.gz" 
        if ( meta.paired ){
            input = "-i ${fastq[0]} -I ${fastq[1]}"     
            output = "-o ${meta.file_id[0]}${suffix}.fastq.gz -O ${meta.file_id[1]}${suffix}.fastq.gz"
        }

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${meta.uid}${suffix}_report.html
        """
}