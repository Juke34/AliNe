/* ------------ KALLISTO -----------

*/


process kallisto_index {
    label 'kallisto'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.kallisto_index')

    script:

        """
        kallisto index ${params.kallisto_index_options} -i ${genome_fasta.baseName}.kallisto_index $genome_fasta
        """
}

// kallisto output sorted bam
process kallisto {
    label 'kallisto'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy', pattern: "*.log"

    input:
        tuple val(meta), path(reads)
        path kallisto_index
        val outpath

    output:
        tuple val(meta), path ("${fileName}/*.bam"), emit: tuple_sample_bam
        path "*.log",  emit: kallisto_summary

    script:
        // options for kallisto
        def kallisto_options = meta.kallisto_options ?: ""

        // catch output file prefix 
        fileName = meta.uid + meta.suffix + "_kallisto"

        // For paired-end reads, Kallisto automatically estimates the fragment length distribution from the data and does not require you to specify it manually
        if (meta.paired){
            """
            kallisto quant ${kallisto_options} \
                -t ${task.cpus} \
                --pseudobam \
                -i ${kallisto_index} \
                ${reads[0]} ${reads[1]} -o ${fileName} 2> ${fileName}.log
            
            mv ${fileName}/pseudoalignments.bam ${fileName}/${fileName}.bam
            # in order that the log file contains the name of the output fastq files (MultiQC)
            sed -i 's/${reads[0]}/${fileName}.fastq.gz/' ${fileName}.log
            """
        } else {

            """
            kallisto quant ${kallisto_options} \
                -t ${task.cpus} \
                --pseudobam \
                -i ${kallisto_index} \
                --single \
                ${reads} -o ${fileName} 2> ${fileName}.log

            mv ${fileName}/pseudoalignments.bam ${fileName}/${fileName}.bam
            # in order that the log file contains the name of the output fastq files (MultiQC)
            sed -i 's/${reads[0]}/${fileName}.fastq.gz/' ${fileName}.log
            """
        }
}