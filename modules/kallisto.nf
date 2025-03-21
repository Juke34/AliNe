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
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "${filename}/*.bam", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library), val(read_length)
        path kallisto_index
        val outpath

    output:
        tuple val(sample), path ("${filename}/*.bam"), emit: tuple_sample_bam
        path "*.log",  emit: kallisto_summary

    script:
       
        // catch filename
        filename = reads[0].baseName.replace('.fastq','') + "_kallisto_sorted"
       
        // deal with library type 
        def read_orientation=""
        if (! params.kallisto_options.contains("--fr-stranded ") && 
            ! params.kallisto_options.contains("--rf-stranded ") && 
            params.read_type == "short_paired" && 
            ! params.skip_libray_usage){ 
            if (library.contains("I") ){
                read_orientation = "--fr-stranded"
            } else if (library.contains("O") ){
                read_orientation = "--rf-stranded"
            } 
        }

        // For paired-end reads, Kallisto automatically estimates the fragment length distribution from the data and does not require you to specify it manually
        if (params.read_type == "short_paired"){
            """
            kallisto quant  ${read_orientation} ${params.kallisto_options} \
                -t ${task.cpus} \
                --pseudobam \
                -i ${kallisto_index} \
                ${reads[0]} ${reads[1]} -o ${filename} 2> ${filename}.log
            
            mv ${filename}/pseudoalignments.bam ${filename}/${filename}.bam
            # in order that the log file contains the name of the output fastq files (MultiQC)
            sed -i 's/${reads[0]}/${filename}.fastq.gz/' ${filename}.log
            """
        } else {
            
            // Use read length (-l) and sd (-s) from params?
            def l_s_params = ""
            def read_length_copy = read_length // to avoid error "Variable read_length already defined in the process scope "
            if ( !params.kallisto_options.contains("-l ") ){
                l_s_params += " -l ${read_length}"
            }
            if ( !params.kallisto_options.contains("-s ") ){
                // 10% of read length will be used as Estimated standard deviation of fragment length
                def tenPercent = (read_length_copy.toInteger() * 10 / 100) as int 
                l_s_params += " -s ${tenPercent}"
            }

            """
            kallisto quant  ${read_orientation} ${params.kallisto_options} \
                ${l_s_params} \
                -t ${task.cpus} \
                --pseudobam \
                -i ${kallisto_index} \
                --single \
                ${reads} -o ${filename} 2> ${filename}.log

            mv ${filename}/pseudoalignments.bam ${filename}/${filename}.bam
            # in order that the log file contains the name of the output fastq files (MultiQC)
            sed -i 's/${reads[0]}/${filename}.fastq.gz/' ${filename}.log
            """
        }
}