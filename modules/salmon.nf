/*
Here are described all processes related to salmon
https://github.com/COMBINE-lab/salmon
*/
 
process salmon_index {
    label 'salmon'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy', enabled: params.aligner.contains('salmon')

    input:
        path genome_fasta
        val outpath

    output:
        path "salmon_index", emit: index

    script:
        """
        salmon index ${params.salmon_index_options} -t ${genome_fasta} -i salmon_index --threads ${task.cpus}
        """
}

process salmon_guess_lib {
    label 'salmon'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.json", mode: 'copy'
   
    input:
        tuple val(meta), path(fastq)
        path salmon_index
        val outpath

    output:
        tuple val(meta), env(LIBTYPE), emit: tuple_id_libtype
        path "*lib_format_counts.json"
   
    script:

        // set input according to read_type parameter
        def input =  "-r ${fastq[0]}"
        if (meta.paired){
            input =  "-1 ${fastq[0]} -2 ${fastq[1]}" // if short reads check paired or not
        }
        // remove the extension from fastq file name
        def output = AlineUtils.getCleanName(fastq)

        """
            salmon quant -i ${salmon_index} -l A ${input} --thread ${task.cpus} -o ${output} --minAssignedFrags 2 
            # extract the result
            LIBTYPE=\$(grep expected_format ${output}/lib_format_counts.json | awk '{print \$2}' | tr -d '",\n')
            # change output name
            mv ${output}/lib_format_counts.json ${meta.file_id}_AliNe_lib_format_counts.json
        """

}

//  Use salmon as aligner - output sorted sam
process salmon {
    label 'salmon'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'
   
    input:
        tuple val(meta), path(fastq)
        path salmon_index
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: salmon_summary
   
    script:
        // options for Salmon
        def salmon_options = meta.salmon_options ?: ""

        // set input according to read_type parameter
        def input =  "-r ${fastq[0]}"
        if ( meta.paired ){
            input =  "-1 ${fastq[0]} -2 ${fastq[1]}" // if short reads check paired or not
        }

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_salmon"
       
        // Salmon automatically estimates the fragment length distribution for paired-end reads (like Kallisto)
        if ( meta.paired ){
            """
                salmon quant -i ${salmon_index} ${salmon_options} \
                    ${input} \
                    --thread ${task.cpus} \
                    --writeMappings \
                    --output ${fileName} > ${fileName}.sam 2> ${fileName}.log
            """
        } else {
            """
                salmon quant -i ${salmon_index} ${salmon_options} \
                    ${input} \
                    --thread ${task.cpus} \
                    --writeMappings \
                    --output ${fileName} > ${fileName}.sam 2> ${fileName}.log
            """
        }
}