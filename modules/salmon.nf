/*
Here are described all processes related to salmon
https://github.com/COMBINE-lab/salmon
*/
 
process salmon_index {
    label 'salmon'

    input:
        path genome_fasta

    output:
        path "salmon_index", emit: index

    script:
        """
        salmon index ${params.salmon_index_options} -t ${genome_fasta} -i salmon_index --threads ${task.cpus}
        """
}

process salmon_guess_lib {
    label 'salmon'
    publishDir "${params.outdir}/${outpath}", pattern: "*/*.json", mode: 'copy'
   
    input:
        tuple val(id), path(fastq)
        path salmon_index
        val outpath

    output:
        tuple val(id), env(LIBTYPE), emit: tuple_id_libtype
        path "*/*lib_format_counts.json"
   
    script:

        // set input according to read_type parameter
        def input =  "-r ${fastq[0]}"
        if (params.read_type == "short_paired"){
            input =  "-1 ${fastq[0]} -2 ${fastq[1]}" // if short reads check paired or not
        }
        def output = "${fastq[0].baseName.replace('.fastq','')}"

        """
            salmon quant -i ${salmon_index} -l A ${input} --thread ${task.cpus} -o ${output} --minAssignedFrags 2 
            # extract the result
            LIBTYPE=\$(grep expected_format ${output}/lib_format_counts.json | awk '{print \$2}' | tr -d '",\n')
            # change output name
            mv ${output}/lib_format_counts.json ${output}/${id}_lib_format_counts.json
        """

}

// Process not related to the tool but to the library guessing step made with salmon
process set_tuple_withUserLib{
    label 'salmon'
   
    input:
        tuple val(id), path(fastq)

    output:
        tuple val(id), val(params.library_type), emit: tuple_id_libtype

   
    script:

        """
        """
}

//  Use salmon as aligner - output sorted sam
process salmon {
    label 'salmon'
    publishDir "${params.outdir}/${outpath}", pattern: "*/*.json", mode: 'copy'
   
    input:
        tuple val(sample), path(fastq), val(library), val(read_length)
        path salmon_index
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: salmon_summary
   
    script:

        // set input according to read_type parameter
        def input =  "-r ${fastq[0]}"
        if (params.read_type == "short_paired"){
            input =  "-1 ${fastq[0]} -2 ${fastq[1]}" // if short reads check paired or not
        }

        // deal with library type 
        def read_orientation=""
        if (! params.salmon_options.contains("-l ") && ! params.salmon_options.contains("--libType ") &&
            ! params.skip_libray_usage){ 
                read_orientation = "-l ${library}"
        }

        // catch filename
        def filename = "${fastq[0].baseName.replace('.fastq','')}_salmon"
       
        // Salmon automatically estimates the fragment length distribution for paired-end reads (like Kallisto)
        if (params.read_type == "short_paired"){
            """
                salmon quant -i ${salmon_index} ${params.salmon_options} \
                    ${read_orientation} \
                    ${input} \
                    --thread ${task.cpus} \
                    --writeMappings \
                    --output ${filename} > ${filename}.sam 2> ${filename}.log
            """
        } else {
            
            // Use read length (--fldMean) and sd (--fldSD) from params?
            def l_s_params = ""
            def read_length_copy = read_length // to avoid error "Variable read_length already defined in the process scope "
            if ( !params.salmon_options.contains("--fldMean ") ){
                l_s_params += " --fldMean ${read_length}"
            }
            if ( !params.salmon_options.contains("--fldSD ") ){
                // 10% of read length will be used as Estimated standard deviation of fragment length
                def tenPercent = (read_length_copy.toInteger() * 10 / 100) as int 
                l_s_params += " --fldSD ${tenPercent}"
            }

            """
                salmon quant -i ${salmon_index} ${params.salmon_options} \
                    ${l_s_params} \
                    ${read_orientation} \
                    ${input} \
                    --thread ${task.cpus} \
                    --writeMappings \
                    --output ${filename} > ${filename}.sam 2> ${filename}.log
            """
        }
}