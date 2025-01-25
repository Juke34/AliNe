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
        salmon index -t ${genome_fasta} -i salmon_index --threads ${task.cpus}
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