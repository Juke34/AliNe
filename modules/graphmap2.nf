// ------------ GRAPHMAP2 -----------
//https://github.com/lbcb-sci/graphmap2#graphmap2---a-highly-sensitive-and-accurate-mapper-for-long-error-prone-reads

/*
* To index
*/ 
process graphmap2_index {
    label 'graphmap2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:

        """
        graphmap2 align -t ${task.cpus} -I -r ${genome_fasta}
        """
}

/*
* To align with graphmap2
*/
process graphmap2 {
    label 'graphmap2'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*graphmap2.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(readtype), val(read_length)
        path genome
        path graphmap2_index_files
        path annotation // needed in case set in the params.graphmap2_options
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam, optional:true
        tuple val(sample), path ("*.mhap"), emit: tuple_sample_mhap, optional:true
        path "*graphmap2.log",  emit: graphmap2_summary

    script:
        // catch filename
        fileName = reads[0].baseName.replace('.fastq','')

        // Check if the owler option is set
        if ( params.graphmap2_options.contains("owler") ){

            if (params.read_type == "short_paired"){
                // For paired-end we concat output 
                """
                graphmap2 ${params.graphmap2_options} -t ${task.cpus} -r ${reads[0]} -d ${reads[0]}  -o ${fileName}_graphmap2.mhap 2> ${fileName}_graphmap2.log
                graphmap2 ${params.graphmap2_options} -t ${task.cpus} -r ${reads[1]} -d ${reads[1]}  -o ${reads[1].baseName}_graphmap2.mhap 2> ${reads[1].baseName}_graphmap2.log
                cat ${fileName}_graphmap2.mhap > ${fileName}_graphmap2_concatR1R2.mhap
                rm ${fileName}_graphmap2.mhap
                cat ${reads[1].baseName}_graphmap2.mhap >> ${fileName}_graphmap2_concatR1R2.mhap
                rm ${reads[1].baseName}_graphmap2.mhap
                """
            } else {
                """
                graphmap2 ${params.graphmap2_options} -t ${task.cpus} -r ${reads[0]} -d ${reads[0]}  -o ${fileName}_graphmap2.mhap 2> ${fileName}_graphmap2.log
                """
            }
        }
        // Align case
        else {
            // add align if absent from the command line
            if (! params.graphmap2_options.contains("align") ){
                graphmap2_options = "align ${params.graphmap2_options}"
            }

            // For paired-end we concat output 
            if (params.read_type == "short_paired"){
                
                """
                graphmap2 ${graphmap2_options} -i ${graphmap2_index_files} -t ${task.cpus} -r ${genome} -d ${reads[0]}  -o ${fileName}_graphmap2.sam 2> ${fileName}_graphmap2.log
                graphmap2 ${graphmap2_options} -i ${graphmap2_index_files} -t ${task.cpus} -r ${genome} -d ${reads[1]}  -o ${reads[1].baseName}_graphmap2.sam 2> ${reads[1].baseName}_graphmap2.log
                
                # Merge sam
                cat ${fileName}_graphmap2.sam > ${fileName}_graphmap2_concatR1R2.sam
                rm ${fileName}_graphmap2.sam
                awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}_graphmap2.sam >> ${fileName}_graphmap2_concatR1R2.sam
                rm ${reads[1].baseName}_graphmap2.sam
                """
            } else {
                """
                
                graphmap2 ${graphmap2_options} -i ${graphmap2_index_files} -t ${task.cpus} -r ${genome} -d ${reads[0]}  -o ${fileName}_graphmap2.sam 2> ${fileName}_graphmap2.log
                """
            }
        }
}