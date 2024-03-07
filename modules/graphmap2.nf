// ------------ GRAPHMAP2 -----------
//https://github.com/lbcb-sci/graphmap2#graphmap2---a-highly-sensitive-and-accurate-mapper-for-long-error-prone-reads

/*
* To index
*/ 
process graphmap2_index {
    label 'graphmap2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/graphmap2_indicies", mode: 'copy'

    input:
        path(genome_fasta)

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
    publishDir "${params.outdir}/graphmap2_alignments", pattern: "*graphmap2.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path genome
        path graphmap2_index_files

    output:
        tuple val(sample), path ("*graphmap2.sam"), emit: tuple_sample_sam, optional:true
        tuple val(sample), path ("*graphmap2.mhap"), emit: tuple_sample_mhap, optional:true
        path "*graphmap2.log",  emit: graphmap2_summary

    script:
        fileName = reads[0].baseName
        if ( params.minimap2_options.contains("owler") ){
            """
            graphmap2 ${params.graphmap2_options} -r ${reads} -d ${reads}  -o ${reads.baseName}_graphmap2.mhap 2> ${fileName}_graphmap2.log
            """
        }
        else if ( params.minimap2_options.contains("align") ){
            """
            graphmap2 ${params.graphmap2_options} -i ${graphmap2_index_files} -r ${genome} -d ${reads}  -o ${reads.baseName}_graphmap2.sam 2> ${fileName}_graphmap2.log
            """
        }
        else {
            """
            graphmap2 align ${params.graphmap2_options} -i ${graphmap2_index_files} -r ${genome} -d ${reads}  -o ${reads.baseName}_graphmap2.sam 2> ${fileName}_graphmap2.log
            """
        }

}