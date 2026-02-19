/* Module related to dragmap
https://github.com/Illumina/DRAGMAP

info:
DRAGEN-GATK is a software-only implementation of Illumina's DRAGEN mapper 
that is freely available and open source. It provides the same accuracy and 
functionality as the FPGA-based DRAGEN Bio-IT Platform, but runs on general 
purpose CPUs.
*/ 

/*
* To index with DRAGMAP
*/
process dragmap_index {
    label 'dragmap'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("dragmap_index")

    script:
        """
        mkdir -p dragmap_index
        dragen-os --build-hash-table true --ht-reference ${genome_fasta} --output-directory dragmap_index
        """
}

/*
* To align with DRAGMAP
*/
process dragmap {
    label 'dragmap'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path dragmap_index
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: dragmap_summary

    script:
        // options for dragmap
        def dragmap_options = meta.dragmap_options ?: ""

        // catch output file prefix 
        def fileName = meta.uid + meta.suffix + "_dragmap"
      
        if (meta.paired){
            """
            dragen-os ${dragmap_options} --num-threads ${task.cpus} -r dragmap_index -1 ${reads[0]} -2 ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            dragen-os ${dragmap_options} --num-threads ${task.cpus} -r dragmap_index -1 ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}
