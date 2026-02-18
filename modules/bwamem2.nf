/* Module related to bwamem2
https://github.com/bwa-mem2/bwa-mem2?tab=readme-ov-file

info:
We are happy to announce that the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. 
For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB. 
There is a substantial reduction in index IO time due to the reduction and hardly any performance impact on read mapping. 
Due to this change in index structure (in commit #4b59796, 10th October 2020), you will need to rebuild the index.
Added MC flag in the output sam file in commit a591e22. Output should match original bwa-mem version 0.7.17.
*/ 

/*
* To index with BWA MEM2
*/
process bwamem2_index {
    label 'bwamem2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:
        """
        bwa-mem2 index $genome_fasta -p ${genome_fasta.baseName}
        """
}

/*
* To align with BWA MEM2
*/
process bwamem2 {
    label 'bwamem2'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwamem2_summary

    script:
        // options for bwa-mem2
        def bwamem2_options = meta.bwamem2_options ?: ""

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_bwamem2"
      
        if (meta.paired){
            """
            bwa-mem2 mem ${bwamem2_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            bwa-mem2 mem ${bwamem2_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}