/* ------------ HISTA2 -----------

*/


process hisat2_index {
    label 'hisat2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.ht2')

    script:

        """
        hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
        """
}

process hisat2 {
    label 'hisat2'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path hisat2_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*_splicesite.txt"
        path "*_summary.txt",  emit: hisat2_summary

    script:
        // options for hisat2
        def hisat2_options = meta.hisat2_options ?: ""

        // catch index basename
        index_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
        
        // catch filename
        def filename = AlineUtils.getCleanName(reads) + "_hisat2"
       
        // deal with library type - default is unstranded.
        def lib_strand=""
        if (! params.hisat2_options.contains("--rna-strandness ") &&
              meta.strandedness && ! meta.strandedness.contains("U")
             ){ // only if -S is not set and if we are not skipping library usage
            if ( meta.paired ) {
                if (meta.strandedness.contains("SR") ) {
                    lib_strand = "--rna-strandness RF"
                } else if (meta.strandedness.contains("SF") ) {
                    lib_strand = "--rna-strandness FR"
                }
            } // unpair
            else { // Lib is not unstranded 
                if ( meta.strandedness.contains("SF") ) {
                    lib_strand = "--rna-strandness F"
                } else if ( meta.strandedness.contains("SR") ) {
                    lib_strand = "--rna-strandness R"
                } 
            }
        }
        // read_orientation is only for paired-end reads
        def read_orientation=""
        if (! params.hisat2_options.contains("--fr ") && 
            ! params.hisat2_options.contains("--rf ") && 
            ! params.hisat2_options.contains("--ff ") && 
            meta.paired && meta.strandedness ){ 
            if ( meta.strandedness.contains("I") ){
                read_orientation = "--fr"
            } else if ( meta.strandedness.contains("O") ){
                read_orientation = "--rf"
            }
            else if ( meta.strandedness.contains("M") ){
                read_orientation = "--ff"
            }  
        }

        if (meta.paired) {
            """
            hisat2 ${lib_strand} ${read_orientation} ${params.hisat2_options} --novel-splicesite-outfile ${filename}_splicesite.txt \\
                --new-summary --summary-file ${filename}_sorted_summary.txt \\
                -p ${task.cpus} -x $index_basename -1 ${reads[0]} -2 ${reads[1]} > ${filename}.sam
            """
            } else {
            """
            hisat2 ${lib_strand} ${read_orientation} ${params.hisat2_options} --novel-splicesite-outfile ${filename}_splicesite.txt \\
                --new-summary --summary-file ${filename}_sorted_summary.txt \\
                -p ${task.cpus} -x $index_basename -U $reads > ${filename}.sam
            """
        }

}