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
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library), val(read_length)
        path hisat2_index_files
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam
        path "${sample}_splicesite.txt"
        path "*hisat2-summary.txt",  emit: hisat2_summary

    script:

        // catch index basename
        index_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
        
        // catch filename
        filename = reads[0].baseName.replace('.fastq','')
       
        // deal with library type - default is unstranded.
        def lib_strand=""
        if (! params.hisat2_options.contains("--rna-strandness ") &&
            ! library.contains("U") &&
            ! params.skip_libray_usage){ // only if -S is not set and if we are not skipping library usage
            if (params.read_type == "short_single" ) { // Lib is not unstranded 
                if ( library.contains("SF") ) {
                    lib_strand = "--rna-strandness F"
                } else if ( library.contains("SR") ) {
                    lib_strand = "--rna-strandness R"
                } 
            } else if ( params.read_type == "short_paired" ) {
                if (library.contains("SR") ) {
                    lib_strand = "--rna-strandness RF"
                } else if (library.contains("SF") ) {
                    lib_strand = "--rna-strandness FR"
                } 
            }
        }
        def read_orientation=""
        if (! params.hisat2_options.contains("--fr ") && 
            ! params.hisat2_options.contains("--rf ") && 
            ! params.hisat2_options.contains("--ff ") && 
            params.read_type == "short_paired" &&
            ! params.skip_libray_usage){ 
            if ( library.contains("I") ){
                read_orientation = "--fr"
            } else if ( library.contains("O") ){
                read_orientation = "--rf"
            }
            else if ( library.contains("M") ){
                read_orientation = "--ff"
            }  
        }

        if (params.read_type == "short_paired"){
            """
            hisat2 ${lib_strand} ${read_orientation} ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
                --new-summary --summary-file ${sample}.hisat2-summary.txt \\
                -p ${task.cpus} -x $index_basename -1 ${reads[0]} -2 ${reads[1]} > ${filename}.sam
            """
            } else {
            """
            hisat2 ${lib_strand} ${read_orientation} ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
                --new-summary --summary-file ${sample}.hisat2-summary.txt \\
                -p ${task.cpus} -x $index_basename -U $reads > ${filename}.sam
            """
        }

}