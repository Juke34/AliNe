/* ------------ TOPHAT2 -----------

*/

// It use bowtie in the background
process tophat2_index {
    label 'tophat2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.bt2')

    script:

        """
        bowtie2-build ${genome_fasta} ${genome_fasta.baseName}
        """
}

process tophat2 {
    label 'tophat2'
    tag "$sample"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${outpath}/${filename}", path: "*",  mode: 'copy', 
            saveAs: { file -> 
                        if (!file.endsWith('.bam')) {
                            return file  // Publish the file
                        } else {
                            return null  // Exclude the file from being published
                        }
                    }
                

    input:
        tuple val(sample), path(reads), val(library)
        path genome
        path index_files
        path annotation // optional -if provided the option is set in tophat2_options in main.nf
        val outpath

    output:
        tuple val(sample), path ("*accepted_hits.bam"), emit: tuple_sample_bam, optional:true
        path "*align_summary.txt",  emit: tophat2_summary
        path ("*")

    script:

        // set input according to read_type parameter
        def input =  "${reads[0]}"
        if (params.read_type == "short_paired"){
            input =  "${reads[0]} ${reads[1]}" // if short reads check paired or not
        }

        // catch filename
        filename = reads[0].baseName.replace('.fastq','')

        // deal with library type - default is unstranded.
        def lib_strand=""
        if (! params.tophat2_options.contains("--library-type ") &&
            ! library.contains("U") &&
            ! params.skip_libray_usage){ // only if -S is not set and if we are not skipping library usage
            if ( library.contains("SR") ) {
                lib_strand = "--library-type=fr-firststrand"
            } else if ( library.contains("SF") ) {
                lib_strand = "--library-type=fr-secondstrand"
            }
        }


            """
            tophat -p ${task.cpus} \\
                   -o ${filename} \\
                   ${lib_strand} ${params.tophat2_options} ${genome.baseName} \\
                   ${input} 

            # mv the output files
            for i in ${filename}/*; do 
                file=\$(basename \${i});
                mv \${i} ${filename}_tophat2_\${file}; 
            done
            """


}