process bowtie2_index {
    label 'bowtie2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.bt2')

    script:

        """
        bowtie2-build --threads ${task.cpus} $genome_fasta ${genome_fasta.baseName}
        """
}

process bowtie2 {
    label 'bowtie2'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*bowtie2.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library), val(read_length)
        path genome
        path index_files
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bowtie2_summary

    script:
    
        def read_orientation=""
        if (! params.bowtie2_options.contains("--fr ") && 
            ! params.bowtie2_options.contains("--rf ") && 
            ! params.bowtie2_options.contains("--ff ") &&
            params.read_type == "short_paired" && 
            ! params.skip_libray_usage){ 
            if (library.contains("I") ){
                read_orientation = "--fr"
            } else if (library.contains("O") ){
                read_orientation = "--rf"
            } else if (library.contains("M") ){
                read_orientation = "--ff"
            }  
        }
        
        // catch filename
        filename = reads[0].baseName.replace('.fastq','') + "_bowtie2"

        if (params.read_type == "short_paired"){
        """
            bowtie2 ${params.bowtie2_options} ${read_orientation}\\
                -p ${task.cpus} \\
                -x ${genome.baseName} \\
                -S ${filename}.sam \\
                -1 ${reads[0]} -2 ${reads[1]}  2> ${filename}_sorted.log
        """
        } else {
        """
            bowtie2 ${params.bowtie2_options} ${read_orientation}\\
                    -p ${task.cpus} \\
                    -x ${genome.baseName} \\
                    -S ${filename}.sam \\
                    -U ${reads} 2> ${filename}_sorted.log
        """
        }

}