process bowtie_index {
    label 'bowtie'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*.ebwt')

    script:

        """
        bowtie-build --threads ${task.cpus} $genome_fasta ${genome_fasta.baseName}
        """
}

process bowtie {
    label 'bowtie'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*bowtie.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library), val(read_length)
        path genome
        path index_files
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam
        path "*bowtie.log",  emit: bowtie_summary

    script:
    
        def read_orientation=""
        if (! params.bowtie_options.contains("--fr ") && 
            ! params.bowtie_options.contains("--rf ") && 
            ! params.bowtie_options.contains("--ff ") &&
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

        if (params.read_type == "short_paired"){
        """
            bowtie ${params.bowtie_options} ${read_orientation}\\
                -p ${task.cpus} \\
                -x ${genome.baseName} \\
                -S ${reads[0].baseName.replace('.fastq','')}_bowtie.sam \\
                -1 ${reads[0]} -2 ${reads[1]}  2> ${reads[0].baseName.replace('.fastq','')}_bowtie.log
        """
        } else {
        """
            bowtie ${params.bowtie_options} ${read_orientation}\\
                    -p ${task.cpus} \\
                    -x ${genome.baseName} \\
                    -S ${reads} > ${reads.baseName.replace('.fastq','')}_bowtie.sam 2> ${reads.baseName}_bowtie.log
        """
        }

}