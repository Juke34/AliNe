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
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bowtie2.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bowtie2_summary

    script:
        // options for bowtie2
        def bowtie2_options = meta.bowtie2_options ?: ""
        
        def read_orientation=""
        if (! params.bowtie2_options.contains("--fr ") && 
            ! params.bowtie2_options.contains("--rf ") && 
            ! params.bowtie2_options.contains("--ff ") &&
              meta.paired && 
            meta.strandedness){ 
            if (meta.strandedness.contains("I") ){
                read_orientation = "--fr"
            } else if (meta.strandedness.contains("O") ){
                read_orientation = "--rf"
            } else if (meta.strandedness.contains("M") ){
                read_orientation = "--ff"
            }  
        }
        
        // catch filename
        def filename = AlineUtils.getCleanName(reads) + "_bowtie2"

        if (meta.paired){
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