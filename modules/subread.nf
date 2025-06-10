/* ------------ subread -----------
https://subread.sourceforge.net/SubreadUsersGuide.pdf
*/

/*
* To index
*/ 
process subread_index {
    label 'subread'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:

        """
        subread-buildindex -o ${genome_fasta.baseName}_index ${genome_fasta}
        """
}

/*
* To align with subread
*/
process subread {
    label 'subread'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*subread.vcf", mode: 'copy'

    input:
        tuple val(meta), path(fastq)
        path genome
        path index
        path annotation // needed in case set in the subread_options
        val outpath

    output:
        tuple val(meta), path ("*.bam"), emit: tuple_sample_bam, optional:true
        path "*subread.vcf", emit: subread_vcf, optional:true
        path "*.log", emit: sublong_log

    script:
        // options for subread
        def subread_options = meta.subread_options ?: ""

        // set input according to short_paired parameter
        def input = "-r ${fastq[0]}"
        if (meta.paired) {
            input =  "-r ${fastq[0]} -R ${fastq[1]}"
        }

        // remove fastq.gz
        def fileName = fastq[0].baseName.replace('.fastq','') + "_subread_sorted"
        
        // prepare index name
        def index_prefix = genome.baseName + "_index"

        // deal with library type
        def read_orientation=""
        if (! subread_options.contains("-S ") &&
            meta.paired && meta.strandedness) { // only if -S is not set and if we are not skipping library usage
            if (meta.strandedness.contains("M") ){
                read_orientation = "-S ff"
            } else if (meta.strandedness.contains("O") ) {
                read_orientation = "-S rf"
            } else if (meta.strandedness.contains("I") ) {
                read_orientation = "-S fr"
            } 
        }
        """
        subread-align -T ${task.cpus} ${read_orientation} -i ${index_prefix} ${input} -o ${fileName}.bam --sortReadsByCoordinates ${subread_options} > ${fileName}_subread_sorted.log 
        """
}

/*
* To index
*/ 
process sublong_index {
    label 'subread'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:

        """
        subread-buildindex -o ${genome_fasta.baseName}_index ${genome_fasta}
        """
}

/*
* To align with sublong
 add  -X to turn on the RNA-seq mode.
*/
process sublong {
    label 'subread'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path index
        val outpath

    output:
        tuple val(meta), path ("*.bam"), emit: tuple_sample_bam, optional:true
        path "*.log", emit: sublong_log

    script:
        // options for sublong
        def sublong_options = meta.sublong_options ?: ""

        // remove fastq.gz
        def fileName = reads[0].baseName.replace('.fastq','') + "_sublong"
        
        // prepare index name
        def index_prefix = genome.baseName + "_index"

        // For paired-end we concat output 
        if (meta.paired){
            """
            sublong -T ${task.cpus} -i ${index_prefix} -r ${reads[0]} -o ${fileName}.bam ${sublong_options} > ${fileName}_sublong.log
            sublong -T ${task.cpus} -i ${index_prefix} -r ${reads[1]} -o ${reads[1].baseName}.bam ${sublong_options} > ${fileName}_sublong.log
            """
        } else {
            """
            sublong -T ${task.cpus} -i ${index_prefix} -r ${reads[0]} -o ${fileName}.bam ${sublong_options} > ${fileName}_sublong.log
            """
        }
}