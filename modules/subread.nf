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
* To align with graphmap2
*/
process subread {
    label 'subread'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*subread.vcf", mode: 'copy'

    input:
        tuple val(sample), path(fastq)
        path genome
        path index
        val outpath

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam, optional:true
        path "*subread.vcf", emit: subread_vcf, optional:true

    script:

        // set input according to single_end parameter
        def input = params.single_end ? "-r ${fastq}" : "-r ${fastq[0]} -R ${fastq[1]}" // if short reads check paired or not

        // remove fastq.gz
        def fileName = fastq[0].baseName.replace('.fastq','')
        
        // prepare index name
        def index_prefix = genome.baseName + "_index"

        """
        subread-align -T ${task.cpus} -i ${index_prefix} ${input} -o ${fileName}.bam --sortReadsByCoordinates ${params.subread_options}
        """
}