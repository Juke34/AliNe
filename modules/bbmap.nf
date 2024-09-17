/* Module related to bbmmap
BBMap is a global aligner. That means it looks for the highest-scoring alignment taking into account all bases in a sequence. A local aligner would look for the best-scoring local alignment, meaning an alignment where the ends are possibly clipped off.
*/ 

process bbmap_index {
    label 'bbmap'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('ref')

    script:
    
        // set memory to 4GB if not specified
        def avail_mem = task.memory ? task.memory.toGiga() : 4    
    
        """
        bbmap.sh ref=${genome_fasta} threads=2 -Xmx${avail_mem}g
        """
}

process bbmap {
    label 'bbmap'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*bbmap.log", mode: 'copy'

    input:
        tuple val(sample), path(fastq)
        path genome_index
        path hisat2_index_files
        val outpath

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam
        path "*bbmap.log",  emit: bbmap_summary

    script:

    // set memory to 4GB if not specified
    def avail_mem = task.memory ? task.memory.toGiga() : 4    
    
    // set tool according to read length
    def tool = params.long_reads ? "mapPacBio.sh" : "bbmap.sh"

    // set input according to single_end parameter
    def input = params.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}" // if short reads check paired or not
    input = params.long_reads ?  "in=${fastq}" : input // else set as expected for long reads samples
   
    """
    ${tool} \\
        ref=$genome_index \\
        $input \\
        out=${fastq[0].baseName}.bam \\
        ${params.bbmap_options} \\
        threads=$task.cpus \\
        -Xmx${avail_mem}g &> ${fastq[0].baseName}.bbmap.log

    """
}