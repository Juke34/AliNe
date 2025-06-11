/* Module related to Novoalign
Novoalign is highly accurate program for mapping next-generation sequencing reads to a reference database.
https://www.novocraft.com/documentation/novoalign-2/
*/ 

process novoalign_index {
    label 'novoalign'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('*nix')

    script:
    
        // set memory to 4GB if not specified
        def avail_mem = task.memory ? task.memory.toGiga() : 4    
    
        """
        novoindex ${genome_fasta.baseName}.nix ${genome_fasta}
        """
}

process novoalign {
    label 'novoalign'
    tag "${meta.id}"
    containerOptions "${novoalign_lic}"
    publishDir "${params.outdir}/${outpath}/stats", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(meta), path(fastq), val(library), val(read_length)
        path genome
        path genome_index
        val novoalign_lic
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.txt",  emit: novoalign_log

    script:
        // options for novoalign
        def novoalign_options = meta.novoalign_options ?: ""

        // set input according to short_paired parameter
        def input = meta.paired ? "${fastq[0]} ${fastq[1]}" : "${fastq}" // if short reads check paired or not

        // deal with library type - Not supported

        // set fileName
        def fileName = AlineUtils.getCleanName(fastq[0])
        """
        novoalign ${novoalign_options} \\
            -d ${genome_index} \\
            -f ${input} \\
            -o SAM > ${fileName}.sam 2> log.txt

        """
}