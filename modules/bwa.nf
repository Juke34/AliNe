process bwa_index {
    label 'bwa'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:
        """
        bwa index $genome_fasta -p ${genome_fasta.baseName}
        """
}

/*
* To align with BWA
* The index used must be the basename of the genome reference file
*/
process bwaaln {
    label 'bwa'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*bwaaln.sam"), emit: tuple_sample_sam
        path "*bwaaln_sam.log",  emit: bwaaln_summary

    script:
        // options for bwa-aln
        def bwaaln_options = meta.bwaaln_options ?: ""

        // catch filename
        fileName = reads[0].baseName.replaceAll(/\.(fastq|fq)$/, '')

        if (meta.paired){
        """
            bwa aln ${params.bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} > ${reads[0].baseName}_bwaaln_r1.sai 2> ${reads[0].baseName}_bwaaln_r1_sai.log 
            bwa aln ${params.bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}_bwaaln_r2.sai 2> ${reads[1].baseName}_bwaaln_r2_sai.log
            bwa sampe ${genome.baseName}  ${reads[0].baseName}_bwaaln_r1.sai ${reads[1].baseName}_bwaaln_r2.sai ${reads[0]} ${reads[1]} > ${fileName}_bwaaln.sam 2> ${fileName}_bwaaln_sam.log 
   
        """
        } else {
        """
            bwa aln ${params.bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwaaln.sai 2> ${fileName}_bwaaln_sai.log 
            bwa samse ${genome.baseName}  ${fileName}_bwaaln.sai ${reads} > ${fileName}_bwaaln.sam 2> ${fileName}_bwaaln_sam.log 
        """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwamem {
    label 'bwa'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*bwamem.sam"), emit: tuple_sample_sam
        path "*bwamem.log",  emit: bwamem_summary

    script:
        fileName = reads[0].baseName.replaceAll(/\.(fastq|fq)$/, '')
      
        if (meta.paired){
            """
            bwa mem ${params.bwamem_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}_bwamem.sam 2> ${fileName}_bwamem.log 
            """
        } else {
            """
            bwa mem ${params.bwamem_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwamem.sam 2> ${fileName}_bwamem.log 
            """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwasw {
    label 'bwa'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*bwasw.sam"), emit: tuple_sample_sam
        path "*bwasw.log",  emit: bwasw_summary

    script:

        fileName = reads[0].baseName.replaceAll(/\.(fastq|fq)$/, '')

        if (meta.paired){
            """
            bwa bwasw ${params.bwasw_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}_bwasw.sam 2> ${fileName}_bwasw.log 
            """
        } else {
            """
            bwa bwasw ${params.bwasw_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwasw.sam 2> ${fileName}_bwasw.log 
            """
        }
}