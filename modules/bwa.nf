process bwa_index {
    label 'bwa'
    tag "$genome_fasta"
    publishDir "${params.outdir}/bwa_indicies", mode: 'copy'

    input:
        path(genome_fasta)

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
    tag "$sample"
    publishDir "${params.outdir}/bwa_alignments", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path genome
        path bwa_index_files

    output:
        tuple val(sample), path ("*bwaaln.sam"), emit: tuple_sample_sam
        path "*bwaaln_sam.log",  emit: bwaaln_summary

    script:
        fileName = reads[0].baseName
        if (params.single_end){
        """
            bwa aln ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwaaln.sai 2> ${fileName}_bwaaln_sai.log 
            bwa samse ${genome.baseName}  ${fileName}_bwaaln.sai ${reads} > ${fileName}_bwaaln.sam 2> ${fileName}_bwaaln_sam.log 
        """
        } else {
        """
            bwa aln ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} > ${reads[0].baseName}_bwaaln_r1.sai 2> ${reads[0].baseName}_bwaaln_r1_sai.log 
            bwa aln ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}_bwaaln_r2.sai 2> ${reads[1].baseName}_bwaaln_r2_sai.log
            bwa sampe ${genome.baseName}  ${reads[0].baseName}_bwaaln_r1.sai ${reads[1].baseName}_bwaaln_r1.sai ${reads[0]} ${reads[1]} > ${fileName}_bwaaln.sam 2> ${fileName}_bwaaln_sam.log 
   
        """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwamem {
    label 'bwa'
    tag "$sample"
    publishDir "${params.outdir}/bwa_alignments", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path genome
        path bwa_index_files

    output:
        tuple val(sample), path ("*bwamem.sam"), emit: tuple_sample_sam
        path "*bwamem.log",  emit: bwamem_summary

    script:
        fileName = reads[0].baseName
        if (params.single_end){
            """
            bwa mem ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwamem.sam 2> ${fileName}_bwamem.log 
            """
        } else {
            """
            bwa mem ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}_bwamem.sam 2> ${fileName}_bwamem.log 
            """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwasw {
    label 'bwa'
    tag "$sample"
    publishDir "${params.outdir}/bwa_alignments", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path genome
        path bwa_index_files

    output:
        tuple val(sample), path ("*bwasw.sam"), emit: tuple_sample_sam
        path "*bwasw.log",  emit: bwasw_summary

    script:
        fileName = reads[0].baseName
        if (params.single_end){
            """
            bwa bwasw ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}_bwasw.sam 2> ${fileName}_bwasw.log 
            """
        } else {
            """
            bwa bwasw ${params.bwa_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}_bwasw.sam 2> ${fileName}_bwasw.log 
            """
        }
}