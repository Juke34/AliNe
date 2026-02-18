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
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*_sam.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*_sam.log",  emit: bwaaln_summary

    script:
        // options for bwa-aln
        def bwaaln_options = meta.bwaaln_options ?: ""

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_bwaaln"

        if (meta.paired){
        """
            bwa aln ${bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} > ${reads[0].baseName}_r1.sai 2> ${reads[0].baseName}_r1_sai.log 
            bwa aln ${bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}_r2.sai 2> ${reads[1].baseName}_r2_sai.log
            bwa sampe ${genome.baseName}  ${reads[0].baseName}_r1.sai ${reads[1].baseName}_r2.sai ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}_sam.log 
   
        """
        } else {
        """
            bwa aln ${bwaaln_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sai 2> ${fileName}_sai.log 
            bwa samse ${genome.baseName}  ${fileName}.sai ${reads} > ${fileName}.sam 2> ${fileName}_sam.log 
        """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwamem {
    label 'bwa'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bwa.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwamem_summary

    script:
        // options for bwa-mem
        def bwamem_options = meta.bwamem_options ?: ""

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_bwamem"
      
        if (meta.paired){
            """
            bwa mem ${bwamem_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            bwa mem ${bwamem_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}

/*
* To align with BWA MEM
* The index used must be the basename of the genome reference file
*/
process bwasw {
    label 'bwa'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwa_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwasw_summary

    script:
        // options for bwa-mem
        def bwasw_options = meta.bwasw_options ?: ""

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_bwasw"

        if (meta.paired){
            """
            bwa bwasw ${bwasw_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            bwa bwasw ${bwasw_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}