/*
* To align with BWA-fastalign
* The index used must be the basename of the genome reference file
*/

process bwafastalign_index {
    label 'bwafastalign'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:
        """
        bwa-fastalign index $genome_fasta -p ${genome_fasta.baseName}
        """
}

/*
* To align with BWA FASTALIGN ALN
* The index used must be the basename of the genome reference file
*/
process bwafastalignaln {
    label 'bwafastalign'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*_sam.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwafastalign_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*_sam.log",  emit: bwafastalignaln_summary

    script:
        // options for bwafastalign-aln
        def bwafastalignaln_options = meta.bwafastalignaln_options ?: ""

        // catch output file prefix 
        def fileName = meta.uid + meta.suffix + "_bwafastalignaln"

        if (meta.paired){
        """
            bwa-fastalign aln ${bwafastalignaln_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} > ${reads[0].baseName}_r1.sai 2> ${reads[0].baseName}_r1_sai.log 
            bwa-fastalign aln ${bwafastalignaln_options} -t ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}_r2.sai 2> ${reads[1].baseName}_r2_sai.log
            bwa-fastalign sampe ${genome.baseName}  ${reads[0].baseName}_r1.sai ${reads[1].baseName}_r2.sai ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}_sam.log 
   
        """
        } else {
        """
            bwa-fastalign aln ${bwafastalignaln_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sai 2> ${fileName}_sai.log 
            bwa-fastalign samse ${genome.baseName}  ${fileName}.sai ${reads} > ${fileName}.sam 2> ${fileName}_sam.log 
        """
        }
}

/*
* To align with BWA FASTALIGN MEM
* The index used must be the basename of the genome reference file
*/
process bwafastalignmem {
    label 'bwafastalign'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*bwafastalign.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwafastalign_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwafastalignmem_summary

    script:
        // options for bwafastalign-mem
        def bwafastalignmem_options = meta.bwafastalignmem_options ?: ""

        // catch output file prefix 
        def fileName = meta.uid + meta.suffix + "_bwafastalignmem"
      
        if (meta.paired){
            """
            bwa-fastalign mem ${bwafastalignmem_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            bwa-fastalign mem ${bwafastalignmem_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}

/*
* To align with BWA FASTALIGN MEM
* The index used must be the basename of the genome reference file
*/
process bwafastalignsw {
    label 'bwafastalign'
    tag "${meta.uid}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path bwafastalign_index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwafastalignsw_summary

    script:
        // options for bwafastalign-mem
        def bwafastalignsw_options = meta.bwafastalignsw_options ?: ""

        // catch output file prefix 
        def fileName = meta.uid + meta.suffix + "_bwafastalignsw"

        if (meta.paired){
            """
            bwa-fastalign bwasw ${bwafastalignsw_options} -t ${task.cpus} ${genome.baseName} ${reads[0]} ${reads[1]} > ${fileName}.sam 2> ${fileName}.log 
            """
        } else {
            """
            bwa-fastalign bwasw ${bwafastalignsw_options} -t ${task.cpus} ${genome.baseName} ${reads} > ${fileName}.sam 2> ${fileName}.log 
            """
        }
}