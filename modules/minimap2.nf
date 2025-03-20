// ------------ MINIMAP2 -----------
// https://github.com/lh3/minimap2


process minimap2_index {
    label 'minimap2'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath
        
    output:
        path("*")

    script:
        """
        minimap2 -d ${genome_fasta.baseName}.mmi ${genome_fasta}
        """
}

/*
* To align with mnimap2
* The index used must be the basename of the genome reference file
*/
process minimap2 {
    label 'minimap2'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(readtype), val(read_length)
        path genome
        path minimap_index_files
        val outpath

    output:
        tuple val(sample), path ("*minimap2.sam"), emit: tuple_sample_sam, optional:true
        tuple val(sample), path ("*minimap2.paf"), emit: tuple_sample_paf, optional:true
        path "*minimap2.log",  emit: minimap2_summary

    script:

        fileName = reads[0].baseName.replace('.fastq','')
        output_format = "paf"
        if ( params.minimap2_options.contains("-a") ){
            output_format = "sam"
        }
    
        if (params.read_type == "short_paired"){
            """
            minimap2 ${params.minimap2_options} -t ${task.cpus} ${genome} ${reads[0]} ${reads[1]} > ${fileName}_minimap2.${output_format} 2> ${fileName}_minimap2.log 
            """
        } else {
            """
            minimap2 ${params.minimap2_options} -t ${task.cpus} ${genome} ${reads} > ${fileName}_minimap2.${output_format} 2> ${fileName}_minimap2.log 
            """
        }
}