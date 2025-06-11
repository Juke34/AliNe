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
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path minimap_index_files
        val outpath

    output:
        tuple val(meta), path ("*minimap2.sam"), emit: tuple_sample_sam, optional:true
        tuple val(meta), path ("*minimap2.paf"), emit: tuple_sample_paf, optional:true
        path "*minimap2.log",  emit: minimap2_summary

    script:
        // options for minimap2
        def minimap2_options = meta.minimap2_options ?: ""
        // catch filename
        def fileName = AlineUtils.getCleanName(reads)
        output_format = "paf"
        if ( minimap2_options.contains("-a") ){
            output_format = "sam"
        }
    
        if (meta.paired){
            """
            minimap2 ${minimap2_options} -t ${task.cpus} ${genome} ${reads[0]} ${reads[1]} > ${fileName}_minimap2.${output_format} 2> ${fileName}_minimap2.log 
            """
        } else {
            """
            minimap2 ${minimap2_options} -t ${task.cpus} ${genome} ${reads} > ${fileName}_minimap2.${output_format} 2> ${fileName}_minimap2.log 
            """
        }
}