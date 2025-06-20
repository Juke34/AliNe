// ------------ LAST -----------
// https://gitlab.com/mcfrith/last/

process last_index {
    label 'last'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path("*")

    script:

        """
        lastdb -P ${task.cpus} -uNEAR -R01 ${genome_fasta.baseName} $genome_fasta 
        """
}

process last {
    label 'last'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*last.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        path index_files
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam
        path "*.log",  emit: last_summary

    script:
        // options for last
        def last_options = meta.last_options ?: ""

        // catch filename
        def fileName = AlineUtils.getCleanName(reads)

        // For paired-end we concat output 
        if (meta.paired){
            """
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[0]} > ${fileName}_last.maf 2> ${fileName}_last.log
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}_last.maf 2> ${fileName}_last.log
            
            # convert to sam
            maf-convert -d sam ${fileName}_last.maf > ${fileName}_last.sam
            maf-convert -d sam ${reads[1].baseName}_last.maf > ${reads[1].baseName}_last.sam

            # Merge sam
            cat ${fileName}_last.sam > ${fileName}_last_concatR1R2.sam
            rm ${fileName}_last.sam
            awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}_last.sam >> ${fileName}_last_concatR1R2.sam
            rm ${reads[1].baseName}_last.sam
            """
        } else {
            """
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[0]} > ${fileName}_last.maf 2> ${fileName}_last.log

            # convert to sam
            maf-convert -d sam ${fileName}_last.maf > ${fileName}_last.sam
            """
        }

}