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
    tag "${meta.uid}"
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

        // catch output file prefix 
        def fileName = meta.uid + meta.suffix + "_last"

        // For paired-end we concat output 
        if (meta.paired){
            """
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[0]} > ${reads[0].baseName}.maf 2> ${reads[0].baseName}.log
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[1]} > ${reads[1].baseName}.maf 2> ${reads[1].baseName}.log
            
            # convert to sam
            maf-convert -d sam ${reads[0].baseName}.maf > ${fileName}.sam
            maf-convert -d sam ${reads[1].baseName}.maf > ${reads[1].baseName}.sam

            # Merge sam
            cat ${fileName}.sam > ${fileName}_concatR1R2.sam
            rm ${fileName}.sam
            awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}.sam >> ${fileName}_concatR1R2.sam
            rm ${reads[1].baseName}.sam
            """
        } else {
            """
            lastal ${last_options} -Q 1 -P ${task.cpus} ${genome.baseName} ${reads[0]} > ${fileName}.maf 2> ${fileName}.log

            # convert to sam
            maf-convert -d sam ${fileName}.maf > ${fileName}.sam
            """
        }

}