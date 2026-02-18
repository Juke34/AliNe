// ------------ NGMLR -----------
// https://github.com/philres/ngmlr

/*
* To align with ngmlr
*/
process ngmlr {
    label 'ngmlr'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam, optional:true
        path "*.log",  emit: ngmlr_summary

    script:
        // options for ngmlr
        def ngmlr_options = meta.ngmlr_options ?: ""

        // catch output file prefix 
        def fileName = meta.file_id + meta.suffix + "_ngmlr"

        // For paired-end we concat output 
        if (meta.paired){
            """
            ngmlr ${ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[0]} -o ${reads[0].baseName}.sam 2> ${reads[0].baseName}_nglmr.log 
            ngmlr ${ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[1]} -o ${reads[1].baseName}.sam 2> ${reads[1].baseName}_nglmr.log 
            
            # Merge sam
            cat ${reads[0].baseName}.sam > ${fileName}_concatR1R2.sam
            rm ${reads[0].baseName}.sam
            awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}.sam >> ${fileName}_concatR1R2.sam
            rm ${reads[1].baseName}.sam
            """
        } else {
            """
            ngmlr ${ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[0]} -o ${fileName}.sam 2> ${fileName}.log 
            """
        }

}