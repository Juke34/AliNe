// ------------ NGMLR -----------
// https://github.com/philres/ngmlr

/*
* To align with ngmlr
*/
process ngmlr {
    label 'ngmlr'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(readtype), val(read_length)
        path genome
        val outpath

    output:
        tuple val(sample), path ("*ngmlr.sam"), emit: tuple_sample_sam, optional:true
        path "*.log",  emit: ngmlr_summary

    script:
        // catch filename
        fileName = reads[0].baseName.replace('.fastq','')

        // For paired-end we concat output 
        if (params.read_type == "short_paired"){
            """
            ngmlr ${params.ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[0]} -o ${fileName}_ngmlr.sam 2> ${fileName}_ngmlr.log 
            ngmlr ${params.ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[1]} -o ${reads[1].baseName}_ngmlr.sam 2> ${fileName}_ngmlr.log 
            
            # Merge sam
            cat ${fileName}_ngmlr.sam > ${fileName}_ngmlr_concatR1R2.sam
            rm ${fileName}_ngmlr.sam
            awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}_ngmlr.sam >> ${fileName}_ngmlr_concatR1R2.sam
            rm ${reads[1].baseName}_ngmlr.sam
            """
        } else {
            """
            ngmlr ${params.ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads[0]} -o ${fileName}_ngmlr.sam 2> ${fileName}_ngmlr.log 
            """
        }

}