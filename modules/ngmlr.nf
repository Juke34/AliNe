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
        tuple val(sample), path(reads)
        path genome
        val outpath

    output:
        tuple val(sample), path ("*ngmlr.sam"), emit: tuple_sample_sam, optional:true
        path "*.log",  emit: ngmlr_summary

    script:

        fileName = reads[0].baseName.replace('.fastq','')

        """
        ngmlr ${params.ngmlr_options} -t ${task.cpus} -r ${genome} -q ${reads} -o ${fileName}_ngmlr.sam 2> ${fileName}_ngmlr.log 
        """

}