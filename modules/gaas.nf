process gaas_fastq_guessMyFormat {
    label 'gaas'
    tag "$sample_id"
    publishDir "${params.outdir}/guessMyFormat",  pattern: "*result.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(sample)

    output:
    env(SCORETYPE)

    script: 
    """
    gaas_fastq_guessMyFormat.pl --fq U2OS_A1_R1_sub100000.fastq 2&> ${sample.baseName}.result.txt
    SCORETYPE=\$(awk '{if(NR==1){ print \$1}}' ${sample.baseName}.result.txt)
    """


}
