process multiqc {
    label 'multiqc'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
        path log_files
        path multiqc_config

    output:
        path "*multiqc_report.html", optional:true
        path "*_data", optional:true

    script:
        """
        multiqc . -c ${multiqc_config}
        """
}