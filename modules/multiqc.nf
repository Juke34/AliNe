process multiqc {
    label 'multiqc'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
        path log_files
        path multiqc_config

    output:
        path "*multiqc_report.html", optional:true, emit: multiqc_report_html
        path "*_data", optional:true, emit: multiqc_report_data

    script:
        """
        multiqc . -c ${multiqc_config}
        """
}