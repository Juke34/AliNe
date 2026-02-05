process r_rendering {
    label 'r_rendering'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
        path multiqc_data_dir

    output:
        path "alignment_comparison.tsv", emit: comparison_table_tsv

    script:
        """
        r_rendering.R -i ${multiqc_data_dir}/multiqc_general_stats.txt -o alignment_comparison.tsv
        """
}