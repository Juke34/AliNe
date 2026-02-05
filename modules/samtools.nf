/* Nucmer sam output is bugged see https://github.com/mummer4/mummer/issues/24
Here a special process to override this problem
*/
process samtools_sam2bam_nucmer {
    label 'samtools'
    tag "${meta.id}"

    input:
        tuple val(meta), file(sam)
        path genome

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_bam

    script:

        """
            cat ${sam} | sed 's/HD\\ /HD/' | sed 's/1.0\\ /1.0/' | sed 's/\tSO:coordinate/SO:coordinate/' | sed s'/VN1/VN:1/' | \
                sed 's/HD/HD\\t/' | sed 's/SO:unsorted/\\tSO:unsorted/' | sed 's/@PG /@PG\\t/' | sed 's/ PN/\\tPN/' | \
                sed 's/ VN/\\tVN/' | sed 's/ CL/\\tCL/' > ${sam.baseName}.fixed;
            samtools view -@ ${task.cpus} --reference  ${genome} -b  ${sam.baseName}.fixed -o ${sam.baseName}.fixed.sam
        """
      

}

process samtools_sam2bam {
    label 'samtools'
    tag "${meta.id}"

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta), path ("*.bam"), emit: tuple_sample_bam

    script:

        """
            samtools view -@ ${task.cpus} ${sam} -b -o ${sam.baseName}.bam 
        """

}
process samtools_merge_bam_if_paired {
    label 'samtools'
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path ("*.bam"), emit: tuple_sample_bam

    script:
        // define the output file name
        def output = AlineUtils.getCleanName(bam[0]) + "_concatR1R2.bam"

        if (meta.paired) {
            """
                samtools merge -@ ${task.cpus} ${output} *.bam
            """
        }
        // For single-end reads, we just link the bam file. The name is misleading but we need to pass the file to next process
        else {
            """
                ln -s \$(realpath ${bam}) ${output}
            """
        }
}
/*
http://www.htslib.org/doc/samtools-sort.html
Sort alignments by leftmost coordinates 
And convert to cram if needed (when samtools_bam2cram can be avoided) to save disk space
*/
process samtools_sort {
    label 'samtools'
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)
        path(genome_fasta)

    output:
        tuple val(meta), path ("*_sorted.{bam,cram}"), emit: tuple_sample_ali

    script:

        if (params.cram) {
            """
                samtools sort -@ ${task.cpus} --reference ${genome_fasta} -o ${bam.baseName}_sorted.cram ${bam}
            """
        } else {
            """
                samtools sort -@ ${task.cpus} -o ${bam.baseName}_sorted.bam ${bam}
            """
        }
}


/*
http://www.htslib.org/doc/samtools-view.html
Convert BAM to CRAM format
*/
process samtools_bam2cram {
    label 'samtools'
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)
        path(genome_fasta)

    output:
        tuple val(meta), path ("*.cram"), emit: tuple_sample_ali

    script:

        """
           samtools view -C -T ${genome_fasta} -o ${bam.baseName}.cram ${bam}
        """
}

/*
http://www.htslib.org/doc/samtools-index.html
Index BAM or CRAM files
*/
process samtools_index {
    label 'samtools'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy', pattern: "{*.bam,*.cram,*.crai,*.bai}"

    input:
        tuple val(meta), path(alignment)
        val outpath

    output:
        tuple val(meta), path(alignment), path ("*.{bai,crai}"), emit: tuple_sample_ali

    script:

        """
            samtools index ${alignment}
        """
}

/*
http://www.htslib.org/doc/samtools-stats.html
Produces comprehensive statistics from alignment file
*/
process samtools_stats {
    label 'samtools'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(bam), path(index)
        path(genome_fasta)
        val outpath
        val suffix

    output:
       path ("*.txt"), emit: bam_samtools_stats

    script:

        """
            samtools stats --reference ${genome_fasta} ${bam} > ${bam.baseName}.txt
        """
}
