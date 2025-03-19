/* Nucmer sam output is bugged see https://github.com/mummer4/mummer/issues/24
Here a special process to override this problem
*/
process samtools_sam2bam_nucmer {
    label 'samtools'
    tag "$sample"

    input:
        tuple val(sample), file(sam)
        path genome

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_bam

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
    tag "$sample"

    input:
        tuple val(sample), path(sam)

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam

    script:

        """
            samtools view -@ ${task.cpus} ${sam} -b -o ${sam.baseName}.bam 
        """

}
process samtools_merge_bam {
    label 'samtools'
    tag "$sample"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam

    script:

        """
            samtools merge -@ ${task.cpus} ${bam[0].baseName}_concatR1R2.bam *.bam
        """

}
/*
http://www.htslib.org/doc/samtools-sort.html
Sort alignments by leftmost coordinates 
*/
process samtools_sort {
    label 'samtools'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(sample), path(bam)
        val outpath

    output:
        tuple val(sample), path ("*_sorted.bam"), emit: tuple_sample_sortedbam

    script:

        """
            samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam  
        """
}


/*
http://www.htslib.org/doc/samtools-stats.html
Produces comprehensive statistics from alignment file
*/
process samtools_stats {
    label 'samtools'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(sample), path(bam)
        path(genome_fasta)
        val outpath
        val suffix

    output:
       path ("*.txt"), emit: bam_samtools_stats

    script:

        """
            samtools stats  --reference ${genome_fasta} ${bam} > ${bam.baseName}.txt
        """
}
