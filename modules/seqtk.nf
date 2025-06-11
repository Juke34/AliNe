/*
Here are described all processes related to seqtk
https://github.com/lh3/seqtk
*/

// A process for subsampling reads
process seqtk_sample {
    label 'seqtk'
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(fastq)
        
    output:
        tuple val(meta), path("*_sampled.fastq.gz"), emit: sampled

    script:

        // get the output base name
        def baseOutFile1 = AlineUtils.getCleanName(fastq[0])

        // set input/output according to short_paired parameter
        if (meta.paired){
            def baseOutFile2 = AlineUtils.getCleanName(fastq[1])
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${baseOutFile1}_sampled.fastq.gz
            seqtk sample -s100 ${fastq[1]} ${params.seqtk_sample_size}\\
                  > ${baseOutFile2}_sampled.fastq.gz
            """
        } else {
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${baseOutFile1}_sampled.fastq.gz
            """
        }
            
}