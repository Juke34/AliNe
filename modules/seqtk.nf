/*
Here are described all processes related to seqtk
https://github.com/lh3/seqtk
*/

// A process for subsampling reads
process seqtk_sample {
    label 'seqtk'
    
    input:
        tuple val(id), path(fastq)
        
    output:
        tuple val(id), path("*_sampled.fastq.gz"), emit: sampled

    script:

        // set input/output according to short_paired parameter
        if (params.read_type == "short_paired"){
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${fastq[0].baseName.replace('.fastq','')}_sampled.fastq.gz
            seqtk sample -s100 ${fastq[1]} ${params.seqtk_sample_size}\\
                  > ${fastq[1].baseName.replace('.fastq','')}_sampled.fastq.gz
            """
        } else {
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${fastq[0].baseName.replace('.fastq','')}_sampled.fastq.gz
            """
        }
            
}