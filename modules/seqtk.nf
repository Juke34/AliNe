/*
Here are described all processes related to seqtk
https://github.com/lh3/seqtk
*/

// A process for subsampling reads
process seqtk_sample {
    label 'seqtk'
    tag "${meta.uid}"
    
    input:
        tuple val(meta), path(fastq)
        
    output:
        tuple val(meta), path("*_sampled.fastq.gz"), emit: sampled

    script:
        // add suffix to meta for output files
        def suffix = meta.suffix ? "${meta.suffix}_sampled" : '_sampled'
        meta = meta + [suffix: suffix]

        // set input/output according to short_paired parameter
        if (meta.paired){
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${meta.file_id[0]}${suffix}.fastq.gz
            seqtk sample -s100 ${fastq[1]} ${params.seqtk_sample_size}\\
                  > ${meta.file_id[1]}${suffix}.fastq.gz
            """
        } else {
            """
            seqtk sample -s100 ${fastq[0]} ${params.seqtk_sample_size}\\
                  > ${meta.file_id[0]}${meta.suffix}.fastq.gz
            """
        }
            
}