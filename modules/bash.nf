/*
Here are described all processes related to bash
*/

// A process to compute the mean read length of a FASTQ
process read_length {
    label 'bash'
    publishDir "${params.outdir}/${outpath}", pattern: "*", mode: 'copy'

    input:
        tuple val(id), path(fastq)
        val outpath
        
    output:
        tuple val(id), env(READLENGTH), emit: tuple_id_readlength
        path "*read_length.txt"

    script:
        // compute made only on first file when is paired
        """
        READLENGTH=\$(cat ${fastq[0]} | awk 'NR % 4 == 2 {sum += length(\$0); count++} END {if (count > 0) print int((sum / count) + 0.5); else print "0"}')
        echo \$READLENGTH > ${fastq[0].baseName.replace('.fastq','')}_read_length.txt
        """
}


// Process to set the read length of a FASTQ from params
process set_tuple_withUserReadLength{
    label 'bash'
   
    input:
        tuple val(id), path(fastq)

    output:
        tuple val(id), val(params.read_length), emit: tuple_id_libtype

   
    script:

        """
        """
}