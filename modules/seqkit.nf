/*
 available values: 'sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'
*/

process seqkit_convert {
    label 'seqkit'
    tag "${meta.file_id}"
    publishDir "${params.outdir}/${outpath}",  pattern: "*result.txt", mode: 'copy'

    input:
    tuple val(meta), path(sample)
    val outpath

    output:
        tuple val(meta), path("*.fastq.gz"), emit: trimmed
        path("*result.txt")

    script: 
        // add suffix to meta for output files
        def suffix = meta.suffix ? "${meta.suffix}_sekqit" : '_sekqit'
        meta = meta + [suffix: suffix]

        // catch filename
        def filebase0 = AlineUtils.getCleanName(sample[0])
        fileout = "${meta.file_id}${suffix}_result.txt"

    if (meta.paired){
        def filebase1 = AlineUtils.getCleanName(sample[1])
    """
        # run seqkit convert
        seqkit convert -d ${sample[0]} 2> ${fileout}
        
        # get the result from las column of second row
        scoring=\$(awk  ' NR==2 { print \$( NF )  } ' ${fileout})

        # \${scoring,,} converts to lowercase
        if [[ \${scoring,,} == "sanger" || \${scoring,,} == "illumina-1.8+" ]];then  
            echo -e "\n keep intact" >> ${fileout}
            # File passes the check, create a symlink to the input file
            ln -s \$(realpath ${sample[0]}) ${filebase0}${suffix}.fastq.gz
            ln -s \$(realpath ${sample[1]}) ${filebase1}${suffix}.fastq.gz
        else
            echo -e "\n converted by seqkit" >> ${fileout}
            seqkit convert ${sample[0]} | gzip > ${filebase0}${suffix}.fastq.gz
            seqkit convert ${sample[1]} | gzip > ${filebase1}${suffix}.fastq.gz
        fi
        """
    } else {
    """
        # run seqkit convert
        seqkit convert -d ${sample[0]} 2> ${fileout}
        
        # get the result from las column of second row
        scoring=\$(awk  ' NR==2 { print \$( NF )  } ' ${fileout})
        
        # \${scoring,,} converts to lowercase
        if [[ \${scoring,,} == "sanger" || \${scoring,,} == "illumina-1.8+" ]];then  
            echo -e "\n keep intact" >> ${fileout}
            # File passes the check, create a symlink to the input file
            ln -s \$(realpath ${sample}) ${filebase0}${suffix}.fastq.gz 
        else
            echo -e "\n converted by seqkit" >> ${fileout}
            seqkit convert ${sample} | gzip > ${filebase0}${suffix}.fastq.gz
        fi
        """
    }
}

/*
 * Clean FASTA headers by removing everything after the first space
 * and create samtools index
 */
process seqkit_clean_fasta_headers {
    label 'seqkit'
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(fasta)
        val outpath

    output:
        path("*_clean.fa"), emit: clean_fasta

    script:
        def filename = AlineUtils.getCleanName(fasta)
        """
            seqkit replace -p " .*" -r "" ${fasta} > ${filename}_clean.fa
        """
}
