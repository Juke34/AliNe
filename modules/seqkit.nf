/*
 available values: 'sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'
*/

process seqkit_convert {
    label 'seqkit'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}",  pattern: "*result.txt", mode: 'copy'

    input:
    tuple val(meta), path(sample)
    val outpath

    output:
        tuple val(meta), path("*.fastq.gz"), emit: trimmed
        path("*result.txt")

    script: 
        def filebase0 = AlineUtils.getCleanName(sample[0])
        fileout = "${meta.id}_result.txt"

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
            ln -s \$(realpath ${sample[0]}) ${filebase0}_seqkit.fastq.gz
            ln -s \$(realpath ${sample[1]}) ${filebase1}_seqkit.fastq.gz
        else
            echo -e "\n converted by seqkit" >> ${fileout}
            seqkit convert ${sample[0]} | gzip > ${filebase0}_seqkit.fastq.gz
            seqkit convert ${sample[1]} | gzip > ${filebase1}_seqkit.fastq.gz
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
            ln -s \$(realpath ${sample}) ${sample.baseName.replace('.fastq','')}_seqkit.fastq.gz 
        else
            echo -e "\n converted by seqkit" >> ${fileout}
            seqkit convert ${sample} | gzip > ${sample[0].baseName.replace('.fastq','')}_seqkit.fastq.gz
        fi
        """
    }
}
