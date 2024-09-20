/*
 available values: 'sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'
*/

process seqkit_convert {
    label 'seqkit'
    tag "$sample_id"
    publishDir "${params.outdir}/${outpath}",  pattern: "*result.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(sample)
    val outpath

    output:
        tuple val(sample_id), path("*.fastq.gz"), emit: trimmed
        path("*result.txt")

    script: 

    if (params.read_type == "short_paired"){
    """
        fileout=${sample[0].baseName.replace('.fastq','')}.result.txt
    
        # run seqkit convert
        seqkit convert -d ${sample[0]} 2> \$fileout
        
        # get the result from las column of second row
        scoring=\$(awk  ' NR==2 { print \$( NF )  } ' \$fileout)

        # \${scoring,,} converts to lowercase
        if [[ \${scoring,,} == "sanger" || \${scoring,,} == "illumina-1.8+" ]];then  
            echo -e "\n keep intact" >> \$fileout
            # File passes the check, create a symlink to the input file
            ln -s ${sample[0]} ${sample[0].baseName.replace('.fastq','')}_seqkit.fastq.gz
            ln -s ${sample[1]} ${sample[1].baseName.replace('.fastq','')}_seqkit.fastq.gz
        else
            echo -e "\n converted by seqkit" >> \$fileout
            seqkit convert ${sample[0]} | gzip > ${sample[0].baseName.replace('.fastq','')}_seqkit.fastq.gz
            seqkit convert ${sample[1]} | gzip > ${sample[1].baseName.replace('.fastq','')}_seqkit.fastq.gz
        fi
        """
    } else {
    """
        fileout=${sample[0].baseName.replace('.fastq','')}.result.txt
    
        # run seqkit convert
        seqkit convert -d ${sample[0]} 2> \$fileout
        
        # get the result from las column of second row
        scoring=\$(awk  ' NR==2 { print \$( NF )  } ' \$fileout)
        
        # \${scoring,,} converts to lowercase
        if [[ \${scoring,,} == "sanger" || \${scoring,,} == "illumina-1.8+" ]];then  
            echo -e "\n keep intact" >> \$fileout
            # File passes the check, create a symlink to the input file
            ln -s ${sample} ${sample.baseName.replace('.fastq','')}_seqkit.fastq.gz 
        else
            echo -e "\n converted by seqkit" >> \$fileout
            seqkit convert ${sample} | gzip > ${sample[0].baseName.replace('.fastq','')}_seqkit.fastq.gz
        fi
        """
    }
}
