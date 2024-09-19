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
    """
    fileout=${sample.baseName.replace('.fastq','')}.result.txt
    seqkit convert -d ${sample} 2> \$fileout
    # get the result from las column of second row
    scoring=\$(awk  ' NR==2 { print \$( NF )  } ' \$fileout)
    # \${scoring,,} converts to lowercase
    if [[ \${scoring,,} == "sanger" || \${scoring,,} == "illumina-1.8+" ]];then  
        echo -e "\n keep intact" >> \$fileout
        # File passes the check, create a symlink to the input file
        ln -s ${sample} ${sample.baseName.replace('.fastq','')}_seqkit.fastq.gz 
    else
        echo -e "\n converted by seqkit" >> \$fileout
        seqkit convert ${sample} | gzip > ${sample.baseName.replace('.fastq','')}_seqkit.fastq.gz
    fi
    """


}
