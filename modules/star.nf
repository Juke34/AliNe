// Will automatically append parameters if annotation is provided other nothing is done 
process prepare_star_index_options {
    label 'bash'

    input:
    tuple val(sample), path(reads)
    val annotation
    val read_length

    output:
        stdout()

    script:
        star_index_options = ""

        if (annotation){
            if (read_length){
               star_index_options = "--sjdbGTFfile ${annotation} --sjdbOverhang " + (read_length.toInteger() - 1)
                """
                #!/bin/bash
                echo -n "${star_index_options}"
                """
            } 
            else {
                file=reads[0] // In paired-end case we take the first file

                // Some bash commands as head or grep -q can interrupt another one, which return an error 141. To avoid the problem we unset pipefail
                """
                #!/bin/bash
                set +euo pipefail
                # use cat or zcat according to extension
                extension=\$(echo ${file} | awk -F . '{print \$NF}')
                if [[ \${extension} == "gz" ]];then
                    command="zcat"
                else
                    command="cat"
                fi
                a=0;b=0; RESULT=\$(\$command ${file} | head -n 40000 | awk '{if(NR%4==2){ b++; a+=length(\$1)}}END{print int(a/b)}')
                echo -n --sjdbGTFfile ${annotation} --sjdbOverhang \$((\${RESULT} - 1 ))
                """

            }  
        } 
        else {
            """
            #!/bin/bash
            echo -n "${star_index_options}"
            """
        }

}
            

process star_index {
    label 'star'
    tag "$genome_fasta"
    publishDir "${params.outdir}/star_indicies", mode: 'copy'

    input:
        path (genome_fasta)
        val star_index_options

    output:
        path ('*star_indicies')

    script:
        //log.info """ star_index_options: ${star_index_options} """
        """
        mkdir -p ${genome_fasta.baseName}_star_indicies
        STAR --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_fasta.baseName}_star_indicies \\
            --genomeFastaFiles  ${genome_fasta} \\
            ${star_index_options}
        """ 
}

process star {
    label 'star'
    tag "$sample"
    publishDir "${params.outdir}/star_alignments", pattern: "*.log", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path star_index

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_sam
        path "*.out",  emit: star_summary

    script:
        if (params.single_end){
        """
            STAR ${params.star_options} --genomeDir ${star_index} \\
                --readFilesIn ${reads}  \\
                --runThreadN ${task.cpus} \\
                --runMode alignReads \\
                --outFileNamePrefix ${reads.baseName}  \\
                --outSAMunmapped Within \\
                --outSAMtype BAM SortedByCoordinate
        """
        } else {
        """

        """
        }

}