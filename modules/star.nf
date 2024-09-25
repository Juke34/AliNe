// Will automatically append parameters if annotation is provided other nothing is done 
process prepare_star_index_options {
    label 'bash'

    input:
    tuple val(sample), path(reads), val(library)
    path annotation
    val read_length

    output:
        stdout()

    script:
        star_index_options = "${params.star_index_options}"

        if ( annotation.toString() != "null.gtf" ){
            // Deal with --sjdbGTFfile option
            if( !star_index_options.contains("--sjdbGTFfile") ){
                star_index_options += " --sjdbGTFfile ${annotation}"
            }
            // Deal with --sjdbOverhang option
            if (read_length && !star_index_options.contains("--sjdbOverhang") ){
               star_index_options += " --sjdbOverhang " + (read_length.toInteger() - 1)
                """
                #!/bin/bash
                echo -n "${star_index_options}"
                """
            } 
            else if ( !star_index_options.contains("--sjdbOverhang") ) {
                file=reads[0] // In paired-end case we take the first file
                log.info """Expected read length parameter not provided (--read_length) for star, using the first file ${reads} from sample ${sample} as reference to deduce this value..."""

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
                echo -n ${star_index_options} --sjdbOverhang \$((\${RESULT} - 1 ))
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
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path (genome_fasta)
        val star_index_options
        path annotation // needed in case the star_index_options contains --sjdbGTFfile option
        val outpath

    output:
        path ('*star_indicies')

    script:
        
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
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library)
        path star_index
        path annotation
        val outpath

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam
        path "*.out",  emit: star_summary
        path "*SJ.out.tab", emit: splice_junctions

    script:

        // set tool according to read length
        def star_tool = "STAR"
        if (params.read_type == "pacbio" || params.read_type == "ont"){
            star_tool = "STARlong"
        }

        // deal with library type - Not supported 

        if (params.read_type == "short_paired"){
        """
            mkfifo pipedRead1
            zcat < ${reads[0]} > pipedRead1 &
            mkfifo pipedRead2
            zcat < ${reads[1]} > pipedRead2 &

            ${star_tool} ${params.star_options} --genomeDir ${star_index} \\
                --readFilesIn pipedRead1 pipedRead2 \\
                --runThreadN ${task.cpus} \\
                --runMode alignReads \\
                --outFileNamePrefix ${reads[0].baseName.replace('.fastq','')}  \\
                --outSAMunmapped Within \\
                --outSAMtype BAM SortedByCoordinate
        """
        } else {
        """
            # the fifo setup (or any other "transient" files) will not work with the 2-pass mode, as STAR need to read the files twice
            mkfifo pipedRead
            zcat < ${reads} > pipedRead &

            ${star_tool} ${params.star_options} --genomeDir ${star_index} \\
                --readFilesIn pipedRead  \\
                --runThreadN ${task.cpus} \\
                --runMode alignReads \\
                --outFileNamePrefix ${reads.baseName.replace('.fastq','')}  \\
                --outSAMunmapped Within \\
                --outSAMtype BAM SortedByCoordinate
        """
        }
}

/* 
your goal is to robustly and accurately identify novel splice junctions for differential splicing analysis and variant discovery, it is highly recommended to use the STAR in 2-pass mode.
In the two-pass mode, the genome indices are re-generated from splice junctions obtained from a 1-pass mode with the usual parameters and then run the mapping step (2-pass mapping). 
Optionally, you can also bypass genome indices re-generation and directly provide a 1-pass splice junction file during the 2-pass mapping step.
For a study with multiple samples, it is recommended to collect 1st pass junctions from all samples.
*/
process star2pass{
    label 'star'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(library)
        path star_index
        path splice_junctions
        path annotation
        val outpath

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam
        path "*.out",  emit: star_summary

    script:

        // set tool according to read length
        def star_tool = "STAR"
        if (params.read_type == "pacbio" || params.read_type == "ont"){
            star_tool = "STARlong"
        }

        // deal with library type - Not supported in STAR

        if (params.read_type == "short_paired"){
        """
            mkfifo pipedRead1
            zcat < ${reads[0]} > pipedRead1 &
            mkfifo pipedRead2
            zcat < ${reads[1]} > pipedRead2 &

            # run STAR
            ${star_tool} ${params.star_options} --genomeDir ${star_index} \\
                --readFilesIn pipedRead1 pipedRead2  \\
                --runThreadN ${task.cpus} \\
                --runMode alignReads \\
                --outFileNamePrefix ${reads.baseName.replace('.fastq','')}_2pass  \\
                --outSAMunmapped Within \\
                --outSAMtype BAM SortedByCoordinate \\
                --sjdbFileChrStartEnd *SJ.out.tab
        """
        } else {
        """
            # Select cat or zcat according to suffif
            suffix=\$(echo ${reads[0]} | awk -F . '{print \$NF}')
            command=cat
            if [[ \${suffix} == "gz" ]];then
                command=zcat
            fi
            
            # the fifo setup (or any other "transient" files) will not work with the 2-pass mode, as STAR need to read the files twice
            mkfifo pipedRead
            command < ${reads} > pipedRead &

            # run STAR
            ${star_tool} ${params.star_options} --genomeDir ${star_index} \\
                --readFilesIn pipedRead  \\
                --runThreadN ${task.cpus} \\
                --runMode alignReads \\
                --outFileNamePrefix ${reads[0].baseName.replace('.fastq','')}_2pass  \\
                --outSAMunmapped Within \\
                --outSAMtype BAM SortedByCoordinate \\
                --sjdbFileChrStartEnd *SJ.out.tab
        """
        }
}

/*
Info: 
--outFilterMismatchNmax 100 : increases the number of allowed mismatches to 100 - need to allow more mismatches for longer reads
--seedSearchStartLmax 30 : increases the number of seed search start position in the read - important for reads with high error rate
--seedSearchLmax 30 : similar to the above, limits the maximum length of the seeds. Presently, I do not recommend changing this parameter
--seedPerReadNmax 100000   --seedPerWindowNmax 100 : increase the number of allowed seeds for each read and alignment window - need to store more seeds for longer reads
--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 : increase the number of allowed alignments for each read and alignment window - need to store more putative alignments for longer reads
*/