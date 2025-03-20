// ------------ nucmer -----------
// https://github.com/mummer4/mummer

/*
* To align with nucmer
*/
process nucmer {
    label 'mummer4'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}", pattern: "*nucmer.log", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(readtype), val(read_length)
        path genome
        val outpath

    output:
        tuple val(sample), path ("*.sam"), emit: tuple_sample_sam

    script:
        
        // deal with library type - Not supported 

        // extract the file name
        fileName = reads[0].baseName.replaceAll(/\.(fastq|fq)$/, '')
        // Name input
        reads0 = reads[0].baseName
        genomeReady = genome.baseName.replaceAll(/\.(fasta|fa)$/, '') + "_nucmer.fa"

        if (params.read_type == "short_paired"){
            reads1 = reads[1].baseName
            """
            # Prepare reference
            extension=\$(echo ${genome} | awk -F . '{print \$NF}')
            if [[ \${extension} == "gz" ]];then
                gzip -c -d ${genome} > ${genomeReady} 
            else
                ln -s ${genome} ${genomeReady} 
            fi
            # Prepare reads 0  
            extension=\$(echo ${reads[0]} | awk -F . '{print \$NF}')
            if [[ \${extension} == "gz" ]];then
                gzip -c -d ${reads[0]} > ${reads0}
            fi
            # Prepare reads 1
            extension=\$(echo ${reads[1]} | awk -F . '{print \$NF}')
            if [[ \${extension} == "gz" ]];then
                gzip -c -d ${reads[1]} > ${reads1}
            fi

            nucmer ${params.nucmer_options} -t ${task.cpus} ${genomeReady} ${reads0}  --sam-long ${fileName}_nucmer.sam
            nucmer ${params.nucmer_options} -t ${task.cpus} ${genomeReady} ${reads1}  --sam-long ${reads[1].baseName}_nucmer.sam

            # Merge sam
            cat ${fileName}_nucmer.sam > ${fileName}_nucmer_concatR1R2.sam
            rm ${fileName}_nucmer.sam
            awk '!/^@HD/ && !/^@SQ/ && !/^@RG/ && !/^@PG/ && !/^@CO/ && NF' ${reads[1].baseName}_nucmer.sam >> ${fileName}_nucmer_concatR1R2.sam
            rm ${reads[1].baseName}_nucmer.sam
            """
        } else{
            """
            # Prepare reference
            extension=\$(echo ${genome} | awk -F . '{print \$NF}')
            if [[ \${extension} == "gz" ]];then
                gzip -c -d ${genome} > ${genomeReady} 
            else
                ln -s ${genome} ${genomeReady} 
            fi

            # Prepare reads 0  
            extension=\$(echo ${reads[0]} | awk -F . '{print \$NF}')
            if [[ \${extension} == "gz" ]];then
                gzip -c -d ${reads[0]} > ${reads0}
            fi 
            nucmer ${params.nucmer_options} -t ${task.cpus} ${genomeReady} ${reads0}  --sam-long ${fileName}_nucmer.sam
            """
        }
}