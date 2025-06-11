// ------------ nucmer -----------
// https://github.com/mummer4/mummer

/*
* To align with nucmer
*/
process nucmer {
    label 'mummer4'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*nucmer.log", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path genome
        val outpath

    output:
        tuple val(meta), path ("*.sam"), emit: tuple_sample_sam

    script:
        // options for nucmer
        def nucmer_options = meta.nucmer_options ?: ""

        // deal with library type - Not supported 

        // catch filename
        def fileName = AlineUtils.getCleanName(reads[0])

        // Name input
        reads0 = AlineUtils.getCleanName(reads[0])
        genomeReady = AlineUtils.getCleanName(genome) + "_nucmer.fa"

        if (meta.paired){
            reads1 = AlineUtils.getCleanName(reads[1])
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

            nucmer ${nucmer_options} -t ${task.cpus} ${genomeReady} ${reads0}  --sam-long ${fileName}_nucmer.sam
            nucmer ${nucmer_options} -t ${task.cpus} ${genomeReady} ${reads1}  --sam-long ${reads[1].baseName}_nucmer.sam

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
            nucmer ${nucmer_options} -t ${task.cpus} ${genomeReady} ${reads0}  --sam-long ${fileName}_nucmer.sam
            """
        }
}