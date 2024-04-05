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
        tuple val(sample), path(reads)
        path genome
        val outpath

    output:
        tuple val(sample), path ("*nucmer.sam"), emit: tuple_sample_sam

    script:


        fileName = reads[0].baseName
        int paired=0
        if (reads[1]){ paired = 1}
        """
        # Prepare reference
        genomeReady="${genome.baseName}_nucmer.fa"
        extension=\$(echo ${genome} | awk -F . '{print \$NF}')
        if [[ \${extension} == "gz" ]];then
             gzip -c -d ${genome} > \$genomeReady 
        else
            ln -s ${genome} \$genomeReady 
        fi
        # Prepare reads
        readsReady="${fileName}.fq"
        extension=\$(echo ${reads[0]} | awk -F . '{print \$NF}')
        if [[ \${extension} == "gz" ]];then
            gzip -c -d ${reads[0]} > \$readsReady 
            # paired reads - we append the fastq
            if [[ ${paired} > 1 ]];then
                gzip -c -d ${reads[1]} >> \$readsReady 
            fi
        else
            cp ${reads[0]} \$readsReady 
            if [[ ${paired} > 1 ]];then
                cat ${reads[1]} >> \$readsReady 
            fi
        fi

        nucmer ${params.nucmer_options} -t ${task.cpus} \$genomeReady \$readsReady  --sam-long ${fileName}nucmer.sam
        """
        


}