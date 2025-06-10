/*
A parallel implementation of gzip for modern
multi-processor, multi-core machines
https://zlib.net/pigz/
*/
process fasta_uncompress {
    tag "$genome"
    label 'pigz'
 
    input:
        path(genome)

    output:
        path genomeFa, emit: genomeFa

    script:
        genomeReady = genome.baseName.replaceAll(/\.(gz)$/, '') // remove .gz
        genomeReady = genomeReady.replaceAll(/\.(fasta|fa)$/, '') // remove .fasta or .fa
        genomeFa = genomeReady + ".fa"
    """
        

        # DEALING WITH GZIPPED FILES to uncompress if needed
        extension=\$(echo ${genome} | awk -F . '{print \$NF}')

        if [[ \${extension} == "gz" ]];then
            pigz -dck -p ${task.cpus} ${genome} > ${genomeFa}
        elif [[ "${genome}" != "${genomeFa}" ]];then
            # link
            ln -s \$(realpath ${genome}) ${genomeFa}
        fi
    """
}