/* Module related to bbmmap
BBMap is a global aligner. That means it looks for the highest-scoring alignment taking into account all bases in a sequence. A local aligner would look for the best-scoring local alignment, meaning an alignment where the ends are possibly clipped off.

info:
To map quickly with very high precision and lower sensitivity, as when removing contaminant reads specific to a genome without risking false-positives:
bbmap.sh minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14
*/ 

process bbmap_index {
    label 'bbmap'
    tag "$genome_fasta"
    publishDir "${params.outdir}/${outpath}", mode: 'copy'

    input:
        path(genome_fasta)
        val outpath

    output:
        path('ref')

    script:
    
        // set memory to 4GB if not specified
        def avail_mem = task.memory ? task.memory.toGiga() : 4    
    
        """
        bbmap.sh ref=${genome_fasta} threads=${task.cpus} -Xmx${avail_mem}g
        """
}

process bbmap {
    label 'bbmap'
    tag "$sample"
    publishDir "${params.outdir}/${outpath}/stats", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(sample), path(fastq)
        path genome_index
        path hisat2_index_files
        val outpath

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam
        path "*.txt",  emit: bbmap_summary

    script:

    // set memory to 4GB if not specified
    def avail_mem = task.memory ? task.memory.toGiga() : 4 
    
    // set tool according to read length
    def tool = "bbmap.sh"
    if (params.read_type == "pacbio" || params.read_type == "ont"){
        tool = "mapPacBio.sh"
    }

    // set input according to read_type parameter
    def input =  "in=${fastq[0]}"
    if (params.read_type == "short_paired"){
        input =  "in=${fastq[0]} in2=${fastq[1]}" // if short reads check paired or not
    }

    // set fileName
    def fileName = fastq[0].baseName.replace('.fastq','')
    """
    ${tool} \\
        ref=$genome_index \\
        $input \\
        out=${fileName}.bam \\
        ${params.bbmap_options} \\
        threads=${task.cpus} \\
        bhist=${fileName}_bhist.txt qhist=${fileName}_qhist.txt aqhist=${fileName}_aqhist.txt lhist=${fileName}_lhist.txt ihist=${fileName}_ihist.txt \\
        ehist=${fileName}_ehist.txt qahist=${fileName}_qahist.txt indelhist=${fileName}_indelhist.txt mhist=${fileName}_mhist.txt \\
        gchist=${fileName}_gchist.txt idhist=${fileName}_idhist.txt scafstats=${fileName}_scafstats.txt \\
        -Xmx${avail_mem}g &> ${fileName}.bbmap.log.txt

    """
}