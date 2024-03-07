#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

// Input/output params
params.reads = "/path/to/reads_{1,2}.fastq.gz"
params.genome = "/path/to/genome.fa"
params.outdir = "alignment_results"
params.reads_extension = ".fastq.gz" // Extension used to detect reads in folder
params.paired_reads_pattern = "_{1,2}"


// Read feature params
params.single_end = false
params.phred_score = "Phred+33"
params.stranded = false
params.strand_type = ""
params.read_length = ""
params.annotation = ""

// Aligner params
align_tools = [ 'bowtie2', 'bwaaln', 'bwamem', 'bwasw', 'graphmap2', 'hisat2', 'minimap2', 'nucmer', 'star' ]
params.aligner = ''
params.bowtie2_options = ''
params.bwa_options = ''
params.bwa_mem_options = ''
params.graphmap2_options = '' // owler option is possible
params.hisat2_options = ''
params.minimap2_options = ''
params.minimap2_index_options = '' //  -k, -w, -H and -I 
params.nucmer_options = ''
params.star_options = ''

// Report params
params.fastqc = false
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

//*************************************************
// STEP 1 - LOG INFO
//*************************************************

log.info """
IRD
.-./`) .-------.     ______
\\ .-.')|  _ _   \\   |    _ `''.
/ `-' \\| ( ' )  |   | _ | ) _  \\
 `-'`\"`|(_ o _) /   |( ''_'  ) |
 .---. | (_,_).' __ | . (_) `. |
 |   | |  |\\ \\  |  ||(_    ._) '
 |   | |  | \\ `'   /|  (_.\\.' /
 |   | |  |  \\    / |       .'
 '---' ''-'   `'-'  '-----'`


AliNe - Alignment in Nextflow
===================================================

General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     single_end                 : ${params.single_end}
     outdir                     : ${params.outdir}
  
Alignment Parameters
 bowtie2 parameters
     bowtie2_options            : ${params.bowtie2_options}
 bwa parameters
     bwa_options                : ${params.bwa_options}
 bwa-mem parameters
     bwa-mem_options            : ${params.bwa_mem_options}
 graphmap2 parameters
     graphmap2_options             : ${params.graphmap2_options}
 hisat2 parameters
     hisat2_options             : ${params.hisat2_options}
 minimap2  parameters
     minimap2_options           : ${params.minimap2_options}
 nucmer parameters
     nucmer_options            : ${params.nucmer_options}
 star parameters
     hisat2_options             : ${params.star_options}
 tophat2 parameters
     tophat2_options            : ${params.tophat2_options}

Report Parameters
 MultiQC parameters
     multiqc_config             : ${params.multiqc_config}

 """

//*************************************************
// STEP 2 - Include needed modules
//*************************************************

include {bowtie2_index; bowtie2} from "$baseDir/modules/bowtie2.nf"
include {bwa_index; bwaaln; bwamem; bwasw} from "$baseDir/modules/bwa.nf"
include {gaas_fastq_guessMyFormat} from "$baseDir/modules/gaas.nf"
include {graphmap2_index; graphmap2} from "$baseDir/modules/graphmap2.nf"
include {fastp} from "$baseDir/modules/fastp.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_ali} from "$baseDir/modules/fastqc.nf"
include {hisat2_index; hisat2} from "$baseDir/modules/hisat2.nf" 
include {minimap2_index; minimap2} from "$baseDir/modules/minimap2.nf" 
include {multiqc} from "$baseDir/modules/multiqc.nf" 
include {nucmer} from "$baseDir/modules/mummer4.nf" 
include { samtools_sam2bam as samtools_sam2bam_bowtie2; samtools_sam2bam as samtools_sam2bam_bwaaln; samtools_sam2bam as samtools_sam2bam_bwamem; samtools_sam2bam as samtools_sam2bam_bwasw; samtools_sam2bam as samtools_sam2bam_graphmap2; samtools_sam2bam as samtools_sam2bam_hisat2; samtools_sam2bam as samtools_sam2bam_minimap2; samtools_sam2bam as samtools_sam2bam_nucmer } from "$baseDir/modules/samtools.nf"
include { samtools_sort as samtools_sort_bowtie2;  samtools_sort as samtools_sort_bwaaln; samtools_sort as samtools_sort_bwamem; samtools_sort as samtools_sort_bwasw; samtools_sort as samtools_sort_graphmap2; samtools_sort as samtools_sort_hisat2; samtools_sort as samtools_sort_minimap2; samtools_sort as samtools_sort_nucmer } from "$baseDir/modules/samtools.nf"

include {prepare_star_index_options; star_index; star} from "$baseDir/modules/star.nf"

//*************************************************
// STEP 3 - Deal with parameters
//*************************************************

// Check aligner params. Can be a list (comma or space separated)
def aligner_list=[]
if( !params.aligner ){
    exit 1, "Error: <aligner> parameter is empty, please provide a aligner(s) among this list ${align_tools}.\n"
} else {
    str_list = params.aligner.tokenize(',')
    str_list.each {
        str_list2 = it.tokenize(' ')
        str_list2.each {
            if ( ! (it.toLowerCase() in align_tools) ){
                exit 1, "Error: <${it}> aligner not acepted, please provide aligner(s) among this list ${align_tools}.\n"
            }
            else{
                aligner_list.add(it.toLowerCase())
            }
        }
    }
}

// check profile
if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker')
  ) { "executer selected" }
else { exit 1, "No executer selected: -profile docker/singularity"}

// check input (file or folder?)
def list_files = []
def pattern_reads
def fromFilePairs_input
if (params.single_end) {
    pattern_reads = "${params.reads_extension}"
    fromFilePairs_input = params.reads_extension
} else {
    pattern_reads = "${fromFilePairs_input}${params.paired_reads_pattern}${params.reads_extension}"
    fromFilePairs_input = "${fromFilePairs_input}${pattern_reads}"
}
File input_reads = new File(params.reads)
if(input_reads.exists()){
    if ( input_reads.isDirectory()) {
       log.info "The input ${params.reads} is a folder!\n"
        //if (! input_reads.name.endsWith("/")) {
        //    params.reads = "${params.reads}/"
        //}
        input_reads.eachFileRecurse(FILES){
            if (it.name =~ ~/\${pattern_reads}/){
                list_files.add(it)
            }
        }
        samples_number = list_files.size()
        log.info "${samples_number} files in ${params.reads} with pattern ${pattern_reads}"

    }
    else {
        log.info "The input ${params.reads} is a file!\n"
        pattern_reads = "${params.reads}"
    }
} else {
    exit 1, "The input ${params.reads} does not exists!\n"
}

//*************************************************
// STEP 4 - Main Workflow
//*************************************************

workflow {

    main:
        Channel.fromFilePairs("${pattern_reads}", size: params.single_end ? 1 : 2, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find reads matching ${params.reads}!\n" }
            .set {reads}
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set {genome}
        align(reads,genome,aligner_list)
}


//*************************************************
// STEP 4 -  Workflow align
//*************************************************

workflow align {

    take:
        reads
        genome
        aligner_list

    main:

        Channel.empty().set{logs}
        
        // ------------------- Phred Score ----------------
        //if(! params.phred_score){
        //    gaas_fastq_guessMyFormat(reads)
        //}
        // ------------------- QC -----------------
        if(params.fastqc){
            fastqc(reads)
            logs.mix(fastqc.out)
        }
        
        // ------------------- BOWTIE2 -----------------
        if ("bowtie2" in aligner_list ){
            bowtie2_index(genome) // index
            bowtie2(reads, bowtie2_index.out.collect(), genome) // align
            logs.concat(bowtie2.out.bowtie2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie2(bowtie2.out.tuple_sample_sam)
            // sort
            samtools_sort_bowtie2(samtools_sam2bam_bowtie2.out.tuple_sample_bam)
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort_bowtie2.out.tuple_sample_bam, "ali_bowtie2")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }

        // ------------------- BWA -----------------
        if ("bwaaln" in aligner_list || "bwamem" in aligner_list || "bwasw" in aligner_list){
            bwa_index(genome) // index
            if ("bwaaln" in aligner_list){
                bwaaln(reads, genome, bwa_index.out.collect()) // align
                logs.concat(bwaaln.out.bwaaln_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwaaln(bwaaln.out.tuple_sample_sam)
                // sort
                samtools_sort_bwaaln(samtools_sam2bam_bwaaln.out.tuple_sample_bam)
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali(samtools_sort_bwaaln.out.tuple_sample_bam, "ali_bwaaln")
                    logs.concat(fastqc_ali.out).set{logs} // save log
                }
            }
            if ("bwamem" in aligner_list){
                bwamem(reads, genome, bwa_index.out.collect()) // align
                logs.concat(bwamem.out.bwamem_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwamem(bwamem.out.tuple_sample_sam)
                // sort
                samtools_sort_bwamem(samtools_sam2bam_bwamem.out.tuple_sample_bam)
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali(samtools_sort_bwamem.out.tuple_sample_bam, "ali_bwamem")
                    logs.concat(fastqc_ali.out).set{logs} // save log
                }
            }
            if ("bwasw" in aligner_list){
                bwasw(reads, genome, bwa_index.out.collect()) // align
                logs.concat(bwasw.out.bwasw_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwasw(bwasw.out.tuple_sample_sam)
                // sort
                samtools_sort_bwasw(samtools_sam2bam_bwasw.out.tuple_sample_bam)
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali(samtools_sort_bwasw.out.tuple_sample_bam, "ali_bwasw")
                    logs.concat(fastqc_ali.out).set{logs} // save log
                }
            }
        }

        // ------------------- GRAPHMAP2 -----------------
        if ("graphmap2" in aligner_list ){
            graphmap2_index(genome)
            graphmap2(reads, genome, graphmap2_index.out.collect()) // align
            logs.concat(graphmap2.out.graphmap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_graphmap2(graphmap2.out.tuple_sample_sam)
            // sort
            samtools_sort_graphmap2(samtools_sam2bam_graphmap2.out.tuple_sample_bam)
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort_graphmap2.out.tuple_sample_bam, "ali_graphmap2")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }

        // ------------------- HISAT2 -----------------
        if ("hisat2" in aligner_list){
            hisat2_index(genome) // index
            hisat2(reads, hisat2_index.out.collect()) // align
            logs.concat(hisat2.out.hisat2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_hisat2(hisat2.out.tuple_sample_sam)
            // sort
            samtools_sort_hisat2(samtools_sam2bam_hisat2.out.tuple_sample_bam)
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort_hisat2.out.tuple_sample_sortedbam, "ali_hisat2")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }

        // ------------------- minimap2 -----------------
        if ("minimap2" in aligner_list ){
            minimap2_index(genome) // index
            minimap2(reads, genome, minimap2_index.out.collect()) // align
            logs.concat(minimap2.out.minimap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_minimap2(minimap2.out.tuple_sample_sam)
            // sort
            samtools_sort_minimap2(samtools_sam2bam_minimap2.out.tuple_sample_bam)
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort_minimap2.out.tuple_sample_bam, "ali_minimap2")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }
        // ------------------- nucmer (mummer4) -----------------
        if ("nucmer" in aligner_list ){
             nucmer(reads, genome) // align
             // convert sam to bam
            samtools_sam2bam_nucmer(nucmer.out.tuple_sample_sam)
            // sort
            samtools_sort_nucmer(samtools_sam2bam_nucmer.out.tuple_sample_bam)
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort_nucmer.out.tuple_sample_bam, "ali_nucmer")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }

        // ------------------- STAR -----------------
        if ("star" in aligner_list){
            // Take read length information from only one sample for --sjdbOverhang option
            // only needed if --sjdbFileChrStartEnd or --sjdbGTFfile option is activated)
            prepare_star_index_options(reads, params.annotation, params.read_length)
            star_index(genome, prepare_star_index_options.out) // index
            star(reads, star_index.out.collect()) // align out is bam and sorted
            logs.concat(star.out.star_summary).set{logs} // save log
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali(samtools_sort.out.tuple_sample_sortedbam, "ali_star")
                logs.concat(fastqc_ali.out).set{logs} // save log
            }
        }

        multiqc(logs.collect(),params.multiqc_config)

}


/**************         onComplete         ****************/

workflow.onComplete {
    log.info ( workflow.success ? "\nAliNe pipeline complete!\n" : "Oops .. something went wrong\n" )

    log.info """
    AliNe Pipeline execution summary
    --------------------------------------
    Completed at : ${workflow.complete}
    UUID         : ${workflow.sessionId}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Exit Status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """
}

    // Move pipeline execution information into the result folder
    int num = 0;
    String save = "${params.pipeline_report}"
    File file = new File("${params.outdir}", save);
    while(file.exists()) {
        save = "${params.pipeline_report}" + "_" + (num++) ;
        file = new File("${params.outdir}", save);
    }
    Files.move(new File("${params.pipeline_report}").toPath(), file.toPath(), StandardCopyOption.REPLACE_EXISTING);
