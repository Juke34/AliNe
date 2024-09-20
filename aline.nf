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
params.paired_reads_pattern = "{1,2}"
read_type_allowed = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type = "short_paired"
params.relax = false // allow to use options that do not reflect expectation according to the read type

// Read feature params
params.stranded = false 
params.strand_type = "" // xxx Not used yet
params.read_length = ""
params.annotation = ""

// Trimming params
params.trimming_fastp = false

// Aligner params
align_tools = [ 'bbmap', 'bowtie2', 'bwaaln', 'bwamem', 'bwasw', 'graphmap2', 'hisat2', 'minimap2', 'novoalign', 'nucmer', 'ngmlr', 'star', 'subread' ]
params.aligner = ''
params.bbmap_options = ''
params.bowtie2_options = ''
params.bwaaln_options = ''
params.bwamem_options = ''
params.bwasw_options = ''
params.graphmap2_options = '' // owler option is possible
params.hisat2_options = ''
params.minimap2_options = '' 
params.minimap2_index_options = '' //  -k, -w, -H and -I
params.ngmlr_options = ''
params.novoalign_options = ''
params.nucmer_options = ''
params.novoalign_license = '' // license. You can ask for one month free trial license at http://www.novocraft.com/products/novoalign/
params.tophat2_options = ''
params.star_options = ''
params.star_2pass = false
params.subread_options = '-t 0'// -t specifes the type of input sequencing data. Possible values include 0, denoting RNA-seq data, or 1, denoting genomic DNA-seq data.

// Report params
params.fastqc = false
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

// other
params.help = null

//*************************************************
// STEP 1 - HELP
//*************************************************

log.info header()
if (params.help) { exit 0, helpMSG() }

log.info """
General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     aligner                    : ${params.aligner}
     read_type                  : ${params.read_type}
     paired_reads_pattern       : ${params.paired_reads_pattern}
     reads_extension            : ${params.reads_extension}
     stranded                   : ${params.stranded}
     outdir                     : ${params.outdir}

Report Parameters
 MultiQC parameters
     fastqc                     : ${params.fastqc}
     multiqc_config             : ${params.multiqc_config}

Alignment Parameters before check
"""
log.info printAlignerOptions()

//*************************************************
// STEP 1 - PARAMS CHECK
//*************************************************
log.info """check aligner provided: ${params.aligner} ..."""
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

// check read type parameter
log.info """check read type parameter: ${params.read_type} ..."""
if( ! params.read_type ){
    exit 1, "Error: <read_type> parameter is empty, please provide a read type among this list ${read_type_allowed}.\n"
} else {
    if ( ! (params.read_type.toLowerCase() in read_type_allowed) ){
        exit 1, "Error: <${params.read_type}> aligner not acepted, please provide a read type among this list ${read_type_allowed}.\n"
    }
}

// ----------------------------------------------------------

log.info """check alinger parameters ..."""
def stop_pipeline = false
// --- bbmap tool ---
def bbmap_tool = "bbmap.sh"
if ("bbmap" in aligner_list && !params.relax){
    if (params.read_type == "pacbio" || params.read_type == "ont"){
        bbmap_tool = "mapPacBio.sh"
        // Function to check and set maxlen in params.bbmap_options when long_reads is set
        // params is supposed to be a immutable. Using params.replace method might not be supported in the future 
        if ( !params.bbmap_options.contains("maxlen") ){
            params.replace("bbmap_options", "${params.bbmap_options} maxlen=5000")
        }
    }
}

// ---- minimap2 tool ---
// Force -a option to be sure to get sam output
if ("minimap2" in aligner_list && !params.relax){

    if (params.read_type == "short_single" || params.read_type == "short_paired"){
        if ( ! params.minimap2_options.contains("--sr ") ){
            params.replace("minimap2_options", "--sr ${params.minimap2_options}")
        }
    }
    if (params.read_type == "pacbio"){
        if ( ! params.minimap2_options.contains(" ava-pb") and ! params.minimap2_options.contains(" splice:hq") and 
            ! params.minimap2_options.contains(" map-hifi") and ! params.minimap2_options.contains(" map-pb") ){
            log.warn("""Error: <${params.minimap2_options}> minimap2 options not accepted for ont data, please provide options among this list, ava-pb, splice:hq, map-hifi, map-pb (see https://github.com/lh3/minimap2).
Otherwise, if you know what you are doing you can activate the aline --relax parameter to use options that do not reflect expectation.\n""")
            stop_pipeline = true
        }
    }
    if (params.read_type == "ont"){
        if ( ! params.minimap2_options.contains(" ava-ont") and ! params.minimap2_options.contains(" splice") and 
            ! params.minimap2_options.contains(" lr:hq") and ! params.minimap2_options.contains(" map-ont") ){
            log.warn("""Error: <${params.minimap2_options}> minimap2 options not accepted for ont data, please provide options among this list, ava-ont, splice, lr:hq, map-ont (see https://github.com/lh3/minimap2).
Otherwise, if you know what you are doing you can activate the aline --relax parameter to use options that do not reflect expectation.\n""")
            stop_pipeline = true
        }
    }
}
// relax or not this option has to be used
if ( ! params.minimap2_options.contains("-a ") ){
     params.replace("minimap2_options", "${params.minimap2_options} -a")
}

// ngmlr tool - check options
if ("ngmlr" in aligner_list ){
    if (!params.relax) {
        // for pacbio reads, set -g 20 and -x 0
        if (params.read_type == "ont"){
            if (! params.ngmlr_options.contains("-x ont")){
                params.replace("ngmlr_options", "${params.ngmlr_options} -x ont")
            }
        }
        else if (params.read_type.contains == "short"){
            log.warn ": ngmlr aligner do not handle ont short reads, please remove it from the list of aligner to use.\nOtherwise, if you know what you are doing you can activate the aline --relax parameter to use options that do not reflect expectation.\n"
            stop_pipeline = true
        }
    }
}       

// novoalign tool - load license into the container
if ("novoalign" in aligner_list ){
    def novoalign_container_options = ""
    if( params.novoalign_license ){
        novoalign_container_options = "-v ${params.novoalign_license}:/usr/local/bin/"
    } else {
        log.warn ": NovoAlign aligner selected but no license provided. Please provide a license to run novoalign.\n"
        stop_pipeline = true
    }

    if (!params.relax) {
        // for pacbio reads, set -g 20 and -x 0
        if (params.read_type == "pacbio"){
            if (! params.novoalign_options.contains("-g ")){
                params.replace("novoalign_options", "${params.novoalign_options} -g 20")
            }
            if (! params.novoalign_options.contains("-x ")){
                params.replace("novoalign_options", "${params.novoalign_options} -x 0")
            }
        }
        if ( params.read_type == "ont" ){
            log.warn ": NovoAlign aligner do not handle ont data, please remove it from the list of aligner to use.\nOtherwise, if you know what you are doing you can activate the aline --relax parameter to use options that do not reflect expectation.\n"
            stop_pipeline = true
        }
    }
}

if(stop_pipeline){
    exit 1, "Please fix previous issues in order to run the pipeline.\n"
}

//*************************************************
// STEP 1 - LOG INFO
//*************************************************

log.info """\nAlignment Parameters after check"""
log.info printAlignerOptions()

//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include {bbmap_index; bbmap} from "$baseDir/modules/bbmap.nf"
include {bowtie2_index; bowtie2} from "$baseDir/modules/bowtie2.nf"
include {bwa_index; bwaaln; bwamem; bwasw} from "$baseDir/modules/bwa.nf"
include {seqkit_convert} from "$baseDir/modules/seqkit.nf"
include {graphmap2_index; graphmap2} from "$baseDir/modules/graphmap2.nf"
include {fastp} from "$baseDir/modules/fastp.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_fastp; fastqc as fastqc_ali_bbmap; fastqc as fastqc_ali_bowtie2; 
         fastqc as fastqc_ali_bwaaln; fastqc as fastqc_ali_bwamem; fastqc as fastqc_ali_bwasw; fastqc as fastqc_ali_graphmap2; 
         fastqc as fastqc_ali_hisat2; fastqc as fastqc_ali_minimap2; fastqc as fastqc_ali_ngmlr; fastqc as fastqc_ali_novoalign; 
         fastqc as fastqc_ali_nucmer; fastqc as fastqc_ali_star; fastqc as fastqc_ali_subread } from "$baseDir/modules/fastqc.nf"
include {hisat2_index; hisat2} from "$baseDir/modules/hisat2.nf" 
include {minimap2_index; minimap2} from "$baseDir/modules/minimap2.nf" 
include {multiqc} from "$baseDir/modules/multiqc.nf" 
include {ngmlr} from "$baseDir/modules/ngmlr.nf" 
include {nucmer} from "$baseDir/modules/mummer4.nf" 
include {novoalign_index; novoalign} from "$baseDir/modules/novoalign.nf"
include {samtools_sam2bam_nucmer; samtools_sam2bam as samtools_sam2bam_bowtie2; samtools_sam2bam as samtools_sam2bam_bwaaln; 
         samtools_sam2bam as samtools_sam2bam_bwamem; samtools_sam2bam as samtools_sam2bam_bwasw; samtools_sam2bam as samtools_sam2bam_graphmap2; 
         samtools_sam2bam as samtools_sam2bam_hisat2; samtools_sam2bam as samtools_sam2bam_minimap2; 
         samtools_sam2bam as samtools_sam2bam_ngmlr; samtools_sam2bam as samtools_sam2bam_novoalign } from "$baseDir/modules/samtools.nf"
include {samtools_sort as samtools_sort_bbmap; samtools_sort as samtools_sort_bowtie2; samtools_sort as samtools_sort_bwaaln; 
         samtools_sort as samtools_sort_bwamem; samtools_sort as samtools_sort_bwasw; samtools_sort as samtools_sort_graphmap2; 
         samtools_sort as samtools_sort_hisat2; samtools_sort as samtools_sort_minimap2; samtools_sort as samtools_sort_ngmlr; 
         samtools_sort as samtools_sort_novoalign;  samtools_sort as samtools_sort_nucmer} from "$baseDir/modules/samtools.nf"
include {subread_index; subread} from "$baseDir/modules/subread.nf"
include {prepare_star_index_options; star_index; star; star2pass} from "$baseDir/modules/star.nf"

//*************************************************
// STEP 3 - CHECK 2 for parameters
//*************************************************

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
def per_pair = false // read the reads per pair
def path_reads = params.reads 

// in case of folder provided, add a trailing slash if missing
File input_reads = new File(path_reads)
if(input_reads.exists()){
    if ( input_reads.isDirectory()) {
        if (! input_reads.name.endsWith("/")) {
            path_reads = "${path_reads}" + "/"
        }
    }
}

if (params.read_type == "short_paired") {
    per_pair = true
    pattern_reads = "${params.paired_reads_pattern}${params.reads_extension}"
    fromFilePairs_input = "${path_reads}*${params.paired_reads_pattern}${params.reads_extension}"
} else{
    pattern_reads = "${params.reads_extension}"
    fromFilePairs_input = "${path_reads}*${params.reads_extension}"
}

if(input_reads.exists()){
    if ( input_reads.isDirectory()) {
        log.info "The input ${path_reads} is a folder!\n"
        input_reads.eachFileRecurse(FILES){
            if (it.name =~ ~/${pattern_reads}/){
                list_files.add(it)
            }
        }
        samples_number = list_files.size()
        log.info "${samples_number} files in ${path_reads} with pattern ${pattern_reads}"
    }
    else {
        log.info "The input ${path_reads} is a file!\n"
        pattern_reads = "${path_reads}"
    }
} else {
    exit 1, "The input ${path_reads} does not exists!\n"
}

//*************************************************
// STEP 4 - Main Workflow
//*************************************************

workflow {

    main:
        def reads
        reads = Channel.fromFilePairs(fromFilePairs_input, size: per_pair ? 2 : 1, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
      
        genome = Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
           
        align(reads,genome,aligner_list)
}


//*************************************************
// STEP 4 -  Workflow align
//*************************************************

workflow align {

    take:
        raw_reads
        genome
        aligner_list

    main:

        Channel.empty().set{logs}
        
        // ------------------------------------------------------------------------------------------------
        //                                          PREPROCESSING 
        // ------------------------------------------------------------------------------------------------

        // ------------------- QC -----------------
        if(params.fastqc){
            fastqc_raw(raw_reads, "fastqc/raw", "raw")
            logs.mix(fastqc_raw.out)
        }

        // ------------------- Standardize Score to Be Phred+33 ----------------
        seqkit_convert(raw_reads, "seqkit_score")
        seqkit_convert.out.trimmed.set{raw_reads_standardized}
        // ------------------------------------------------------------------------------------------------
        //                                          TRIMMING 
        // ------------------------------------------------------------------------------------------------        

        // ------------------- FASTP -----------------
        if (params.trimming_fastp){
            fastp(raw_reads_standardized, "fastp")
            fastp.out.trimmed.set{reads}
            logs.concat(fastp.out.report).set{logs} // save log
            if(params.fastqc){
                fastqc_fastp(reads, "fastqc/trimming_fastp", "trimmed")
                logs.concat(fastqc_fastp.out).set{logs} // save log
            }
        } else {
            reads = raw_reads_standardized
        }

        // ------------------------------------------------------------------------------------------------
        //                                          ALIGNEMENT 
        // ------------------------------------------------------------------------------------------------

        // ------------------- BBMAP -----------------
        if ("bbmap" in aligner_list ){
            bbmap_index(genome.collect(), "alignment/bbmap/indicies") // index
            bbmap(reads, genome.collect(), bbmap_index.out.collect(), "alignment/bbmap") // align
            logs.concat(bbmap.out.bbmap_summary).set{logs} // save log
            // sort
            samtools_sort_bbmap(bbmap.out.tuple_sample_bam, "alignment/bbmap")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bbmap(samtools_sort_bbmap.out.tuple_sample_sortedbam, "fastqc/bbmap", "bbmap")
                logs.concat(fastqc_ali_bbmap.out).set{logs} // save log
            }
        }

        // ------------------- BOWTIE2 -----------------
        if ("bowtie2" in aligner_list ){
            bowtie2_index(genome.collect(), "alignment/bowtie2/indicies") // index
            bowtie2(reads, genome.collect(), bowtie2_index.out.collect(), "alignment/bowtie2") // align
            logs.concat(bowtie2.out.bowtie2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie2(bowtie2.out.tuple_sample_sam)
            // sort
            samtools_sort_bowtie2(samtools_sam2bam_bowtie2.out.tuple_sample_bam, "alignment/bowtie2")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bowtie2(samtools_sort_bowtie2.out.tuple_sample_sortedbam, "fastqc/bowtie2", "bowtie2")
                logs.concat(fastqc_ali_bowtie2.out).set{logs} // save log
            }
        }

        // ------------------- BWA -----------------
        if ("bwaaln" in aligner_list || "bwamem" in aligner_list || "bwasw" in aligner_list){
            bwa_index(genome.collect(), "alignment/bwa/indicies") // index
            if ("bwaaln" in aligner_list){
                bwaaln(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwaaln") // align
                logs.concat(bwaaln.out.bwaaln_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwaaln(bwaaln.out.tuple_sample_sam)
                // sort
                samtools_sort_bwaaln(samtools_sam2bam_bwaaln.out.tuple_sample_bam, "alignment/bwa/bwaaln")
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwaaln(samtools_sort_bwaaln.out.tuple_sample_sortedbam, "fastqc/bwaaln", "bwaaln")
                    logs.concat(fastqc_ali_bwaaln.out).set{logs} // save log
                }
            }
            if ("bwamem" in aligner_list){
                bwamem(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwamem") // align
                logs.concat(bwamem.out.bwamem_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwamem(bwamem.out.tuple_sample_sam)
                // sort
                samtools_sort_bwamem(samtools_sam2bam_bwamem.out.tuple_sample_bam, "alignment/bwa/bwamem")
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwamem(samtools_sort_bwamem.out.tuple_sample_sortedbam, "fastqc/bwamem", "bwamem")
                    logs.concat(fastqc_ali_bwamem.out).set{logs} // save log
                }
            }
            if ("bwasw" in aligner_list){
                bwasw(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwasw") // align
                logs.concat(bwasw.out.bwasw_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwasw(bwasw.out.tuple_sample_sam)
                // sort
                samtools_sort_bwasw(samtools_sam2bam_bwasw.out.tuple_sample_bam, "alignment/bwa/bwasw")
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwasw(samtools_sort_bwasw.out.tuple_sample_sortedbam, "fastqc/bwasw", "bwasw")
                    logs.concat(fastqc_ali_bwasw.out).set{logs} // save log
                }
            }
        }

        // ------------------- GRAPHMAP2 -----------------
        if ("graphmap2" in aligner_list ){
            graphmap2_index(genome.collect(), "alignment/graphmap2/indicies")
            graphmap2(reads, genome.collect(), graphmap2_index.out.collect(), "alignment/graphmap2") // align
            logs.concat(graphmap2.out.graphmap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_graphmap2(graphmap2.out.tuple_sample_sam)
            // sort
            samtools_sort_graphmap2(samtools_sam2bam_graphmap2.out.tuple_sample_bam, "alignment/graphmap2")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_graphmap2(samtools_sort_graphmap2.out.tuple_sample_sortedbam, "fastqc/graphmap2", "graphmap2")
                logs.concat(fastqc_ali_graphmap2.out).set{logs} // save log
            }
        }

        // ------------------- HISAT2 -----------------
        if ("hisat2" in aligner_list){
            hisat2_index(genome.collect(),  "alignment/hisat2/indicies") // index
            hisat2(reads, hisat2_index.out.collect(), "alignment/hisat2") // align
            logs.concat(hisat2.out.hisat2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_hisat2(hisat2.out.tuple_sample_sam)
            // sort
            samtools_sort_hisat2(samtools_sam2bam_hisat2.out.tuple_sample_bam, "alignment/hisat2")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_hisat2(samtools_sort_hisat2.out.tuple_sample_sortedbam, "fastqc/hisat2", "hisat2")
                logs.concat(fastqc_ali_hisat2.out).set{logs} // save log
            }
        }

        // ------------------- minimap2 -----------------
        if ("minimap2" in aligner_list ){
            minimap2_index(genome.collect(), "alignment/minimap2/indicies") // index
            minimap2(reads, genome.collect(), minimap2_index.out.collect(), "alignment/minimap2") // align
            logs.concat(minimap2.out.minimap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_minimap2(minimap2.out.tuple_sample_sam)
            // sort
            samtools_sort_minimap2(samtools_sam2bam_minimap2.out.tuple_sample_bam, "alignment/minimap2")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_minimap2(samtools_sort_minimap2.out.tuple_sample_sortedbam, "fastqc/minimap2", "minimap2")
                logs.concat(fastqc_ali_minimap2.out).set{logs} // save log
            }
        }
        // --------------------- NGMLR --------------------
        if ("ngmlr" in aligner_list ){
            ngmlr(reads, genome.collect(), "alignment/ngmlr") // align
            logs.concat(ngmlr.out.ngmlr_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_ngmlr(ngmlr.out.tuple_sample_sam)
            // sort
            samtools_sort_ngmlr(samtools_sam2bam_ngmlr.out.tuple_sample_bam, "alignment/ngmlr")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_ngmlr(samtools_sort_ngmlr.out.tuple_sample_sortedbam, "fastqc/ngmlr", "ngmlr")
                logs.concat(fastqc_ali_ngmlr.out).set{logs} // save log
            }
        }

        // ------------------- novoalign  -----------------
        if ("novoalign" in aligner_list ){
            novoalign_index(genome.collect(), "alignment/minimap2/indicies") // index
            novoalign(reads, genome.collect(), novoalign_index.out.collect(), "alignment/novoalign") // align
            // convert sam to bam
            samtools_sam2bam_novoalign(novoalign.out.tuple_sample_sam)
            // sort
            samtools_sort_novoalign(samtools_sam2bam_novoalign.out.tuple_sample_bam, "alignment/novoalign")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_novoalign(samtools_sort_novoalign.out.tuple_sample_sortedbam, "fastqc/novoalign", "novoalign")
                logs.concat(fastqc_ali_nucmer.out).set{logs} // save log
            }
        }

        // ------------------- nucmer (mummer4) -----------------
        if ("nucmer" in aligner_list ){
            nucmer(reads, genome.collect(), "alignment/nucmer") // align
            // No summary availbe. To get one we could run show-coords see https://mummer4.github.io/tutorial/tutorial.html
            // convert sam to bam
            samtools_sam2bam_nucmer(nucmer.out.tuple_sample_sam, genome.collect())
            // sort
            samtools_sort_nucmer(samtools_sam2bam_nucmer.out.tuple_sample_bam, "alignment/nucmer")
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_nucmer(samtools_sort_nucmer.out.tuple_sample_sortedbam, "fastqc/nucmer", "nucmer")
                logs.concat(fastqc_ali_nucmer.out).set{logs} // save log
            }
        }

        // ------------------- STAR -----------------
        if ("star" in aligner_list){
            // Take read length information from only one sample for --sjdbOverhang option
            // only needed if --sjdbFileChrStartEnd or --sjdbGTFfile option is activated)
            first_file = reads.first()
            prepare_star_index_options(first_file, params.annotation, params.read_length)
            star_index(genome.collect(), prepare_star_index_options.out, "alignment/star/indicies") // index
            star(reads, star_index.out.collect(), "alignment/star") // align out is bam and sorted
            logs.concat(star.out.star_summary).set{logs} // save log
            star.out.splice_junctions.collect().set{splice_junctions} // save splice junction files
            // If  2pass mode
            if(params.star_2pass){
                star2pass(reads, star_index.out.collect(), splice_junctions, "alignment/star") // align out is bam and sorted
                logs.concat(star2pass.out.star_summary).set{logs} // save log
                star2pass.out.tuple_sample_bam.set{star_result}
            } else {
                star.out.tuple_sample_bam.set{star_result}
            }

            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_star(star_result, "fastqc/star", "star")
                logs.concat(fastqc_ali_star.out).set{logs} // save log
            }
        }

        // ---------------- subread -----------------
        if ("subread" in aligner_list ){
            subread_index(genome.collect(), "alignment/subread/indicies")
            subread(reads, genome.collect(), subread_index.out.collect(), "alignment/subread") // align
            // stat on sorted aligned reads
            if(params.fastqc){
                fastqc_ali_subread(subread.out.tuple_sample_bam, "fastqc/subread", "subread")
                logs.concat(fastqc_ali_subread.out).set{logs} // save log
            }
        }

        // ------------------- MULTIQC -----------------
        multiqc(logs.collect(),params.multiqc_config)
}


//*************************************************
def header(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    return """
    -${c_dim}--------------------------------------------------${c_reset}-
    ${c_blue}.-./`) ${c_white}.-------.    ${c_red} ______${c_reset}
    ${c_blue}\\ .-.')${c_white}|  _ _   \\  ${c_red} |    _ `''.${c_reset}     French National   
    ${c_blue}/ `-' \\${c_white}| ( ' )  |  ${c_red} | _ | ) _  \\${c_reset}    
    ${c_blue} `-'`\"`${c_white}|(_ o _) /  ${c_red} |( ''_'  ) |${c_reset}    Research Institute for    
    ${c_blue} .---. ${c_white}| (_,_).' __ ${c_red}| . (_) `. |${c_reset}
    ${c_blue} |   | ${c_white}|  |\\ \\  |  |${c_red}|(_    ._) '${c_reset}    Sustainable Development
    ${c_blue} |   | ${c_white}|  | \\ `'   /${c_red}|  (_.\\.' /${c_reset}
    ${c_blue} |   | ${c_white}|  |  \\    / ${c_red}|       .'${c_reset}
    ${c_blue} '---' ${c_white}''-'   `'-'  ${c_red}'-----'`${c_reset}
    ${c_purple} AliNe - Alignment in Nextflow - v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

// Help Message
def helpMSG() {
    log.info """
    AliNe - Alignment in Nextflow - v${workflow.manifest.version}

        Workflow: The input fastq files are standardized to Phred+33 via seqkit, trimmed if the step is activated,  
        aligned to a reference genome using one or more aligners, the ouput is converted in bam in needed, and finaly 
        sorted. If the fastqc option is activated, fastqc is run on raw and aligned reads. A multiqc report is generated containing
        fastqc output and supported aligner output.

        Usage example:
        nextflow run aline.nf --reads /path/to/reads_{1,2}.fastq.gz --genome /path/to/genome.fa --outdir alignment_results --aligner bbmap,bowtie2 --fastqc true

        --help                      prints the help section

    General Parameters
        --reads                     path to the reads file or folder
        --reads_extension           extension of the reads files (default: .fastq.gz)
        --genome                    path to the genome file
        --aligner                   aligner(s) to use among this list (comma or space separated) ${align_tools}
        --outdir                    path to the output directory (default: alignment_results)

    Type of input reads
        --read_type                 type of reads among this list ${read_type_allowed} (default: short_paired)
        --paired_reads_pattern      pattern to detect paired reads (default: {1,2})
        --stranded                  set to true if reads are stranded (default: false)

    Extra steps 
        --trimming_fastp            run fastp for trimming (default: false)
        --fastqc                    run fastqc on raw and aligned reads (default: false)
        --multiqc_config            path to the multiqc config file (default: config/multiqc_conf.yml)

    Aligner specific options
        --bbmap_options             additional options for bbmap
        --bowtie2_options           additional options for bowtie2
        --bwaaln_options            additional options for bwaaln
        --bwamem_options            additional options for bwamem
        --bwasw_options             additional options for bwasw
        --graphmap2_options         additional options for graphmap2
        --hisat2_options            additional options for hisat2
        --minimap2_options          additional options for minimap2 (default: -a (to get sam output))
        --minimap2_index_options    additional options for minimap2 index
        --ngmlr_options             additional options for ngmlr
        --novoalign_options         additional options for novoalign
        --novoalign_license         license for novoalign. You can ask for one month free trial license at http://www.novocraft.com/products/novoalign/
        --nucmer_options            additional options for nucmer
        --star_options              additional options for star
        --star_2pass                  set to true to run STAR in 2pass mode (default: false)
            --read_length               [Optional][used by STAR] length of the reads, if none provided it is automatically deduced
            --annotation                [Optional][used by STAR] path to the annotation file (gtf or gff3)
        --subread_options           additional options for subread

    """
}

def printAlignerOptions() {
    return """
    bbmap parameters
        bbmap_options              : ${params.bbmap_options}
    bowtie2 parameters
        bowtie2_options            : ${params.bowtie2_options}
    bwaaln parameters
        bwa_options                : ${params.bwaaln_options}
    bwamem parameters
        bwamem_options             : ${params.bwamem_options}
    bwasw parameters
        bwasw_options              : ${params.bwasw_options}
    graphmap2 parameters
        graphmap2_options          : ${params.graphmap2_options}
    hisat2 parameters
        hisat2_options             : ${params.hisat2_options}
    minimap2 parameters
        minimap2_options           : ${params.minimap2_options}
    ngmlr parameters
        ngmlr_options              : ${params.ngmlr_options}
    novalign parameters
        novalign_options           : ${params.novoalign_options}
        novoalign_license          : ${params.novoalign_license}
    nucmer parameters
        nucmer_options             : ${params.nucmer_options}
    star parameters
        star_options               : ${params.star_options}
        star_2pass                 : ${params.star_2pass}
    subread parameters
        subread_options            : ${params.subread_options}

    """
}


/**************         onComplete         ***************/

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