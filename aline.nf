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
params.relax = false // Avoid to automatically set option specific to ready type (e.g minimap, bwa-mem for long reads.).

// Read feature params
libtype_allowed = [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR' ]
params.library_type = "auto" 
params.skip_libray_usage = false // Avoid to use library type provided by library_type or auto
params.read_length = "" // Use by star to set the sjdbOverhang parameter
// annotation is used by different aligner (star, etc.). To avoid to duplicate processes according to the presence of the annotation file, a specific process is dedicated to create a fake file is none provided. 
// If process receive a file wich is not the fake one it includes the file in the command. To append the options of aligner we will use the annotation_file variable
// While the processes will be called sending the "annotation" channel created by the prepare_annotation process.
params.annotation = ""

// Trimming params 
params.trimming_fastp = false

// Aligner params
align_tools = [ 'bbmap', 'bowtie', 'bowtie2', 'bwaaln', 'bwamem', 'bwamem2', 'bwasw', 'graphmap2', 'hisat2', 'kallisto', 'minimap2', 'novoalign', 'nucmer', 'ngmlr', 'salmon', 'star', 'subread', 'sublong' ]
params.aligner = ''
params.bbmap_options = ''
params.bowtie_options = ''
params.bowtie2_options = ''
params.bwaaln_options = ''
params.bwamem_options = ''
params.bwamem2_options = ''
params.bwasw_options = ''
params.graphmap2_options = '' // owler option is possible
params.hisat2_options = ''
params.kallisto_options = ''
params.kallisto_index_options = '' // e.g. to use --distinguish, --make-unique, etc...
params.minimap2_options = '' 
params.minimap2_index_options = '' //  -k, -w, -H and -I
params.ngmlr_options = ''
params.novoalign_options = ''
params.novoalign_license = '' // license. You can ask for one month free trial license at http://www.novocraft.com/products/novoalign/
params.nucmer_options = ''
params.salmon_options = ''
params.salmon_index_options = ''
params.star_options = ''
params.star_index_options = ''
params.star_2pass = false
params.subread_options = '-t 0'// -t specifes the type of input sequencing data. Possible values include 0, denoting RNA-seq data, or 1, denoting genomic DNA-seq data.
params.sublong_options = '-X'// -X turn on the RNA-seq mode.

// Report params
params.fastqc = false
params.samtools_stats = false
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

// other
params.help = null
params.seqtk_sample_size = 10000 // number of reads to sample for seqtk - used to determnine the library type
bbmap_tool = "bbmap.sh"
star_tool = "STAR"
params.debug = false

//*************************************************
// STEP 1 - HELP
//*************************************************

log.info header()
if (params.help) { exit 0, helpMSG() }

//*************************************************
// STEP 1 - PARAMS CHECK
//*************************************************
def stop_pipeline = false
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

// check read library type parameter
log.info """check readlibrary type parameter: ${params.library_type} ..."""
if( ! params.library_type.contains("auto") ){
    if ( ! (params.library_type.toUpperCase() in libtype_allowed) ){
        exit 1, "Error: <${params.library_type}> library type is not accepted, please provide a library type among this list ${libtype_allowed}.\n"
    }
}

// ---- check annotation file
def annotation_file
if (params.annotation){
    if (params.annotation.startsWith('http:') || params.annotation.startsWith('https:') || params.annotation.startsWith('s3:') || params.annotation.startsWith('az:') || params.annotation.startsWith('gs:') || params.annotation.startsWith('ftp:') ) {
        annotation_file = file(params.annotation)
        annotation_file = annotation_file.getName()
    } else {
        File f = new File( "${params.annotation}" );
        if (! f.exists() ){
            log.error "Warning: Annotation file <${params.annotation}> does not exist.\n"
            stop_pipeline = true
        }
        annotation_file = f.getName()
    }
}

// ----------------------------------------------------------
// Add annotation file within the tool options if annotation provided
// Add specific options for aligner according to the read type

log.info """check alinger parameters ..."""

// --- bbmap tool ---
if ( "bbmap" in aligner_list && !params.relax ){
    if (params.read_type == "pacbio" || params.read_type == "ont"){
        bbmap_tool = "mapPacBio.sh"
        log.warn "Long reads being used, using mapPacBio.sh to align with bbmap!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter to use bbmap.sh anyway."
        // Function to check and set maxlen in params.bbmap_options when long_reads is set
        // params is supposed to be a immutable. Using params.replace method might not be supported in the future 
        if ( !params.bbmap_options.contains("maxlen") ){
            params.replace("bbmap_options", "${params.bbmap_options} maxlen=5000")
        }
    }
}

// --- bwa aln tool ---
if ( "bwaaln" in aligner_list ){
    if ( params.read_type == "pacbio" || params.read_type == "ont"){
        log.warn "Bwaaln aligner is not recommended to align long reads!\n"
    }
}

// --- bwa mem tool ---
if ( "bwamem" in aligner_list && !params.relax ){
    if (params.read_type == "pacbio"){
        if ( !params.bwamem_options.contains(" pacbio") ){
            params.replace("bwamem_options", "${params.bwamem_options} -x pacbio")
            log.warn "Pacbio reads being used, setting -x pacbio to bwamem!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
        }
    }
    if (params.read_type == "ont"){
        if ( !params.bwamem_options.contains(" ont2d") ){
            params.replace("bwamem_options", "${params.bwamem_options} -x ont2d")
            log.warn "Ont reads being used, setting -x ont2d to bwamem!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
        }
    }
}

// --- bwa mem2 tool ---
if ("bwamem2" in aligner_list && !params.relax){
    if (params.read_type == "pacbio"){
        if ( !params.bwamem2_options.contains(" pacbio") ){
            params.replace("bwamem2_options", "${params.bwamem2_options} -x pacbio")
            log.warn "Pacbio reads being used, setting -x pacbio to bwamem2!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
        }
    }
    if (params.read_type == "ont"){
        if ( !params.bwamem2_options.contains(" ont2d") ){
            params.replace("bwamem2_options", "${params.bwamem2_options} -x ont2d")
            log.warn "Ont reads being used, setting -x ont2d to bwamem2!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
        }
    }
}

// --- bwa sw tool ---
if ( "bwasw" in aligner_list ){
    if (params.read_type == "pacbio" || params.read_type == "ont"){
        log.warn "Bwasw aligner is not recommended to align long reads!\n"
    }
}

// --- graphmap2 tool ---
if ( "graphmap2" in aligner_list ){
    if (! params.read_type == "pacbio" && ! params.read_type == "ont"){
        log.warn "Graphmap2 aligner is not recommended to align short reads!\n"
    }
    if (annotation_file && !params.graphmap2_options.contains("--gtf ") ){
        params.replace("graphmap2_options", "${params.graphmap2_options} --gtf ${annotation_file}")
    }
}

// --- kallisto tool ---
if ( "kallisto" in aligner_list ){
    if ( params.read_type == "ont" || params.read_type == "pacbio"){
        log.warn "Kallisto aligner is not recommended to align long reads!\n"
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
            log.warn("""Warn: <${params.minimap2_options}> minimap2 options missing or not accepted for pacbio data. 
We set the default <map-pb> parameter. If you do not agree, please provide options among this list:
    ava-pb, splice:hq, map-hifi, map-pb (see https://github.com/lh3/minimap2).
If you wish to use parameter not intended for pacbio data use the --relax parameter to skip this warning message.\n""")
            params.replace("minimap2_options", "${params.minimap2_options} -x map-pb")
        }
    }
    if (params.read_type == "ont"){
        if ( ! params.minimap2_options.contains(" ava-ont") and ! params.minimap2_options.contains(" splice") and 
            ! params.minimap2_options.contains(" lr:hq") and ! params.minimap2_options.contains(" map-ont") ){
            log.warn("""Warn: <${params.minimap2_options}> minimap2 options missing or not accepted for ont data.
We set the default <map-ont> option. If you do not agree, please provide options among this list:
    ava-ont, splice, lr:hq, map-ont (see https://github.com/lh3/minimap2).
If you wish to use parameter not intended for pacbio data use the --relax parameter to skip this warning message.\n""")
            params.replace("minimap2_options", "${params.minimap2_options} -x map-ont")
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
                log.warn "Ont reads being used, setting -x ont to ngmlr!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
            }
        }
    }
    if (! params.read_type == "pacbio" && ! params.read_type == "ont"){
        log.warn "Ngmlr aligner is not recommended to align short reads!\n"
    }
}

// novoalign tool - load license into the container
def novoalign_lic = ""
if ("novoalign" in aligner_list ){
    if( params.novoalign_license ){
        File f = new File( "${params.novoalign_license}" );
        if (! f.exists() ){
            log.warn ": NovoAlign aligner selected but license file <${params.novoalign_license}> does not exist.\n"
            stop_pipeline = true
        }
        license_file = f.getName()
        novoalign_lic = "-v ${params.novoalign_license}:/usr/local/bin/${license_file}"
    } else {
        log.warn ": NovoAlign aligner selected but no license provided. Please provide a license to run novoalign.\n"
        stop_pipeline = true
    }

    if (!params.relax) {
        // for pacbio reads, set -g 20 and -x 0
        if (params.read_type == "pacbio" || params.read_type == "ont"){
            if (! params.novoalign_options.contains("-g ")){
                params.replace("novoalign_options", "${params.novoalign_options} -g 20")
                log.warn "Long reads being used, setting -g 20 to Novoalign!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
            }
            if (! params.novoalign_options.contains("-x ")){
                params.replace("novoalign_options", "${params.novoalign_options} -x 0")
                log.warn "Long reads being used, setting -x 0 to Novoalign!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
            }
        }
    }
}

// --- salmon tool ---
if ( "salmon" in aligner_list ){
    if ( params.read_type == "ont" || params.read_type == "pacbio"){
        log.warn "Salmon aligner is not recommended to align long reads!\n"
    }
}

// --- star tool ---
if ( "star" in aligner_list ){
    if (annotation_file && !params.star_options.contains("--sjdbGTFfile ") ){
         params.replace("star_options", "${params.star_options} --sjdbGTFfile ${annotation_file}")
    }
    if (!params.relax){
        if (params.read_type == "pacbio" || params.read_type == "ont"){
            star_tool = "STARlong"
            log.warn "Long reads being used, using STARlong to align with star!\nHowever, if you know what you are doing you can activate the AliNe --relax parameter to use star anyway."
        }
    }
}

// --- subread tool ---
if ( "subread" in aligner_list ){
    if (annotation_file && !params.subread_options.contains("-a ") ){
        params.replace("subread_options", "${params.subread_options} -a ${annotation_file}")
    }
}

if ( "sublong" in aligner_list ){
    if ( params.read_type == "short_single" || params.read_type == "short_paired"){
        log.warn "Sublong aligner is not recommended to align short reads!\n"
    }
}

if(stop_pipeline){
    exit 1, "Please fix previous issues in order to run the pipeline.\n"
}

//*************************************************
// STEP 1 - LOG INFO
//*************************************************
log.info """
General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     annotation                 : ${params.annotation}
     aligner                    : ${params.aligner}
     read_type                  : ${params.read_type}
     paired_reads_pattern       : ${params.paired_reads_pattern}
     reads_extension            : ${params.reads_extension}
     library_type               : ${params.library_type}
     skip_libray_usage          : ${params.skip_libray_usage}
     outdir                     : ${params.outdir}

Report Parameters
    fastqc                     : ${params.fastqc}
    samtools_stats             : ${params.samtools_stats}
    multiqc_config             : ${params.multiqc_config}

"""
log.info printAlignerOptions(aligner_list, annotation_file, params.star_index_options)

//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include {read_length; set_tuple_withUserReadLength} from "$baseDir/modules/bash.nf" 
include {bbmap_index; bbmap} from "$baseDir/modules/bbmap.nf"
include {bowtie_index; bowtie} from "$baseDir/modules/bowtie.nf"
include {bowtie2_index; bowtie2} from "$baseDir/modules/bowtie2.nf"
include {bwa_index; bwaaln; bwamem; bwasw} from "$baseDir/modules/bwa.nf"
include {bwamem2_index; bwamem2} from "$baseDir/modules/bwamem2.nf"
include {seqkit_convert} from "$baseDir/modules/seqkit.nf"
include {graphmap2_index; graphmap2} from "$baseDir/modules/graphmap2.nf"
include {fastp} from "$baseDir/modules/fastp.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_fastp; fastqc as fastqc_ali_bbmap; fastqc as fastqc_ali_bowtie ; fastqc as fastqc_ali_bowtie2 ; 
         fastqc as fastqc_ali_bwaaln; fastqc as fastqc_ali_bwamem; fastqc as fastqc_ali_bwamem2; fastqc as fastqc_ali_bwasw; fastqc as fastqc_ali_graphmap2 ; 
         fastqc as fastqc_ali_hisat2; fastqc as fastqc_ali_kallisto; fastqc as fastqc_ali_minimap2; fastqc as fastqc_ali_ngmlr; 
         fastqc as fastqc_ali_novoalign ; fastqc as fastqc_ali_nucmer;  fastqc as fastqc_ali_salmon; fastqc as fastqc_ali_star; fastqc as fastqc_ali_subread ; 
         fastqc as fastqc_ali_sublong } from "$baseDir/modules/fastqc.nf"
include {hisat2_index; hisat2} from "$baseDir/modules/hisat2.nf"
include {kallisto_index; kallisto} from "$baseDir/modules/kallisto.nf" 
include {minimap2_index; minimap2} from "$baseDir/modules/minimap2.nf" 
include {multiqc} from "$baseDir/modules/multiqc.nf" 
include {ngmlr} from "$baseDir/modules/ngmlr.nf" 
include {nucmer} from "$baseDir/modules/mummer4.nf" 
include {novoalign_index; novoalign} from "$baseDir/modules/novoalign.nf"
include {salmon_index; salmon_guess_lib; set_tuple_withUserLib; salmon} from "$baseDir/modules/salmon.nf" 
include {samtools_sam2bam_nucmer; samtools_sam2bam as samtools_sam2bam_bowtie; samtools_sam2bam as samtools_sam2bam_bowtie2; samtools_sam2bam as samtools_sam2bam_bwaaln; 
         samtools_sam2bam as samtools_sam2bam_bwamem; samtools_sam2bam as samtools_sam2bam_bwamem2; samtools_sam2bam as samtools_sam2bam_bwasw; samtools_sam2bam as samtools_sam2bam_graphmap2; 
         samtools_sam2bam as samtools_sam2bam_hisat2; samtools_sam2bam as samtools_sam2bam_minimap2; 
         samtools_sam2bam as samtools_sam2bam_ngmlr; samtools_sam2bam as samtools_sam2bam_novoalign; samtools_sam2bam as samtools_sam2bam_salmon } from "$baseDir/modules/samtools.nf"
include {samtools_sort as samtools_sort_bbmap; samtools_sort as samtools_sort_bowtie; samtools_sort as samtools_sort_bowtie2; samtools_sort as samtools_sort_bwaaln; 
         samtools_sort as samtools_sort_bwamem; samtools_sort as samtools_sort_bwamem2; samtools_sort as samtools_sort_bwasw; samtools_sort as samtools_sort_graphmap2; 
         samtools_sort as samtools_sort_hisat2; samtools_sort as samtools_sort_minimap2; samtools_sort as samtools_sort_ngmlr; 
         samtools_sort as samtools_sort_novoalign;  samtools_sort as samtools_sort_nucmer; samtools_sort as samtools_sort_salmon;
         samtools_sort as samtools_sort_sublong; } from "$baseDir/modules/samtools.nf"
include {samtools_stats as samtools_stats_ali_bbmap; samtools_stats as samtools_stats_ali_bowtie ; samtools_stats as samtools_stats_ali_bowtie2 ; 
         samtools_stats as samtools_stats_ali_bwaaln; samtools_stats as samtools_stats_ali_bwamem; samtools_stats as samtools_stats_ali_bwamem2; samtools_stats as samtools_stats_ali_bwasw; samtools_stats as samtools_stats_ali_graphmap2 ; 
         samtools_stats as samtools_stats_ali_hisat2; samtools_stats as samtools_stats_ali_kallisto; samtools_stats as samtools_stats_ali_minimap2; samtools_stats as samtools_stats_ali_ngmlr; 
         samtools_stats as samtools_stats_ali_novoalign ; samtools_stats as samtools_stats_ali_nucmer; samtools_stats as samtools_stats_ali_salmon; samtools_stats as samtools_stats_ali_star; 
         samtools_stats as samtools_stats_ali_subread; samtools_stats as samtools_stats_ali_sublong } from "$baseDir/modules/samtools.nf"
include {samtools_merge_bam} from "$baseDir/modules/samtools.nf"
include {seqtk_sample} from "$baseDir/modules/seqtk.nf" 
include {subread_index; subread; sublong_index; sublong} from "$baseDir/modules/subread.nf"
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

// --------- handle read input (file or folder / local or remote / paired or not) --------
def list_files = []
def pattern_reads
def fromFilePairs_input
def path_reads = params.reads 

// check if paired reads
def per_pair = false // read the reads per pair
if (params.read_type == "short_paired") {
        per_pair = true
}

// Case of remote data 
def via_URL = false
def read_list=[]
if (path_reads.startsWith('https:') || path_reads.startsWith('s3:') || path_reads.startsWith('az:') || path_reads.startsWith('gs:') || path_reads.startsWith('ftp:') ) {
    log.info "The input is a based on URLs! ${path_reads}\n"
    via_URL = true
    str_list = path_reads.tokenize(',')
    str_list.each {
        str_list2 = it.tokenize(' ')
        str_list2.each {
            read_list.add(file(it)) // use file insted of File for URL
        }
    }
    fromFilePairs_input = read_list
}
// Case of local data
else {
    File input_reads = new File(path_reads)
    if(input_reads.exists()){
        if ( input_reads.isDirectory()) {
            // in case of folder provided, add a trailing slash if missing
            path_reads = "${input_reads}" + "/"
        }
    }

    if (params.read_type == "short_paired") {
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
            fromFilePairs_input = "${path_reads}"
            if (params.read_type == "short_paired") {
                log.error "Providing a file is not authorized for (local) paired data! Please provide a folder path or change <read_type> parameter to <short_single>.\n"
            }
        }
    } else {
        exit 1, "The input ${path_reads} does not exists!\n"
    }
}
//*************************************************
// STEP 4 - Main Workflow
//*************************************************

workflow {

    main:
        // In case of URL paired data, we cannot use fromFilePairs because matching pattern impossible. We must recreate manually a structure similar 
        if (via_URL && per_pair){
            my_samples = Channel.of(read_list)
            reads = my_samples.flatten().map { it -> [it.name.split('_')[0], it] }
                                .groupTuple()
                                .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
        } else {
            reads = Channel.fromFilePairs(fromFilePairs_input, size: per_pair ? 2 : 1, checkIfExists: true)
                .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
        }

        genome = Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }

        if(annotation_file){
            annotation = Channel.fromPath(params.annotation, checkIfExists: true)
                .ifEmpty { exit 1, "Cannot find annotation matching ${annotation_file}!\n" }
        } else {
            annotation = Channel.of("$baseDir/config/aline_null.gtf") // use the fake file (not used by tools just for the processes to be called)
        }

        align(reads, genome, annotation, aligner_list)
}


//*************************************************
// STEP 4 -  Workflow align
//*************************************************

workflow align {

    take:
        raw_reads
        genome
        annotation
        aligner_list

    main:

        // Initialize channels
        Channel.empty().set{logs}
        Channel.empty().set{sorted_bam}

        // ------------------------------------------------------------------------------------------------
        //                                          PREPROCESSING 
        // ------------------------------------------------------------------------------------------------

        // ------------------- QC -----------------
        if(params.fastqc){
            fastqc_raw(raw_reads, "fastqc/raw", "raw")
            logs.concat(fastqc_raw.out).set{logs} // save log
        }

        // ------------------- Standardize Score to Be Phred+33 ----------------
        params.debug && log.info('Standardize Score to Be Phred+33')
        seqkit_convert(raw_reads, "seqkit_score")
        seqkit_convert.out.trimmed.set{raw_reads_standardized}
        // ------------------------------------------------------------------------------------------------
        //                                          TRIMMING 
        // ------------------------------------------------------------------------------------------------        

        // ------------------- FASTP -----------------
        if (params.trimming_fastp){
            params.debug && log.info('fastp trimming')
            fastp(raw_reads_standardized, "fastp")
            fastp.out.trimmed.set{raw_reads_trim}
            logs.concat(fastp.out.report).set{logs} // save log
            if(params.fastqc){
                fastqc_fastp(raw_reads_trim, "fastqc/trimming_fastp", "trimmed")
                logs.concat(fastqc_fastp.out).set{logs} // save log
            }
        } else {
            raw_reads_trim = raw_reads_standardized
        }

        // ------------------------------------------------------------------------------------------------
        //                              LIBRARY TYPE GUESSING - 
        // here we subsample and do read length and library type guessing only if needed.
        // ------------------------------------------------------------------------------------------------
        
        // In which case I need the mean read length (kallisto, salmon, star)
        // Need to list here all aligner that need the read length
        params.debug && log.info('read_length guessing')
        def guess_read_length = null
        if ( (!params.read_length) && 
           ( ("kallisto" in aligner_list && params.read_type == "short_single" && ( !params.kallisto_options.contains("-l ") || !params.kallisto_options.contains("--fragment-length  ") ) ) || 
           ( annotation.toString() != "aline_null.gtf" && !params.star_index_options.contains("--sjdbOverhang") ) || 
           ( "salmon" in aligner_list && params.read_type == "short_single" && !params.salmon.contains("--fldMean ") ) ) ){
                guess_read_length=1
        }
        params.debug && log.info('subsample guessing')
        // Should I subsample for library type guessing?
        def subsample_lg = false
        if (params.library_type.contains("auto")){
            subsample_lg = true
        }

        // Subsample if needed (can be for library guessing or read length guessing)
        if (subsample_lg || guess_read_length){
            // ------------------- subsample -----------------
            seqtk_sample(raw_reads_trim)
        }

        // ------------------- set libtype -----------------
        if (params.library_type.contains("auto")){
            // ------------------- guess libtype -------------------
            salmon_index(genome.collect())
            salmon_guess_lib(seqtk_sample.out.sampled, salmon_index.out.index, "salmon_libtype")
            salmon_guess_lib.out.tuple_id_libtype.set{tuple_id_lib}
        } else {
             set_tuple_withUserLib(raw_reads_trim)
             set_tuple_withUserLib.out.set{tuple_id_lib}
        }
        reads = raw_reads_trim.join(tuple_id_lib)

        // ------------------- set read length -----------------
        if (guess_read_length){
            read_length(seqtk_sample.out.sampled, "mean_read_length")
            read_length.out.tuple_id_readlength.set{tuple_id_readle}
        }else{
            set_tuple_withUserReadLength(raw_reads_trim)
            set_tuple_withUserReadLength.out.set{tuple_id_readle}
        }
        reads = reads.join(tuple_id_readle)
        // ------------------------------------------------------------------------------------------------
        //                                          ALIGNEMENT 
        // ------------------------------------------------------------------------------------------------

        // ------------------- BBMAP -----------------
        if ("bbmap" in aligner_list ){
            // index
            bbmap_index(genome.collect(), "alignment/bbmap/indicies")
            // align
            bbmap(reads, genome.collect(), bbmap_index.out.collect(), "alignment/bbmap")
            logs.concat(bbmap.out.bbmap_summary).set{logs} // save log
            // sort
            samtools_sort_bbmap(bbmap.out.tuple_sample_bam, "alignment/bbmap")
            samtools_sort_bbmap.out.tuple_sample_sortedbam.set{bbmap_ali} // set name
            // save aligned reads
            sorted_bam.concat(bbmap_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bbmap(bbmap_ali, "fastqc/bbmap", "bbmap")
                logs.concat(fastqc_ali_bbmap.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bbmap(bbmap_ali, genome.collect(), "samtools_stats/bbmap", "bbmap")
                logs.concat(samtools_stats_ali_bbmap.out).set{logs} // save log
            }
        }

        // ------------------- BOWTIE -----------------
        if ( "bowtie" in aligner_list ){ // &&
            bowtie_index(genome.collect(), "alignment/bowtie/indicies") // index
            bowtie(reads, genome.collect(), bowtie_index.out.collect(), "alignment/bowtie") // align
            logs.concat(bowtie.out.bowtie_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie(bowtie.out.tuple_sample_sam)
            // sort
            samtools_sort_bowtie(samtools_sam2bam_bowtie.out.tuple_sample_bam, "alignment/bowtie")
            samtools_sort_bowtie.out.tuple_sample_sortedbam.set{bowtie_ali} // set name
            // save aligned reads
            sorted_bam.concat(bowtie_ali).set{sorted_bam}
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bowtie(bowtie_ali, "fastqc/bowtie", "bowtie")
                logs.concat(fastqc_ali_bowtie.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bowtie(bowtie_ali, genome.collect(), "samtools_stats/bowtie", "bowtie")
                logs.concat(samtools_stats_ali_bowtie.out).set{logs} // save log
            }
        }

        // ------------------- BOWTIE2 -----------------
        if ( "bowtie2" in aligner_list ){ 
            // index
            bowtie2_index(genome.collect(), "alignment/bowtie2/indicies")
            // align
            bowtie2(reads, genome.collect(), bowtie2_index.out.collect(), "alignment/bowtie2")
            logs.concat(bowtie2.out.bowtie2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie2(bowtie2.out.tuple_sample_sam)
            // sort
            samtools_sort_bowtie2(samtools_sam2bam_bowtie2.out.tuple_sample_bam, "alignment/bowtie2")
            samtools_sort_bowtie2.out.tuple_sample_sortedbam.set{bowtie2_ali} // set name
            // save aligned reads
            sorted_bam.concat(bowtie2_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bowtie2(bowtie2_ali, "fastqc/bowtie2", "bowtie2")
                logs.concat(fastqc_ali_bowtie2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bowtie2(bowtie2_ali, genome.collect(), "samtools_stats/bowtie2", "bowtie2")
                logs.concat(samtools_stats_ali_bowtie2.out).set{logs} // save log
            }
        }

        // ------------------- BWA ALN/MEM/SW -----------------
        if ("bwaaln" in aligner_list || "bwamem" in aligner_list || "bwasw" in aligner_list){
            // index
            bwa_index(genome.collect(), "alignment/bwa/indicies")
            if ("bwaaln" in aligner_list){
                // align
                bwaaln(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwaaln") 
                logs.concat(bwaaln.out.bwaaln_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwaaln(bwaaln.out.tuple_sample_sam)
                // sort
                samtools_sort_bwaaln(samtools_sam2bam_bwaaln.out.tuple_sample_bam, "alignment/bwa/bwaaln")
                samtools_sort_bwaaln.out.tuple_sample_sortedbam.set{bwaaln_ali} // set name
                // save aligned reads
                sorted_bam.concat(bwaaln_ali).set{sorted_bam} 
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwaaln(bwaaln_ali, "fastqc/bwaaln", "bwaaln")
                    logs.concat(fastqc_ali_bwaaln.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwaaln(bwaaln_ali, genome.collect(), "samtools_stats/bwaaln", "bwaaln")
                    logs.concat(samtools_stats_ali_bwaaln.out).set{logs} // save log
                }
            }
            if ("bwamem" in aligner_list){
                // align
                bwamem(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwamem")
                logs.concat(bwamem.out.bwamem_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwamem(bwamem.out.tuple_sample_sam)
                // sort
                samtools_sort_bwamem(samtools_sam2bam_bwamem.out.tuple_sample_bam, "alignment/bwa/bwamem")
                samtools_sort_bwamem.out.tuple_sample_sortedbam.set{bwamem_ali} // set name
                // save aligned reads
                sorted_bam.concat(bwamem_ali).set{sorted_bam} 
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwamem(bwamem_ali, "fastqc/bwamem", "bwamem")
                    logs.concat(fastqc_ali_bwamem.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwamem(bwamem_ali, genome.collect(), "samtools_stats/bwamem", "bwamem")
                    logs.concat(samtools_stats_ali_bwamem.out).set{logs} // save log
                }
            }
            if ("bwasw" in aligner_list){
                // align
                bwasw(reads, genome.collect(), bwa_index.out.collect(), "alignment/bwa/bwasw") 
                logs.concat(bwasw.out.bwasw_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwasw(bwasw.out.tuple_sample_sam)
                // sort
                samtools_sort_bwasw(samtools_sam2bam_bwasw.out.tuple_sample_bam, "alignment/bwa/bwasw")
                samtools_sort_bwasw.out.tuple_sample_sortedbam.set{bwasw_ali} // set name
                // save aligned reads
                sorted_bam.concat(bwasw_ali).set{sorted_bam} 
                // stat on aligned reads
                if(params.fastqc){
                    fastqc_ali_bwasw(bwasw_ali, "fastqc/bwasw", "bwasw")
                    logs.concat(fastqc_ali_bwasw.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwasw(bwasw_ali, genome.collect(), "samtools_stats/bwasw", "bwasw")
                    logs.concat(samtools_stats_ali_bwasw.out).set{logs} // save log
                }
            }
        }

        // ------------------- BWA MEM2 -----------------
        if ("bwamem2" in aligner_list){
            // index
            bwamem2_index(genome.collect(), "alignment/bwamem2/indicies")
            // align
            bwamem2(reads, genome.collect(), bwamem2_index.out.collect(), "alignment/bwamem2/") 
            logs.concat(bwamem2.out.bwamem2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bwamem2(bwamem2.out.tuple_sample_sam)
            // sort
            samtools_sort_bwamem2(samtools_sam2bam_bwamem2.out.tuple_sample_bam, "alignment/bwamem2/")
            samtools_sort_bwamem2.out.tuple_sample_sortedbam.set{bwamem2_ali} // set name
            // save aligned reads
            sorted_bam.concat(bwamem2_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_bwamem2(bwamem2_ali, "fastqc/bwamem2", "bwamem2")
                logs.concat(fastqc_ali_bwamem2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bwamem2(bwamem2_ali, genome.collect(), "samtools_stats/bwamem2", "bwamem2")
                logs.concat(samtools_stats_ali_bwamem2.out).set{logs} // save log
            }
        }

        // ------------------- GRAPHMAP2 -----------------
        if ("graphmap2" in aligner_list ){
            // index
            graphmap2_index(genome.collect(), "alignment/graphmap2/indicies")
            // align
            graphmap2(reads, genome.collect(), graphmap2_index.out.collect(), annotation.collect(), "alignment/graphmap2")
            logs.concat(graphmap2.out.graphmap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_graphmap2(graphmap2.out.tuple_sample_sam)
            // sort
            samtools_sort_graphmap2(samtools_sam2bam_graphmap2.out.tuple_sample_bam, "alignment/graphmap2")
            samtools_sort_graphmap2.out.tuple_sample_sortedbam.set{graphmap2_ali} // set name
            // save aligned reads
            sorted_bam.concat(graphmap2_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_graphmap2(graphmap2_ali, "fastqc/graphmap2", "graphmap2")
                logs.concat(fastqc_ali_graphmap2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_graphmap2(graphmap2_ali, genome.collect(), "samtools_stats/graphmap2", "graphmap2")
                logs.concat(samtools_stats_ali_graphmap2.out).set{logs} // save log
            }
        }

        // ------------------- HISAT2 -----------------
        if ("hisat2" in aligner_list){
            // index
            hisat2_index(genome.collect(),  "alignment/hisat2/indicies")
            // align
            hisat2(reads, hisat2_index.out.collect(), "alignment/hisat2")
            logs.concat(hisat2.out.hisat2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_hisat2(hisat2.out.tuple_sample_sam)
            // sort
            samtools_sort_hisat2(samtools_sam2bam_hisat2.out.tuple_sample_bam, "alignment/hisat2")
            samtools_sort_hisat2.out.tuple_sample_sortedbam.set{hisat2_ali} // set name
            // save aligned reads
            sorted_bam.concat(hisat2_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_hisat2(hisat2_ali, "fastqc/hisat2", "hisat2")
                logs.concat(fastqc_ali_hisat2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_hisat2(hisat2_ali, genome.collect(), "samtools_stats/hisat2", "hisat2")
                logs.concat(samtools_stats_ali_hisat2.out).set{logs} // save log
            }
        }

        // ------------------- KALLISTO -----------------
        if ("kallisto" in aligner_list){
            // index
            kallisto_index(genome.collect(),  "alignment/kallisto/indicies")
            // align
            kallisto(reads, kallisto_index.out.collect(), "alignment/kallisto")
            logs.concat(kallisto.out.kallisto_summary).set{logs} // save log
            kallisto.out.tuple_sample_bam.set{kallisto_ali} // set name
            // save aligned reads
            sorted_bam.concat(kallisto_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_kallisto(kallisto_ali, "fastqc/kallisto", "kallisto")
                logs.concat(fastqc_ali_kallisto.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_kallisto(kallisto_ali, genome.collect(), "samtools_stats/kallisto", "kallisto")
                logs.concat(samtools_stats_ali_kallisto.out).set{logs} // save log
            }
        }
        // ------------------- minimap2 -----------------
        if ("minimap2" in aligner_list ){
            // index
            minimap2_index(genome.collect(), "alignment/minimap2/indicies")
            // align
            minimap2(reads, genome.collect(), minimap2_index.out.collect(), "alignment/minimap2")
            logs.concat(minimap2.out.minimap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_minimap2(minimap2.out.tuple_sample_sam)
            // sort
            samtools_sort_minimap2(samtools_sam2bam_minimap2.out.tuple_sample_bam, "alignment/minimap2")
            samtools_sort_minimap2.out.tuple_sample_sortedbam.set{minimap2_ali} // set name
            // save aligned reads
            sorted_bam.concat(minimap2_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_minimap2(minimap2_ali, "fastqc/minimap2", "minimap2")
                logs.concat(fastqc_ali_minimap2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_minimap2(minimap2_ali, genome.collect(), "samtools_stats/minimap2", "minimap2")
                logs.concat(samtools_stats_ali_minimap2.out).set{logs} // save log
            }
        }
        // --------------------- NGMLR --------------------
        if ("ngmlr" in aligner_list ){
            // align
            ngmlr(reads, genome.collect(), "alignment/ngmlr") 
            logs.concat(ngmlr.out.ngmlr_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_ngmlr(ngmlr.out.tuple_sample_sam)
            // sort
            samtools_sort_ngmlr(samtools_sam2bam_ngmlr.out.tuple_sample_bam, "alignment/ngmlr")
            samtools_sort_ngmlr.out.tuple_sample_sortedbam.set{ngmlr_ali} // set name
            // save aligned reads
            sorted_bam.concat(ngmlr_ali).set{sorted_bam} 
            // stat on aligned reads
            if (params.fastqc){
                fastqc_ali_ngmlr(ngmlr_ali, "fastqc/ngmlr", "ngmlr")
                logs.concat(fastqc_ali_ngmlr.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_ngmlr(ngmlr_ali, genome.collect(), "samtools_stats/ngmlr", "ngmlr")
                logs.concat(samtools_stats_ali_ngmlr.out).set{logs} // save log
            }
        }

        // ------------------- novoalign  -----------------
        if ("novoalign" in aligner_list ){
            // index
            novoalign_index(genome.collect(), "alignment/minimap2/indicies")
            // align
            novoalign(reads, genome.collect(), novoalign_index.out.collect(), novoalign_lic, "alignment/novoalign") 
            // convert sam to bam
            samtools_sam2bam_novoalign(novoalign.out.tuple_sample_sam)
            // sort
            samtools_sort_novoalign(samtools_sam2bam_novoalign.out.tuple_sample_bam, "alignment/novoalign")
            samtools_sort_novoalign.out.tuple_sample_sortedbam.set{novoalign_ali} // set name
            // save aligned reads
            sorted_bam.concat(novoalign_ali).set{sorted_bam} 
            // stat on aligned reads
            if (params.fastqc){
                fastqc_ali_novoalign(novoalign_ali, "fastqc/novoalign", "novoalign")
                logs.concat(fastqc_ali_novoalign.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_novoalign(novoalign_ali, genome.collect(), "samtools_stats/novoalign", "novoalign")
                logs.concat(samtools_stats_ali_novoalign.out).set{logs} // save log
            }
        }

        // ------------------- nucmer (mummer4) -----------------
        if ("nucmer" in aligner_list ){
            // align
            nucmer(reads, genome.collect(), "alignment/nucmer") 
            // No summary available. To get one we could run show-coords see https://mummer4.github.io/tutorial/tutorial.html
            // convert sam to bam
            samtools_sam2bam_nucmer(nucmer.out.tuple_sample_sam, genome.collect())
            // sort
            samtools_sort_nucmer(samtools_sam2bam_nucmer.out.tuple_sample_bam, "alignment/nucmer")
            samtools_sort_nucmer.out.tuple_sample_sortedbam.set{nucmer_ali} // set name
            // save aligned reads
            sorted_bam.concat(nucmer_ali).set{sorted_bam} 
            // stat on aligned reads
            if (params.fastqc){
                fastqc_ali_nucmer(nucmer_ali, "fastqc/nucmer", "nucmer")
                logs.concat(fastqc_ali_nucmer.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_nucmer(nucmer_ali, genome.collect(), "samtools_stats/nucmer", "nucmer")
                logs.concat(samtools_stats_ali_nucmer.out).set{logs} // save log
            }
        }
        
        if ("salmon" in aligner_list ){
            // index
            if (! params.library_type.contains("auto")){ // run salmon index only if library type is not provided otherwise it has been already ran
                salmon_index(genome.collect(), "alignment/salmon/indicies")
            }
            // align
            salmon(reads, salmon_index.out.collect(), "alignment/salmon")
            logs.concat(salmon.out.salmon_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_salmon(salmon.out.tuple_sample_sam)
            // sort
            samtools_sort_salmon(samtools_sam2bam_salmon.out.tuple_sample_bam, "alignment/salmon")
            samtools_sort_salmon.out.tuple_sample_sortedbam.set{salmon_ali} // set name
            // save aligned reads
            sorted_bam.concat(salmon_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_salmon(salmon_ali, "fastqc/salmon", "salmon")
                logs.concat(fastqc_ali_salmon.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_salmon(salmon_ali, genome.collect(), "samtools_stats/salmon", "salmon")
                logs.concat(samtools_stats_ali_salmon.out).set{logs} // save log
            }
        }

        // ------------------- STAR -----------------
        if ( "star" in aligner_list ){
            // Take read length information from only one sample for --sjdbOverhang option
            // only needed if --sjdbFileChrStartEnd or --sjdbGTFfile option is activated)
            first_file = reads.first()
            prepare_star_index_options(first_file, annotation.collect())
            star_index(genome.collect(), prepare_star_index_options.out, annotation, "alignment/star/indicies") // index
            star(reads, star_index.out.collect(), annotation.collect(), "alignment/star") // align out is bam and sorted
            logs.concat(star.out.star_summary).set{logs} // save log
            star.out.splice_junctions.collect().set{splice_junctions} // save splice junction files
            // If  2pass mode
            if(params.star_2pass){
                star2pass(reads, star_index.out.collect(), splice_junctions, annotation.collect(), "alignment/star") // align out is bam and sorted
                logs.concat(star2pass.out.star_summary).set{logs} // save log
                star2pass.out.tuple_sample_bam.set{star_ali} // save aligned reads
            } else {
                star.out.tuple_sample_bam.set{star_ali} // save aligned reads
            }
            // save aligned reads
            sorted_bam.concat(star_ali).set{sorted_bam} 
            // stat on aligned reads
            if(params.fastqc){
                fastqc_ali_star(star_ali, "fastqc/star", "star")
                logs.concat(fastqc_ali_star.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_star(star_ali, genome.collect(), "samtools_stats/star", "star")
                logs.concat(samtools_stats_ali_star.out).set{logs} // save log
            }
        }

        // ---------------- subread -----------------
        if ( "subread" in aligner_list ){
            // index
            subread_index(genome.collect(), "alignment/subread/indicies")
            // align
            subread(reads, genome.collect(), subread_index.out.collect(), annotation.collect(), "alignment/subread")
            subread.out.tuple_sample_bam.set{subread_ali} // set name
            // save aligned reads
            sorted_bam.concat(subread_ali).set{sorted_bam}
            // stat on sorted aligned reads
            if(params.fastqc){
                fastqc_ali_subread(subread_ali, "fastqc/subread", "subread")
                logs.concat(fastqc_ali_subread.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_subread(subread_ali, genome.collect(), "samtools_stats/subread", "subread")
                logs.concat(samtools_stats_ali_subread.out).set{logs} // save log
            }
        }

        // ---------------- sublong -----------------
        if ( "sublong" in aligner_list ){
            // index
            sublong_index(genome.collect(), "alignment/sublong/indicies")
            // align
            sublong(reads, genome.collect(), sublong_index.out.collect(), "alignment/sublong")
            sublong.out.tuple_sample_bam.set{sublong_ali_tmp} // set name
            // merge bam if paired
            if (params.read_type == "short_paired"){ 
                samtools_merge_bam(sublong.out.tuple_sample_bam)
                samtools_merge_bam.out.tuple_sample_bam.set{sublong_ali_tmp} // set name
            }
            // sort
            samtools_sort_sublong(sublong_ali_tmp, "alignment/sublong")
            samtools_sort_sublong.out.tuple_sample_sortedbam.set{sublong_ali} // set name
            // save aligned reads
            sorted_bam.concat(sublong_ali).set{sorted_bam}
            // stat on sorted aligned reads
            if(params.fastqc){
                fastqc_ali_sublong(sublong_ali, "fastqc/sublong", "sublong")
                logs.concat(fastqc_ali_sublong.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_sublong(sublong_ali, genome.collect(), "samtools_stats/sublong", "sublong")
                logs.concat(samtools_stats_ali_sublong.out).set{logs} // save log
            }
        }

        // ------------------- MULTIQC -----------------
        multiqc(logs.collect(),params.multiqc_config)

    emit:
        sorted_bam                // channel: [ val(meta), pdf ]

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
        --reads                     path to the reads file or folder. If a folder is provided, all the files with <reads_extension> in the folder will be used. You can provide remote files (commma separated list).
        --reads_extension           extension of the reads files (default: .fastq.gz)
        --genome                    path to the genome file
        --aligner                   aligner(s) to use among this list (comma or space separated) ${align_tools}
        --outdir                    path to the output directory (default: alignment_results)
        --annotation                [Optional][used by STAR, Tophat2] Absolute path to the annotation file (gtf or gff3)

    Type of input reads
        --read_type                 type of reads among this list ${read_type_allowed} (default: short_paired)
        --paired_reads_pattern      pattern to detect paired reads (default: {1,2})
        --library_type              Set the library_type of your reads (default: auto). In auto mode salmon will guess the library type for each sample.
                                    If you know the library type you can set it to one of the following: ${libtype_allowed}. See https://salmon.readthedocs.io/en/latest/library_type.html for more information.
                                    In such case the sample library type will be used for all the samples.
        --skip_libray_usage         Skip the usage of library type provided by the user or guessed by salmon. 

    Extra steps 
        --trimming_fastp            run fastp for trimming (default: false)
        --fastqc                    run fastqc on raw and aligned reads (default: false)
        --samtools_stats            run samtools stats on aligned reads (default: false)
        --multiqc_config            path to the multiqc config file (default: config/multiqc_conf.yml)

    Aligner specific options
        --bbmap_options             additional options for bbmap
        --bowtie_options            additional options for bowtie
        --bowtie2_options           additional options for bowtie2
        --bwaaln_options            additional options for bwaaln
        --bwamem_options            additional options for bwamem
        --bwamem2_options           additional options for bwamem2
        --bwasw_options             additional options for bwasw
        --graphmap2_options         additional options for graphmap2
        --hisat2_options            additional options for hisat2
        --kallisto_options          additional options for kallisto
        --kallisto_index_options    additional options for kallisto index
        --minimap2_options          additional options for minimap2 (default: -a (to get sam output))
        --minimap2_index_options    additional options for minimap2 index
        --ngmlr_options             additional options for ngmlr
        --novoalign_options         additional options for novoalign
        --novoalign_license         license for novoalign. You can ask for one month free trial license at http://www.novocraft.com/products/novoalign/
        --nucmer_options            additional options for nucmer
        --salmon_options            additional options for salmon
        --salmon_index_options      additional options for salmon index
        --star_options              additional options for star
        --star_index_options        additional options for star index
        --star_2pass                set to true to run STAR in 2pass mode (default: false)
        --read_length               [Optional][used by STAR] length of the reads, if none provided it is automatically deduced
        --subread_options           additional options for subread
        --sublong_options           additional options for sublong

    """
}

def printAlignerOptions(aligner_list, annotation_file, star_index_options) {
    def sentence = ""
    if ("bbmap" in aligner_list){ 
        sentence += """
    bbmap parameters
        bbmap_tool                 : ${bbmap_tool}
        bbmap_options              : ${params.bbmap_options}
    """} 
    if ("bowtie" in aligner_list){ 
        sentence += """       
    bowtie parameters
        bowtie_options             : ${params.bowtie_options}
    """} 
    if ("bowtie2" in aligner_list){ 
        sentence += """       
    bowtie2 parameters
        bowtie2_options            : ${params.bowtie2_options}
    """} 
    if ("bwaaln" in aligner_list){
        sentence += """
    bwaaln parameters
        bwa_options                : ${params.bwaaln_options}
    """} 
    if ("bwamem" in aligner_list){
        sentence += """
    bwamem parameters
        bwamem_options             : ${params.bwamem_options}
    """} 
    if ("bwamem2" in aligner_list){
        sentence += """
    bwamem2 parameters
        bwamem2_options            : ${params.bwamem2_options}
    """} 
    if ("bwasw" in aligner_list){
        sentence += """
    bwasw parameters
        bwasw_options              : ${params.bwasw_options}
    """} 
    if ("graphmap2" in aligner_list){
        sentence += """
    graphmap2 parameters
        graphmap2_options          : ${params.graphmap2_options}
    """} 
    if ("hisat2" in aligner_list){
        sentence += """
    hisat2 parameters
        hisat2_options             : ${params.hisat2_options}
    """}
    if ("kallisto" in aligner_list){
        def new_kallisto_sentence = "${params.kallisto_options}"
        def extra_info="" 
        if (!params.kallisto_options.contains("-l ") && !params.kallisto_options.contains("--fragment-length  ")){
            new_kallisto_sentence += " -l XXX" 
            extra_info = " xxx computed by AliNe"
        }
        if (!params.kallisto_options.contains("-s ")){
            new_kallisto_sentence += " -s YYY" 
             extra_info += " yyy computed by AliNe"
        }
        if (extra_info){
            new_kallisto_sentence += " # ${extra_info}"
        }
        sentence += """
    kallisto parameters
        kallisto_options           : ${new_kallisto_sentence}
        kallisto_index_options     : ${params.kallisto_index_options}
    """}
    if ("minimap2" in aligner_list){
        sentence += """
    minimap2 parameters
        minimap2_options           : ${params.minimap2_options}
        minimap2_index_options     : ${params.minimap2_index_options}
    """} 
    if ("ngmlr" in aligner_list){
        sentence += """
    ngmlr parameters
        ngmlr_options              : ${params.ngmlr_options}
    """} 
    if ("novoalign" in aligner_list){
        sentence += """
    novalign parameters
        novalign_options           : ${params.novoalign_options}
        novoalign_license          : ${params.novoalign_license}
    """} 
    if ("nucmer" in aligner_list){
        sentence += """
    nucmer parameters
        nucmer_options             : ${params.nucmer_options}
    """}
     if ("salmon" in aligner_list){
        def new_salmon_sentence = "${params.salmon_options}"
        def extra_info="" 
        if ( !params.salmon_options.contains("--fldMean ") ){
            new_salmon_sentence += " --fldMean XXX" 
            extra_info = " xxx computed by AliNe"
        }
        if (!params.salmon_options.contains("--fldSD ")){
            new_salmon_sentence += " --fldSD YYY" 
             extra_info += " yyy computed by AliNe"
        }
        if (extra_info){
            new_salmon_sentence += " # ${extra_info}"
        }
        sentence += """
    salmon parameters
        salmon_options             : ${new_salmon_sentence}
        salmon_index_options       : ${params.salmon_index_options}
    """}
    if ("star" in aligner_list){
        def new_index_sentence = "${star_index_options}"
        if(annotation_file){      
            if( !star_index_options.contains("--sjdbGTFfile") ){
                new_index_sentence += " --sjdbGTFfile ${annotation_file}"
            }
            if (!star_index_options.contains("--sjdbOverhang") ){
                new_index_sentence += " --sjdbOverhang XXX # xxx computed by AliNe" // to be replaced by the read length
            }
        }
        sentence += """
    star parameters
        star_index_options         : ${new_index_sentence}
        star_tool                  : ${star_tool}
        star_options               : ${params.star_options}      
        star_2pass                 : ${params.star_2pass}
    """}
    if ("subread" in aligner_list){
        sentence += """
    subread parameters
        subread_options            : ${params.subread_options}
    """}
    
    return sentence
}


/**************         onComplete         ***************/

workflow.onComplete {

    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.success) {
        log.info "\n${c_green}    AliNe pipeline complete!${c_reset}"
    } else {
        log.error "${c_red}Oops .. something went wrong${c_reset}"
    }

    log.info "    The results are available in the ‘${params.outdir}’ directory."
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

/*
Other aligner?
 - to map ont?
 - Check salmon can be used as aligner?
 - dragmap 
 - bowtie1
 
Annotation:
Star: When annotation is provided star will need read length information to index the genome. 
If no read length provided by the user and several fastq files are provided, only the first one will be used to get the read length (we perform only one index).

How to add an aliger ?
Add a module,
In aline.nf 
    Add the aligner in the aligner_list
    Import module here
    If behaves differently depending sequencing techniology add condition in the PARAMS CHECK section
    If tool need read length guessing add it in the LIBRARY TYPE GUESSING to activate the guessing if read length not provided by user.
    Add in help
    Add in printAlignerOptions
    Add process in the workflow align
    Think to convert sam to bam if necessary
    Think to sort bam output if necessary 
Add info in README.md
Add tool in multiqc config file (at least for fastqc)
*/