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
params.reference = "/path/to/reference.fa"
params.outdir = "alignment_results"
read_type_allowed = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type = ""
data_type_allowed = [ 'DNA', 'RNA' ]
params.data_type = ""
params.relax = false // Avoid to automatically set option specific to ready type (e.g minimap, bwa-mem for long reads.).

// Read feature params
strandedness_allowed = [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]
params.strandedness = ""
params.read_length = "" // Used by star to set the sjdbOverhang parameter but also by some other aligners
// annotation is used by different aligner (star, etc.). To avoid to duplicate processes according to the presence of the annotation file, a specific process is dedicated to create a fake file is none provided. 
// If process receive a file wich is not the fake one it includes the file in the command. To append the options of aligner we will use the annotation_file variable
// While the processes will be called sending the "annotation" channel created by the prepare_annotation process.
params.annotation = ""

// Trimming params 
params.trimming_fastp = false

// Aligner params
align_tools = [ 'bbmap', 'bowtie', 'bowtie2', 'bwaaln', 'bwamem', 'bwamem2', 'bwasw', 'dragmap', 'graphmap2', 'hisat2', 'kallisto', 'last', 'minimap2', 'novoalign', 'nucmer', 'ngmlr', 'salmon', 'star', 'subread', 'sublong' ]
params.aligner = ''
params.bbmap_options      = ''
params.bowtie_options     = ''
params.bowtie2_options    = ''
params.bwaaln_options     = ''
params.bwamem_options     = ''
params.bwamem2_options    = ''
params.bwasw_options      = ''
params.dragmap_options    = ''
params.graphmap2_options  = '' // owler option is possible
params.hisat2_options     = ''
params.kallisto_options   = ''
params.kallisto_index_options = '' // e.g. to use --distinguish, --make-unique, etc...
params.last_options       = ''
params.last_index_options = ''
params.minimap2_options   = '' 
params.minimap2_index_options = '' //  -k, -w, -H and -I
params.ngmlr_options      = ''
params.novoalign_options  = ''
params.novoalign_license  = '' // license. You can ask for one month free trial license at http://www.novocraft.com/products/novoalign/
params.nucmer_options     = ''
params.salmon_options     = ''
params.salmon_index_options = ''
params.star_options       = ''
params.star_index_options = ''
params.star_2pass         = false
params.subread_options    = ''
params.sublong_options    = ''

// Report params
params.fastqc = false
params.samtools_stats = false
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"
params.cram = false
params.filter_unmapped = false

// other
params.help = null
params.seqtk_sample_size = 10000 // number of reads to sample for seqtk - used to determnine the library type
params.debug = false

//*************************************************
// STEP 1 - HELP
//*************************************************

println header()
if (params.help) { exit 0, helpMSG() }

//*************************************************
// STEP 1 - PARAMS CHECK
//*************************************************
def aline_processed_params = "aline_updated_params"
def path_reads = params.reads 
def via_csv = false
if ( path_reads.endsWith('.csv') ){  
    via_csv = true
    println  "Using CSV input file: ${path_reads}"
}

println """check aligner provided: ${params.aligner} ..."""
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

// check read library type parameter
println """check strandedness parameter: ..."""
def strandednessp = params.strandedness == "null" ? null : params.strandedness // if read_type is "null" it will be set to null, otherwise it will be the value of params.strandedness
// params.strandedness can be a string "null" when coming from daisy chaining. To be really null it must be set to null here. It is why we transfer the params.strandedness into strandednessp (p for prameter)
if (! strandednessp ) {
    println """    No value provided by --strandedness"""
    if( via_csv ) {
        println """    value will be taken from the csv file."""
    }
}
else {
    if (strandednessp instanceof Boolean) {
        exit 1, "Error: --strandedness parameter needs a value among this list when used (${strandedness_allowed}), currently seen as a boolean <${strandednessp}>."
    } else {
        if ( strandednessp ){
            if ( ! (strandednessp.toUpperCase() in strandedness_allowed*.toUpperCase()) ){
                exit 1, "Error: <${strandednessp}> strandedness is not accepted, please provide a strandedness among this list ${strandedness_allowed}."
            } else{
                println """    strandedness set to : ${strandednessp}"""
                if( via_csv ) {
                    println """    This value will replace any strandedness value found in your csv!"""
                }
            }
        }
    }
}

// check data type parameter
println """check read type parameter: ${params.data_type} ..."""
if( via_csv ) {
    if( params.data_type ){
        if ( ! (params.data_type.toLowerCase() in data_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${params.data_type}> data_type not acepted, please provide a data type among this list ${data_type_allowed}."
        }
        println """    This value will replace any data_type value found in your csv!"""
    } else {
        println """    No data_type provided by --data_type parameter, value will be taken from the csv file."""
    }
} else {
    if( ! params.data_type ){
        exit 1, "Error: <data_type> parameter is empty, please provide a data type among this list ${data_type_allowed}."
    } else {
        if ( ! (params.data_type.toLowerCase() in data_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${params.data_type}> data_type not acepted, please provide a data type among this list ${data_type_allowed}."
        }
    }
}

// check read type parameter
println """check read type parameter: ${params.read_type} ..."""
def read_typep = params.read_type == "null" ? null : params.read_type // if read_type is "null" it will be set to null, otherwise it will be the value of params.read_type

if( via_csv ) {
    if( read_typep ){
        if ( ! (read_typep.toLowerCase() in read_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${read_typep}> read_type not acepted, please provide a read type among this list ${read_type_allowed}."
        }
        println """    This value will replace any read_type value found in your csv!"""
    } else {
        println """    No read_type provided by --read_type parameter, value will be taken from the csv file."""
    }
} else {
    if( ! read_typep ){
        exit 1, "Error: <read_type> parameter is empty it is allowed only when CSV provided as input. Please provide a read type among this list ${read_type_allowed}."
    } else {
        if ( ! (read_typep.toLowerCase() in read_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${read_typep}> read_type not acepted, please provide a read type among this list ${read_type_allowed}."
        }
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
            exit 1, "Error: Annotation file <${params.annotation}> does not exist.\n"
        }
        annotation_file = f.getName()
    }
}

// check license file for novoalign
if ("novoalign" in aligner_list ){
    if( params.novoalign_license ){
        File f = new File( "${params.novoalign_license}" );
        if (! f.exists() ){
            exit 1, "Error: NovoAlign aligner selected but license file <${params.novoalign_license}> does not exist.\n"
        }
        license_file = f.getName()
        novoalign_lic = "-v ${params.novoalign_license}:/usr/local/bin/${license_file}"
    } else {
        exit 1, "Error: NovoAlign aligner selected but no license provided. Please provide a license to run novoalign.\n"
    }
}

// check fastqc and cram compatibility
if (params.fastqc && params.cram) {
    println """
    WARNING: FastQC cannot process CRAM files directly. 
    Since --cram is enabled, --fastqc will be automatically disabled for alignment output files.
    FastQC will only run on raw and trimmed FASTQ files (if applicable).
    """
}

//*************************************************
// STEP 1 - LOG INFO
//*************************************************
println """
General Parameters
     reference                  : ${params.reference}
     reads                      : ${params.reads}
     annotation                 : ${params.annotation}
     aligner                    : ${params.aligner}
     data_type                  : ${params.data_type}
     read_type                  : ${read_typep}
     strandedness               : ${strandednessp}
     outdir                     : ${params.outdir}

Optional Parameters
     seqtk_sample_size          : ${params.seqtk_sample_size}
     read_length                : ${params.read_length}
     relax                      : ${params.relax}

Extra step paramesters
    trimming_fastp              : ${params.trimming_fastp}
    fastqc                      : ${params.fastqc}
    samtools_stats              : ${params.samtools_stats}
    cram                        : ${params.cram}
    filter_unmapped             : ${params.filter_unmapped}

Report Parameters
    multiqc_config              : ${params.multiqc_config}

Aligner Parameters (provided by user)
"""
println printAlignerOptions(aligner_list)
println """
Aligner parameters (updated by Aline)
    Available at                : ${params.outdir}/${aline_processed_params}
"""
//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include {read_length; check_aligner; check_aligner_params} from "$baseDir/modules/bash.nf" 
include {bbmap_index; bbmap} from "$baseDir/modules/bbmap.nf"
include {bowtie_index; bowtie} from "$baseDir/modules/bowtie.nf"
include {bowtie2_index; bowtie2} from "$baseDir/modules/bowtie2.nf"
include {bwa_index; bwaaln; bwamem; bwasw} from "$baseDir/modules/bwa.nf"
include {bwamem2_index; bwamem2} from "$baseDir/modules/bwamem2.nf"
include {dragmap_index; dragmap} from "$baseDir/modules/dragmap.nf"
include {seqkit_convert; seqkit_clean_fasta_headers} from "$baseDir/modules/seqkit.nf"
include {graphmap2_index; graphmap2} from "$baseDir/modules/graphmap2.nf"
include {fastp} from "$baseDir/modules/fastp.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_fastp} from "$baseDir/modules/fastqc.nf"
include {fastqc_ali as fastqc_ali_bbmap; fastqc_ali as fastqc_ali_bowtie ; fastqc_ali as fastqc_ali_bowtie2 ; 
    fastqc_ali as fastqc_ali_bwaaln; fastqc_ali as fastqc_ali_bwamem; fastqc_ali as fastqc_ali_bwamem2; fastqc_ali as fastqc_ali_bwasw; fastqc_ali as fastqc_ali_dragmap; fastqc_ali as fastqc_ali_graphmap2 ; 
    fastqc_ali as fastqc_ali_hisat2; fastqc_ali as fastqc_ali_kallisto; fastqc_ali as fastqc_ali_last; fastqc_ali as fastqc_ali_minimap2; fastqc_ali as fastqc_ali_ngmlr; 
    fastqc_ali as fastqc_ali_novoalign ; fastqc_ali as fastqc_ali_nucmer;  fastqc_ali as fastqc_ali_salmon; fastqc_ali as fastqc_ali_star; fastqc_ali as fastqc_ali_subread ; 
    fastqc_ali as fastqc_ali_sublong } from "$baseDir/modules/fastqc.nf"
include {hisat2_index; hisat2} from "$baseDir/modules/hisat2.nf"
include {kallisto_index; kallisto} from "$baseDir/modules/kallisto.nf"
include {last_index; last} from "$baseDir/modules/last.nf"
include {minimap2_index; minimap2} from "$baseDir/modules/minimap2.nf"
include {multiqc} from "$baseDir/modules/multiqc.nf"
include {ngmlr} from "$baseDir/modules/ngmlr.nf"
include {nucmer} from "$baseDir/modules/mummer4.nf"
include {novoalign_index; novoalign} from "$baseDir/modules/novoalign.nf"
include {fasta_uncompress} from "$baseDir/modules/pigz.nf"
include {r_rendering} from "$baseDir/modules/r.nf"
include {salmon_index; salmon_guess_lib; salmon} from "$baseDir/modules/salmon.nf" 
include {samtools_sam2bam_nucmer; samtools_sam2bam as samtools_sam2bam_bowtie; samtools_sam2bam as samtools_sam2bam_bowtie2; 
         samtools_sam2bam as samtools_sam2bam_bwaaln; samtools_sam2bam as samtools_sam2bam_bwamem; samtools_sam2bam as samtools_sam2bam_bwamem2; 
         samtools_sam2bam as samtools_sam2bam_bwasw; samtools_sam2bam as samtools_sam2bam_dragmap; samtools_sam2bam as samtools_sam2bam_graphmap2; samtools_sam2bam as samtools_sam2bam_hisat2;
         samtools_sam2bam as samtools_sam2bam_last; samtools_sam2bam as samtools_sam2bam_minimap2; 
         samtools_sam2bam as samtools_sam2bam_ngmlr; samtools_sam2bam as samtools_sam2bam_novoalign; samtools_sam2bam as samtools_sam2bam_salmon } from "$baseDir/modules/samtools.nf"
include {samtools_bam2cram as samtools_bam2cram_star; samtools_bam2cram as samtools_bam2cram_subread} from "$baseDir/modules/samtools.nf"
include {samtools_view_filter as samtools_view_filter_star; samtools_view_filter as samtools_view_filter_subread} from "$baseDir/modules/samtools.nf"
include {samtools_sort as samtools_sort_bbmap; samtools_sort as samtools_sort_bowtie; samtools_sort as samtools_sort_bowtie2; samtools_sort as samtools_sort_bwaaln; 
         samtools_sort as samtools_sort_bwamem; samtools_sort as samtools_sort_bwamem2; samtools_sort as samtools_sort_bwasw; samtools_sort as samtools_sort_dragmap; samtools_sort as samtools_sort_graphmap2; 
         samtools_sort as samtools_sort_hisat2; samtools_sort as samtools_sort_kallisto; samtools_sort as samtools_sort_last; samtools_sort as samtools_sort_minimap2; samtools_sort as samtools_sort_ngmlr; 
         samtools_sort as samtools_sort_novoalign;  samtools_sort as samtools_sort_nucmer; samtools_sort as samtools_sort_salmon;
         samtools_sort as samtools_sort_sublong; } from "$baseDir/modules/samtools.nf"
include {samtools_stats as samtools_stats_ali_bbmap; samtools_stats as samtools_stats_ali_bowtie; samtools_stats as samtools_stats_ali_bowtie2 ;
         samtools_stats as samtools_stats_ali_bwaaln; samtools_stats as samtools_stats_ali_bwamem; samtools_stats as samtools_stats_ali_bwamem2;
         samtools_stats as samtools_stats_ali_bwasw; samtools_stats as samtools_stats_ali_dragmap; samtools_stats as samtools_stats_ali_graphmap2; samtools_stats as samtools_stats_ali_hisat2; 
         samtools_stats as samtools_stats_ali_kallisto; samtools_stats as samtools_stats_ali_last; samtools_stats as samtools_stats_ali_minimap2; samtools_stats as samtools_stats_ali_ngmlr; 
         samtools_stats as samtools_stats_ali_novoalign ; samtools_stats as samtools_stats_ali_nucmer; samtools_stats as samtools_stats_ali_salmon; samtools_stats as samtools_stats_ali_star; 
         samtools_stats as samtools_stats_ali_subread; samtools_stats as samtools_stats_ali_sublong } from "$baseDir/modules/samtools.nf"
include {samtools_merge_bam_if_paired} from "$baseDir/modules/samtools.nf"
include {samtools_index as samtools_index_bbmap; samtools_index as samtools_index_bowtie; samtools_index as samtools_index_bowtie2; samtools_index as samtools_index_bwaaln;
         samtools_index as samtools_index_bwamem; samtools_index as samtools_index_bwamem2; samtools_index as samtools_index_bwasw; samtools_index as samtools_index_dragmap; samtools_index as samtools_index_graphmap2;
         samtools_index as samtools_index_hisat2; samtools_index as samtools_index_kallisto; samtools_index as samtools_index_last; samtools_index as samtools_index_minimap2;
         samtools_index as samtools_index_ngmlr; samtools_index as samtools_index_novoalign; samtools_index as samtools_index_nucmer; samtools_index as samtools_index_salmon;
         samtools_index as samtools_index_star; samtools_index as samtools_index_subread; samtools_index as samtools_index_sublong} from "$baseDir/modules/samtools.nf"
include {seqtk_sample; seqtk_sample as seqtk_sample2} from "$baseDir/modules/seqtk.nf" 
include {subread_index; subread; sublong_index; sublong} from "$baseDir/modules/subread.nf"
include {prepare_star_index_options; star_index; star; star2pass} from "$baseDir/modules/star.nf"

//*************************************************
// STEP 3 - CHECK 2 for parameters
//*************************************************

// check profile
if (
    workflow.containerEngine == 'singularity' ||
    workflow.containerEngine == 'docker'
  ) { println "executer selected: ${workflow.containerEngine}" }
else { exit 1, "No containerEngine selected: you must use a profile that activate a docker or singularity engine (-profile docker/singularity/itrop)"}

// --------- handle read input (file or folder / local or remote / paired or not) --------
def list_files = []
def pattern_reads = ""
def fromFilePairs_input
def via_URL = false
def read_list=[]
def per_pair = false // read the reads per pair
if (read_typep == "short_paired") {
        per_pair = true
}

// Case of csv file
if( ! via_csv ) {
    // Case of local data
    if( path_reads.indexOf(',') >= 0) {
        println "The input is a list!"
        via_URL = true
        // Cut into list with coma separator
        str_list = path_reads.tokenize(',')
        // loop over elements
        str_list.each {
            str_list2 = it.tokenize(' ')
            str_list2.each {
                if (  AlineUtils.is_url(it) ) {
                    println "This input is an URL: ${it}"
                }
                read_list.add(file(it)) // use file insted of File for URL
            }
        }
        // check if the list is a paired list
        if ( read_list.size() % 2 != 0 && read_typep == "short_paired") {
            exit 1, "The list has an odd number of elements which is not in line with read type <${read_typep}>."
        }
        fromFilePairs_input = read_list
    }
    else {

        File input_reads = new File(path_reads)
        if(input_reads.exists()){
            if ( input_reads.isDirectory()) {
                // in case of folder provided, add a trailing slash if missing
                path_reads = "${input_reads}" + "/"
            }
        }

        if (read_typep == "short_paired") {
            pattern_reads = "_R?[12](_\\d+)?(\\.fastq|\\.fq)(\\.gz)?\$"
            fromFilePairs_input = "${path_reads}*_{,R}{1,2}{,_*}.{fastq,fq}{,.gz}"
        } else{
            pattern_reads = "(\\.fastq|\\.fq)(\\.gz)?\$"
            fromFilePairs_input = "${path_reads}*.{fastq,fq}{,.gz}"
        }

        if(input_reads.exists()){
            if ( input_reads.isDirectory()) {
                println "The input ${path_reads} is a folder!"
                input_reads.eachFileRecurse(FILES){
                    if (it.name =~ ~/${pattern_reads}/){
                        list_files.add(it)
                        println "Found file ${it}"
                    }
                }
                samples_number = list_files.size()
                println "${samples_number} files in ${path_reads} with pattern ${pattern_reads}"
            }
            else {
                println "The input ${path_reads} is a file!\n"
                fromFilePairs_input = "${path_reads}"
                if (read_typep == "short_paired") {
                    println "Providing a file is not authorized for (local) paired data! Please provide a folder path or change <read_type> parameter to <short_single>.\n"
                }
            }
        } 
        else if ( AlineUtils.is_url(path_reads) ) {
            println "The input is a based on URLs! ${path_reads}\n"
            via_URL = true
            read_list = path_reads
        }   
        else {
            exit 1, "The input ${path_reads} does not exists!\n"
        }

    }
}
//*************************************************
// STEP 4 - Main Workflow
//*************************************************

workflow {

    main:

        // ------------------------------------------------------------------------------------------------
        //                                  HANDEL INPUT DATA (FILE AND PARAMETERS)
        // ------------------------------------------------------------------------------------------------

        // In case of URL paired data, we cannot use fromFilePairs because matching pattern impossible. We must recreate manually a structure similar 
        if (via_csv){
            File input_csv = new File(path_reads)
            if(!input_csv.exists()){ 
                error "The input ${path_reads} file does not exist!\n" 
            }
            params.debug && log.info("Using CSV input file: ${path_reads}")
            raw_reads = Channel.fromPath(path_reads)
                                .splitCsv(header: true, sep: ',')
                                .map { row ->
                                    // Check sample column
                                    if ( row.sample == null ){ 
                                            error "The input ${input_csv} file does not contain a 'sample' column!\n" 
                                    } 
                                    def sample_id    = row.sample
                                    
                                    if(row.input_1 == null && row.fastq_1 == null){ 
                                            error "The input ${input_csv} file does not contain a 'input_1' or 'fastq_1' column!\n" 
                                    }
                                    // Check input_1/fastq_1 column
                                    def fastq1;
                                    if(row.input_1) {
                                        fastq1 = file(row.input_1.trim())
                                    } else {
                                        fastq1 = file(row.fastq_1.trim())
                                    }
                                    if(! fastq1.toString().endsWith('bam')) {
                                        if (! AlineUtils.is_url(fastq1) ) {
                                                    if (! fastq1.exists() ) {
                                                        error "The input ${fastq1} file does not does not exits!\n"
                                                    }
                                        } else {
                                            log.info "This fastq input is an URL: ${fastq1}"
                                        }
                                        // Check input_2/fastq_2 column
                                        def fastq2;
                                        if(row.input_2) {
                                            fastq2 = file(row.input_2.trim())
                                        } else if (row.fastq_2) {
                                            fastq2 = file(row.fastq_2.trim())
                                        }
                                        if (fastq2){
                                            if ( ! AlineUtils.is_url(fastq2) ) {
                                                if (! fastq2.exists() ) {
                                                    error "The input ${fastq2} file does not does not exits!\n"
                                                }
                                            } else {
                                                log.info "This fastq input is an URL: ${fastq2}"
                                            }
                                        }
                                        // strandedness
                                        def libtype = "auto"
                                        if ( ! strandednessp ) { // this two parameters have priority over strand found in the csv
                                            if (row.strandedness != null) {
                                                libtype_value = row.strandedness.trim()
                                                if(libtype_value){
                                                    if ( ! ( libtype_value.toUpperCase() in strandedness_allowed*.toUpperCase()) ){
                                                        error "The input ${input_csv} file contains an invalid strandedness value: ${libtype_value}. Please provide one of the following values: ${strandedness_allowed}."
                                                    } else {
                                                        libtype = libtype_value
                                                    }
                                                } else {
                                                    log.info "The input ${input_csv} file contains an empty strandedness for sample ${sample_id}! Setting strandedness to <auto> for this sample."
                                                }
                                            } else {
                                                log.info "The input ${input_csv} file does not contain a strandedness column! Seting strandedness to <auto> to all csv samples." 
                                            }
                                        }
                                        // data type
                                        def data_type = null
                                        if ( !params.data_type ) {
                                            if (row.data_type != null) {
                                                data_type_value = row.data_type.trim().toUpperCase()
                                                if (data_type_value){
                                                    if ( ! ( data_type_value in data_type_allowed*.toUpperCase()) ){
                                                        error "The input ${input_csv} file contains an invalid data type value: ${data_type_value}. Please provide one of the following values: ${data_type_allowed}."
                                                    } else {
                                                        data_type = data_type_value
                                                    }
                                                } else {
                                                    error "The input ${input_csv} file contains an empty data_type value for sample ${sample_id}!"
                                                }
                                            } else {
                                                error """Error: The input file ${input_csv} does not contain a data_type column, and the --data_type parameter was not provided.
    Please specify the read type either by including a data_type column in the input file or by using the --data_type option."""
                                            }
                                        } else {
                                            data_type = params.data_type
                                        }

                                        // read type
                                        def read_type = null
                                        def pair = false
                                        if ( !read_typep ) {
                                            if (row.read_type != null) {
                                                read_type_value = row.read_type.trim().toLowerCase()
                                                if (read_type_value){
                                                    if ( ! ( read_type_value in read_type_allowed*.toLowerCase()) ){
                                                        error "The input ${input_csv} file contains an invalid read type value: ${read_type_value}. Please provide one of the following values: ${read_type_allowed}."
                                                    } else {
                                                        read_type = read_type_value
                                                    }
                                                } else {
                                                    error "The input ${input_csv} file contains an empty read_type value for sample ${sample_id}!"
                                                }
                                            } else {
                                                error """Error: The input file ${input_csv} does not contain a read_type column, and the --read_type parameter was not provided.
    Please specify the read type either by including a read_type column in the input file or by using the --read_type option."""
                                            }
                                        } else {
                                            read_type = read_typep
                                        }

                                        // check its is paired or not
                                        if ( fastq2 ) {
                                            if (read_type == "short_paired") {
                                                pair = true
                                                file_id2 = AlineUtils.cleanPrefix(fastq2)
                                            } else {
                                                log.info "The input ${input_csv} file contains a second fastq file for sample ${sample_id} but the read_type is set to <${read_type}>! R2 will not be taken into account! paired set to false."
                                            }
                                        } else {
                                            if (read_type == "short_paired") {
                                                error "The input ${input_csv} file does not contain a second fastq file for sample ${sample_id} but the read_type is set to <short_paired>!"
                                            }
                                        }
                                        // Create a tuple with metadata and reads
                                        def file_id1 = AlineUtils.cleanPrefix(fastq1)
                                        def uid = AlineUtils.get_file_uid(fastq1)
                                        def file_id = pair ? [file_id1, file_id2] : [file_id1]
                                        def meta = [ uid: uid, file_id: file_id, sample_id: sample_id, strandedness: libtype, read_type: read_type, data_type: data_type, paired: pair ]
                                        def reads = pair ? [fastq1, fastq2] : fastq1
                                        // Return only if the fastq file(s) extension are valid
                                        if ( AlineUtils.is_fastq( fastq1.toString() ) and ( ! fastq2 || AlineUtils.is_fastq( fastq2.toString() ) ) ){
                                            return tuple(meta, reads)
                                        }
                                    } else {
                                         log.info "The input file for ${sample_id} sample is a BAM file! Discarded from AliNe alignment..."
                                    }
                                }
        }
        else {
            if (via_URL && per_pair){
            my_samples = Channel.of(read_list)
            raw_reads = my_samples.flatten().map { it -> 
                                                def uid = AlineUtils.get_file_uid(it)
                                                [uid, it] }
                                .groupTuple()
                                .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
            } else {
                log.info "Equivalent fromFilePairs regex: ${fromFilePairs_input} \n"
                raw_reads = Channel.fromFilePairs(fromFilePairs_input, size: per_pair ? 2 : 1, checkIfExists: true)
                    .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
            }
            raw_reads = raw_reads.map { uid, files -> 
                def file_id = files.size() > 1 ? files.collect { AlineUtils.cleanPrefix(it) } : [AlineUtils.cleanPrefix(files[0])]
                [ [ uid: uid, file_id: file_id ], files ] 
            }
        }

        // ------------------------ deal with reference file ------------------------
        Channel.fromPath(params.reference, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find reference matching ${params.reference}!\n" }
            .set{reference_raw}
        // uncompress it if needed because some aligner need the reference to be uncompressed (e.g. histat2)
        fasta_uncompress(reference_raw)
        // clean fasta headers to avoid issues with CRAM and other tools
        seqkit_clean_fasta_headers(fasta_uncompress.out.genomeFa, "reference")
        seqkit_clean_fasta_headers.out.clean_fasta.set{reference} // set reference to the cleaned fasta
       
        // ------------------------ deal with annotation file ------------------------
        def annotation_provided = null
        if(annotation_file){
            annotation = Channel.fromPath(params.annotation, checkIfExists: true)
                .ifEmpty { exit 1, "Cannot find annotation matching ${annotation_file}!\n" }
            annotation_provided = true
        } else {
            annotation = Channel.of("$baseDir/config/aline_null.gtf") // use the fake file (not used by tools just for the processes to be called)
        }
        raw_reads = raw_reads.map { meta, files -> [ meta + [ annotation: annotation_provided ], files ] }
        params.debug && log.info("Set annotation")
        params.debug && raw_reads.view()

        // ------------------------ data_type ------------------------ 
        // // By default priority to data_type from --data_type over csv value (the same for all csv params e.g. strandedness, read_type)
        if (params.data_type) {
            params.debug && log.info("Set data_type value from parameter: ${params.data_type}")
            raw_reads = raw_reads.map { meta, files -> [ meta + [ data_type: params.data_type ], files ] }
            params.debug && raw_reads.view()
        }
        
        // ------------------------ read_type ------------------------ 
        // // By default priority to read_type from --read_type over csv value (the same for all csv params e.g. strandedness, data_type)
        if (read_typep) {
            params.debug && log.info("Set read_type and paired meta value from parameter: ${read_typep}")
            raw_reads = raw_reads.map { meta, files -> [ meta + [ read_type: read_typep ], files ] }
            pair = (read_typep == "short_paired") ? true : false
            raw_reads = raw_reads.map { meta, files -> [ meta + [ paired: pair ], files ] }
            params.debug && raw_reads.view()
        }
        
        // -------------- Warning aligner read_type usage ---------------
        log.info """Check aligner ..."""
        check_aligner( raw_reads, aligner_list )

        // -------------- Other ---------------
        // Initialize channels
        Channel.empty().set{logs}
        // extra params
        salmon_index_done = false // to avoid multiple calls to salmon_index
        // add AliNe suffix for output files - can be used to extract uid in the output file name 
        raw_reads = raw_reads.map { meta, files -> [ meta + [ suffix: "_AliNe" ], files ] }

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
        //                              - READ LENGTH GUESSING - 
        // ------------------------------------------------------------------------------------------------
        Channel.empty().set{subsampled}
        params.debug && log.info('read_length guessing')
       
        if ( params.read_length ) {
            // add read length in meta from params
            params.debug && log.info('read_length provided by parameter')
            raw_reads_trim.map { meta, files -> [ meta + [ read_length: params.read_length ], files ] }
                          .set{raw_reads_trim_length}
            raw_reads_trim.map { meta, files -> [ meta + [ read_length: params.read_length ], files ] }
                          .set{raw_reads_trim_length}
        } // STAR case (touch all type of sample) - need to set sjdbOverhang when we do have an annotation
        else if ( "star" in aligner_list && params.annotation && !params.star_index_options.contains("--sjdbOverhang") ) {
            params.debug && log.info('read_length star case')
            subsampled = seqtk_sample(raw_reads_trim)
            read_length(subsampled, "mean_read_length")
                
            // add read length in meta
            raw_reads_trim.map { meta, fastq -> tuple(meta.uid, meta, fastq) }
                        .join(read_length.out.tuple_id_readlength.map { meta, length -> tuple(meta.uid, length) })
                        .map { uid, meta, fastq, length ->
                            def updated_meta = meta + [read_length: length, subsampled: true]
                            tuple(updated_meta, fastq)
                        }
                        .set { raw_reads_trim_length }

        } // kallisto, salmon case that touch only short_single samples
        else if ( ("kallisto" in aligner_list && ( !params.kallisto_options.contains("-l ") || !params.kallisto_options.contains("--fragment-length  ") ) ) || 
                   ( "salmon" in aligner_list && !params.salmon_options.contains("--fldMean ") ) 
                ) {
                params.debug && log.info('kallisto or salmon case')
                // For kallisto and salmon, we need to guess the read length only if short_single
                raw_reads_trim_short_single = raw_reads_trim.filter { meta, reads -> !meta.paired }
                raw_reads_trim_others = raw_reads_trim.filter { meta, reads -> meta.paired }

                params.debug && log.info('subsample short_single reads for read length guessing')
                subsampled = seqtk_sample(raw_reads_trim_short_single)
                read_length(subsampled, "mean_read_length")

                // add read length in meta
                raw_reads_trim_short_single.map { meta, fastq -> tuple(meta.uid, meta, fastq) }
                                           .join(read_length.out.tuple_id_readlength.map { meta, length -> tuple(meta.uid, length) })
                                           .map { uid, meta, fastq, length ->
                                                def updated_meta = meta + [read_length: length, subsampled: true]
                                                tuple(updated_meta, fastq)
                                                }
                                            .set{raw_reads_trim_short_single_length}

                raw_reads_trim_length = raw_reads_trim_short_single_length.concat(raw_reads_trim_others)
        }
        // else we do not add read length in meta
        else {
            params.debug && log.info('read_length not needed case')
            raw_reads_trim_length = raw_reads_trim
        }

        params.debug && log.info("raw_reads_trim_length channel output: \n")
        params.debug && raw_reads_trim_length.view()

        // ------------------------------------------------------------------------------------------------
        //                              - LIBRARY TYPE GUESSING - 
        // ------------------------------------------------------------------------------------------------
        params.debug && log.info('library type guessing')

        // If strandednessp is empty, we do not guess strandedness
        if ( ! strandednessp ) {
            params.debug && log.info('No value provided for strandedness parameter => strandedness set to null')
            // add set strandedness to null only if not already set (e.g. specific case of salmon short_single reads)
            raw_reads_trim_length.map { meta, files ->
                                        def updated_meta = meta.containsKey('strandedness') ? meta : meta + [ strandedness: null ]
                                        [ updated_meta, files ]
                                    }
                                 .set{raw_reads_trim_length_strandedness}
        } else {
              if ( strandedness_allowed.contains(strandednessp)){
                params.debug && log.info("Parameter strandedness in use => strandedness set to ${strandednessp}.")
                raw_reads_trim_length.map { meta, files -> [ meta + [ strandedness: strandednessp ], files ] }
                                     .set{raw_reads_trim_length_strandedness}
              } else {
                raw_reads_trim_length_strandedness = raw_reads_trim_length
              }
        }
        
        // Separate samples to guess strandedness and not guess strandedness
        sample_to_guess = raw_reads_trim_length_strandedness.filter { meta, reads -> meta.strandedness == 'auto' }
        sample_to_notguess = raw_reads_trim_length_strandedness.filter { meta, reads -> meta.strandedness != 'auto' }

        // catch what is already subsampled
        sample_to_guess_already_subsampled = sample_to_guess.filter { meta, reads -> meta.subsampled }
        subsample_sample_to_guess_already_subsampled = sample_to_guess_already_subsampled.map { meta, reads -> tuple(meta.uid, meta, reads) }
                                                                                        .join( subsampled.map { meta2, subreads -> tuple(meta2.uid, meta2, subreads) } )
                                                                                        .map { uid, meta, reads, meta2, subreads ->
                                                                                            tuple(meta2, subreads)
                                                                                        }
        // subsample whath has to be subsampled
        sample_to_guess_to_subsampled = sample_to_guess.filter { meta, reads ->  !meta.subsampled }
        subsampled_sample_to_guess_to_subsampled = seqtk_sample2(sample_to_guess_to_subsampled)

        // Merge all subsampled reads for strandedness guessing
        all_subsampled_read_guessing = subsample_sample_to_guess_already_subsampled.concat(subsampled_sample_to_guess_to_subsampled)

        // ------------------- guess strandedness -------------------
        salmon_index_ch = salmon_index(reference.collect(), "alignment/salmon/indicies" )
        salmon_index_done = true // set to true to avoid multiple calls to salmon_index
        salmon_guess_lib(all_subsampled_read_guessing, salmon_index_ch, "salmon_strandedness")

        // add strandedness in meta
        sample_to_guess.map { meta, fastq -> tuple(meta.uid, meta, fastq) }
        .join(salmon_guess_lib.out.tuple_id_libtype.map { meta, libtype -> tuple(meta.uid, libtype) })
        .map { uid, meta, fastq, libtype ->
            def updated_meta = meta + [strandedness: libtype, subsampled: true]
            tuple(updated_meta, fastq)
        }
        .set { sample_to_guess_done }

        // ------------------- merge sample_to_guess and sample_to_notguess -------------------
        raw_reads_trim_length_lib = sample_to_notguess.concat(sample_to_guess_done)
        params.debug && raw_reads_trim_length_lib.view()

        // ------------------------------------------------------------------------------------------------
        //                                          ADAPT ALIGNER PARAMETERS 
        // ------------------------------------------------------------------------------------------------
        log.info """Adapt aligner parameters ..."""
        check_aligner_params( raw_reads_trim_length_lib, aligner_list, annotation.collect(), aline_processed_params )
        reads = check_aligner_params.out[0]
        params.debug && reads.view()

        // ------------------------------------------------------------------------------------------------
        //                                          ALIGNEMENT 
        // ------------------------------------------------------------------------------------------------
        // Initialize sorted_ali channel
        Channel.empty().set{sorted_ali}
        
        params.debug && log.info('library type alignment')
        // ------------------- BBMAP -----------------
        if ("bbmap" in aligner_list ){
            // index
            bbmap_index(reference.collect(), "alignment/bbmap/indicies")
            // align
            bbmap(reads, reference.collect(), bbmap_index.out.collect(), "alignment/bbmap")
            logs.concat(bbmap.out.bbmap_summary).set{logs} // save log
            // sort and convert to cram
            samtools_sort_bbmap(bbmap.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_bbmap(samtools_sort_bbmap.out, "alignment/bbmap")
            samtools_index_bbmap.out.tuple_sample_ali.set{bbmap_ali} // set name
            // save aligned reads
            sorted_ali.concat(bbmap_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_bbmap(bbmap_ali, "fastqc/bbmap", "bbmap")
                logs.concat(fastqc_ali_bbmap.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bbmap(bbmap_ali, reference.collect(), "samtools_stats/bbmap", "bbmap")
                logs.concat(samtools_stats_ali_bbmap.out).set{logs} // save log
            }
        }

        // ------------------- BOWTIE -----------------
        if ( "bowtie" in aligner_list ){ // &&
            bowtie_index(reference.collect(), "alignment/bowtie/indicies") // index
            bowtie(reads, reference.collect(), bowtie_index.out.collect(), "alignment/bowtie") // align
            logs.concat(bowtie.out.bowtie_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie(bowtie.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_bowtie(samtools_sam2bam_bowtie.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_bowtie(samtools_sort_bowtie.out, "alignment/bowtie")
            samtools_index_bowtie.out.tuple_sample_ali.set{bowtie_ali} // set name
            // save aligned reads
            sorted_ali.concat(bowtie_ali).set{sorted_ali}
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_bowtie(bowtie_ali, "fastqc/bowtie", "bowtie")
                logs.concat(fastqc_ali_bowtie.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bowtie(bowtie_ali, reference.collect(), "samtools_stats/bowtie", "bowtie")
                logs.concat(samtools_stats_ali_bowtie.out).set{logs} // save log
            }
        }

        // ------------------- BOWTIE2 -----------------
        if ( "bowtie2" in aligner_list ){ 
            // index
            bowtie2_index(reference.collect(), "alignment/bowtie2/indicies")
            // align
            bowtie2(reads, reference.collect(), bowtie2_index.out.collect(), "alignment/bowtie2")
            logs.concat(bowtie2.out.bowtie2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bowtie2(bowtie2.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_bowtie2(samtools_sam2bam_bowtie2.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_bowtie2(samtools_sort_bowtie2.out, "alignment/bowtie2")
            samtools_index_bowtie2.out.tuple_sample_ali.set{bowtie2_ali} // set name
            // save aligned reads
            sorted_ali.concat(bowtie2_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_bowtie2(bowtie2_ali, "fastqc/bowtie2", "bowtie2")
                logs.concat(fastqc_ali_bowtie2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bowtie2(bowtie2_ali, reference.collect(), "samtools_stats/bowtie2", "bowtie2")
                logs.concat(samtools_stats_ali_bowtie2.out).set{logs} // save log
            }
        }

        // ------------------- BWA ALN/MEM/SW -----------------
        if ("bwaaln" in aligner_list || "bwamem" in aligner_list || "bwasw" in aligner_list){
            // index
            bwa_index(reference.collect(), "alignment/bwa/indicies")
            if ("bwaaln" in aligner_list){
                // align
                bwaaln(reads, reference.collect(), bwa_index.out.collect(), "alignment/bwa/bwaaln") 
                logs.concat(bwaaln.out.bwaaln_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwaaln(bwaaln.out.tuple_sample_sam)
                // sort and convert to cram
                samtools_sort_bwaaln(samtools_sam2bam_bwaaln.out.tuple_sample_bam, reference.collect())
                // index
                samtools_index_bwaaln(samtools_sort_bwaaln.out, "alignment/bwa/bwaaln")
                samtools_index_bwaaln.out.tuple_sample_ali.set{bwaaln_ali} // set name
                // save aligned reads
                sorted_ali.concat(bwaaln_ali).set{sorted_ali} 
                // stat on aligned reads
                if(params.fastqc && !params.cram){
                    fastqc_ali_bwaaln(bwaaln_ali, "fastqc/bwaaln", "bwaaln")
                    logs.concat(fastqc_ali_bwaaln.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwaaln(bwaaln_ali, reference.collect(), "samtools_stats/bwaaln", "bwaaln")
                    logs.concat(samtools_stats_ali_bwaaln.out).set{logs} // save log
                }
            }
            if ("bwamem" in aligner_list){
                // align
                bwamem(reads, reference.collect(), bwa_index.out.collect(), "alignment/bwa/bwamem")
                logs.concat(bwamem.out.bwamem_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwamem(bwamem.out.tuple_sample_sam)
                // sort and convert to cram
                samtools_sort_bwamem(samtools_sam2bam_bwamem.out.tuple_sample_bam, reference.collect())
                // index
                samtools_index_bwamem(samtools_sort_bwamem.out, "alignment/bwa/bwamem")
                samtools_index_bwamem.out.tuple_sample_ali.set{bwamem_ali} // set name
                // save aligned reads
                sorted_ali.concat(bwamem_ali).set{sorted_ali} 
                // stat on aligned reads
                if(params.fastqc && !params.cram){
                    fastqc_ali_bwamem(bwamem_ali, "fastqc/bwamem", "bwamem")
                    logs.concat(fastqc_ali_bwamem.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwamem(bwamem_ali, reference.collect(), "samtools_stats/bwamem", "bwamem")
                    logs.concat(samtools_stats_ali_bwamem.out).set{logs} // save log
                }
            }
            if ("bwasw" in aligner_list){
                // align
                bwasw(reads, reference.collect(), bwa_index.out.collect(), "alignment/bwa/bwasw") 
                logs.concat(bwasw.out.bwasw_summary).set{logs} // save log
                // convert sam to bam
                samtools_sam2bam_bwasw(bwasw.out.tuple_sample_sam)
                // sort and convert to cram
                samtools_sort_bwasw(samtools_sam2bam_bwasw.out.tuple_sample_bam, reference.collect())
                // index
                samtools_index_bwasw(samtools_sort_bwasw.out, "alignment/bwa/bwasw")
                samtools_index_bwasw.out.tuple_sample_ali.set{bwasw_ali} // set name
                // save aligned reads
                sorted_ali.concat(bwasw_ali).set{sorted_ali} 
                // stat on aligned reads
                if(params.fastqc && !params.cram){
                    fastqc_ali_bwasw(bwasw_ali, "fastqc/bwasw", "bwasw")
                    logs.concat(fastqc_ali_bwasw.out).set{logs} // save log
                }
                if(params.samtools_stats){
                    samtools_stats_ali_bwasw(bwasw_ali, reference.collect(), "samtools_stats/bwasw", "bwasw")
                    logs.concat(samtools_stats_ali_bwasw.out).set{logs} // save log
                }
            }
        }

        // ------------------- BWA MEM2 -----------------
        if ("bwamem2" in aligner_list){
            // index
            bwamem2_index(reference.collect(), "alignment/bwamem2/indicies")
            // align
            bwamem2(reads, reference.collect(), bwamem2_index.out.collect(), "alignment/bwamem2/") 
            logs.concat(bwamem2.out.bwamem2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_bwamem2(bwamem2.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_bwamem2(samtools_sam2bam_bwamem2.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_bwamem2(samtools_sort_bwamem2.out, "alignment/bwamem2/")
            samtools_index_bwamem2.out.tuple_sample_ali.set{bwamem2_ali} // set name
            // save aligned reads
            sorted_ali.concat(bwamem2_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_bwamem2(bwamem2_ali, "fastqc/bwamem2", "bwamem2")
                logs.concat(fastqc_ali_bwamem2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_bwamem2(bwamem2_ali, reference.collect(), "samtools_stats/bwamem2", "bwamem2")
                logs.concat(samtools_stats_ali_bwamem2.out).set{logs} // save log
            }
        }

        // ------------------- DRAGMAP -----------------
        if ("dragmap" in aligner_list){
            // index
            dragmap_index(reference.collect(), "alignment/dragmap/indicies")
            // align
            dragmap(reads, reference.collect(), dragmap_index.out.collect(), "alignment/dragmap")
            logs.concat(dragmap.out.dragmap_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_dragmap(dragmap.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_dragmap(samtools_sam2bam_dragmap.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_dragmap(samtools_sort_dragmap.out, "alignment/dragmap")
            samtools_index_dragmap.out.tuple_sample_ali.set{dragmap_ali} // set name
            // save aligned reads
            sorted_ali.concat(dragmap_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_dragmap(dragmap_ali, "fastqc/dragmap", "dragmap")
                logs.concat(fastqc_ali_dragmap.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_dragmap(dragmap_ali, reference.collect(), "samtools_stats/dragmap", "dragmap")
                logs.concat(samtools_stats_ali_dragmap.out).set{logs} // save log
            }
        }

        // ------------------- GRAPHMAP2 -----------------
        if ("graphmap2" in aligner_list ){
            // index
            graphmap2_index(reference.collect(), "alignment/graphmap2/indicies")
            // align
            graphmap2(reads, reference.collect(), graphmap2_index.out.collect(), annotation.collect(), "alignment/graphmap2")
            logs.concat(graphmap2.out.graphmap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_graphmap2(graphmap2.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_graphmap2(samtools_sam2bam_graphmap2.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_graphmap2(samtools_sort_graphmap2.out, "alignment/graphmap2")
            samtools_index_graphmap2.out.tuple_sample_ali.set{graphmap2_ali} // set name
            // save aligned reads
            sorted_ali.concat(graphmap2_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_graphmap2(graphmap2_ali, "fastqc/graphmap2", "graphmap2")
                logs.concat(fastqc_ali_graphmap2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_graphmap2(graphmap2_ali, reference.collect(), "samtools_stats/graphmap2", "graphmap2")
                logs.concat(samtools_stats_ali_graphmap2.out).set{logs} // save log
            }
        }

        // ------------------- HISAT2 -----------------
        if ("hisat2" in aligner_list){
            // index
            hisat2_index(reference.collect(),  "alignment/hisat2/indicies")
            // align
            hisat2(reads, hisat2_index.out.collect(), "alignment/hisat2")
            logs.concat(hisat2.out.hisat2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_hisat2(hisat2.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_hisat2(samtools_sam2bam_hisat2.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_hisat2(samtools_sort_hisat2.out, "alignment/hisat2")
            samtools_index_hisat2.out.tuple_sample_ali.set{hisat2_ali} // set name
            // save aligned reads
            sorted_ali.concat(hisat2_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_hisat2(hisat2_ali, "fastqc/hisat2", "hisat2")
                logs.concat(fastqc_ali_hisat2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_hisat2(hisat2_ali, reference.collect(), "samtools_stats/hisat2", "hisat2")
                logs.concat(samtools_stats_ali_hisat2.out).set{logs} // save log
            }
        }

        // ------------------- KALLISTO ----------------- output is bam but not sorted and bai missing
        if ("kallisto" in aligner_list){
            // index
            kallisto_index(reference.collect(),  "alignment/kallisto/indicies")
            // align
            kallisto(reads, kallisto_index.out.collect(), "alignment/kallisto")
            logs.concat(kallisto.out.kallisto_summary).set{logs} // save log
            // sort and convert to cram
            samtools_sort_kallisto(kallisto.out.tuple_sample_bam, reference.collect())
            // index bam/cram files
            samtools_index_kallisto(samtools_sort_kallisto.out, "alignment/kallisto")
            samtools_index_kallisto.out.tuple_sample_ali.set{kallisto_ali} // set name
            // save aligned reads
            sorted_ali.concat(kallisto_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_kallisto(kallisto_ali, "fastqc/kallisto", "kallisto")
                logs.concat(fastqc_ali_kallisto.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_kallisto(kallisto_ali, reference.collect(), "samtools_stats/kallisto", "kallisto")
                logs.concat(samtools_stats_ali_kallisto.out).set{logs} // save log
            }
        }

        // --------------------- LAST --------------------
        if ("last" in aligner_list ){
            // index
            last_index(reference.collect(), "alignment/last/indicies")
            // align
            last(reads, reference.collect(), last_index.out.collect(), "alignment/last") 
            logs.concat(last.out.last_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_last(last.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_last(samtools_sam2bam_last.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_last(samtools_sort_last.out, "alignment/last")
            samtools_index_last.out.tuple_sample_ali.set{last_ali} // set name
            // save aligned reads
            sorted_ali.concat(last_ali).set{sorted_ali} 
            // stat on aligned reads
            if (params.fastqc && !params.cram){
                fastqc_ali_last(last_ali, "fastqc/last", "last")
                logs.concat(fastqc_ali_last.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_last(last_ali, reference.collect(), "samtools_stats/last", "last")
                logs.concat(samtools_stats_ali_last.out).set{logs} // save log
            }
        }

        // ------------------- minimap2 -----------------
        if ("minimap2" in aligner_list ){
            // index
            minimap2_index(reference.collect(), "alignment/minimap2/indicies")
            // align
            minimap2(reads, reference.collect(), minimap2_index.out.collect(), "alignment/minimap2")
            logs.concat(minimap2.out.minimap2_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_minimap2(minimap2.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_minimap2(samtools_sam2bam_minimap2.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_minimap2(samtools_sort_minimap2.out, "alignment/minimap2")
            samtools_index_minimap2.out.tuple_sample_ali.set{minimap2_ali} // set name
            // save aligned reads
            sorted_ali.concat(minimap2_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_minimap2(minimap2_ali, "fastqc/minimap2", "minimap2")
                logs.concat(fastqc_ali_minimap2.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_minimap2(minimap2_ali, reference.collect(), "samtools_stats/minimap2", "minimap2")
                logs.concat(samtools_stats_ali_minimap2.out).set{logs} // save log
            }
        }

        // --------------------- NGMLR --------------------
        if ("ngmlr" in aligner_list ){
            // align
            ngmlr(reads, reference.collect(), "alignment/ngmlr") 
            logs.concat(ngmlr.out.ngmlr_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_ngmlr(ngmlr.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_ngmlr(samtools_sam2bam_ngmlr.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_ngmlr(samtools_sort_ngmlr.out, "alignment/ngmlr")
            samtools_index_ngmlr.out.tuple_sample_ali.set{ngmlr_ali} // set name
            // save aligned reads
            sorted_ali.concat(ngmlr_ali).set{sorted_ali} 
            // stat on aligned reads
            if (params.fastqc && !params.cram){
                fastqc_ali_ngmlr(ngmlr_ali, "fastqc/ngmlr", "ngmlr")
                logs.concat(fastqc_ali_ngmlr.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_ngmlr(ngmlr_ali, reference.collect(), "samtools_stats/ngmlr", "ngmlr")
                logs.concat(samtools_stats_ali_ngmlr.out).set{logs} // save log
            }
        }

        // ------------------- novoalign  -----------------
        if ("novoalign" in aligner_list ){
            // index
            novoalign_index(reference.collect(), "alignment/minimap2/indicies")
            // align
            novoalign(reads, reference.collect(), novoalign_index.out.collect(), novoalign_lic, "alignment/novoalign") 
            // convert sam to bam
            samtools_sam2bam_novoalign(novoalign.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_novoalign(samtools_sam2bam_novoalign.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_novoalign(samtools_sort_novoalign.out, "alignment/novoalign")
            samtools_index_novoalign.out.tuple_sample_ali.set{novoalign_ali} // set name
            // save aligned reads
            sorted_ali.concat(novoalign_ali).set{sorted_ali} 
            // stat on aligned reads
            if (params.fastqc && !params.cram){
                fastqc_ali_novoalign(novoalign_ali, "fastqc/novoalign", "novoalign")
                logs.concat(fastqc_ali_novoalign.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_novoalign(novoalign_ali, reference.collect(), "samtools_stats/novoalign", "novoalign")
                logs.concat(samtools_stats_ali_novoalign.out).set{logs} // save log
            }
        }

        // ------------------- nucmer (mummer4) -----------------
        if ("nucmer" in aligner_list ){
            // align
            nucmer(reads, reference.collect(), "alignment/nucmer") 
            // No summary available. To get one we could run show-coords see https://mummer4.github.io/tutorial/tutorial.html
            // convert sam to bam
            samtools_sam2bam_nucmer(nucmer.out.tuple_sample_sam, reference.collect())
            // sort and convert to cram
            samtools_sort_nucmer(samtools_sam2bam_nucmer.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_nucmer(samtools_sort_nucmer.out, "alignment/nucmer")
            samtools_index_nucmer.out.tuple_sample_ali.set{nucmer_ali} // set name
            // save aligned reads
            sorted_ali.concat(nucmer_ali).set{sorted_ali} 
            // stat on aligned reads
            if (params.fastqc && !params.cram){
                fastqc_ali_nucmer(nucmer_ali, "fastqc/nucmer", "nucmer")
                logs.concat(fastqc_ali_nucmer.out).set{logs} // save log
            }
            if (params.samtools_stats){
                samtools_stats_ali_nucmer(nucmer_ali, reference.collect(), "samtools_stats/nucmer", "nucmer")
                logs.concat(samtools_stats_ali_nucmer.out).set{logs} // save log
            }
        }
        
        // ------------------- salmon -----------------
        if ("salmon" in aligner_list ){
            // index
            if (! salmon_index_done){ // run salmon index only if not already done when libray type is guessed
                salmon_index(reference.collect(), "alignment/salmon/indicies")
            }
            // align
            salmon(reads, salmon_index.out.collect(), "alignment/salmon")
            logs.concat(salmon.out.salmon_summary).set{logs} // save log
            // convert sam to bam
            samtools_sam2bam_salmon(salmon.out.tuple_sample_sam)
            // sort and convert to cram
            samtools_sort_salmon(samtools_sam2bam_salmon.out.tuple_sample_bam, reference.collect())
            // index
            samtools_index_salmon(samtools_sort_salmon.out, "alignment/salmon")
            samtools_index_salmon.out.tuple_sample_ali.set{salmon_ali} // set name
            // save aligned reads
            sorted_ali.concat(salmon_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_salmon(salmon_ali, "fastqc/salmon", "salmon")
                logs.concat(fastqc_ali_salmon.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_salmon(salmon_ali, reference.collect(), "samtools_stats/salmon", "salmon")
                logs.concat(samtools_stats_ali_salmon.out).set{logs} // save log
            }
        }

        // ------------------- STAR ----------------- output is sorted bam but bai missing
        if ( "star" in aligner_list ){
            // Take read length information from only one sample for --sjdbOverhang option
            // only needed if --sjdbFileChrStartEnd or --sjdbGTFfile option is activated)
            first_file = reads.first()
            prepare_star_index_options(first_file, annotation.collect())
            star_index(reference.collect(), prepare_star_index_options.out, annotation.collect(), "alignment/star/indicies") // index
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
            // filter unmapped reads if requested
            if(params.filter_unmapped){
                samtools_view_filter_star(star_ali)
                samtools_view_filter_star.out.tuple_sample_bam.set{star_ali}
            }
            // convert to cram if requested
            if(params.cram){
                samtools_bam2cram_star(star_ali, reference.collect())
                samtools_bam2cram_star.out.tuple_sample_ali.set{star_ali}
            } 
            // index bam files
            samtools_index_star(star_ali, "alignment/star")
            samtools_index_star.out.tuple_sample_ali.set{star_ali}
            // save aligned reads
            sorted_ali.concat(star_ali).set{sorted_ali} 
            // stat on aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_star(star_ali, "fastqc/star", "star")
                logs.concat(fastqc_ali_star.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_star(star_ali, reference.collect(), "samtools_stats/star", "star")
                logs.concat(samtools_stats_ali_star.out).set{logs} // save log
            }
        }

        // ---------------- subread ----------------- output is sorted bam with bai
        if ( "subread" in aligner_list ){
            // index
            subread_index(reference.collect(), "alignment/subread/indicies")
            // align
            subread(reads, reference.collect(), subread_index.out.collect(), annotation.collect(), "alignment/subread")
            subread.out.tuple_sample_bam.set{subread_ali} // set name
            // filter unmapped reads if requested
            if(params.filter_unmapped){
                samtools_view_filter_subread(subread_ali)
                samtools_view_filter_subread.out.tuple_sample_bam.set{subread_ali}
            }
            // convert to cram if requested
            if(params.cram){
                samtools_bam2cram_subread(subread_ali, reference.collect())
                samtools_bam2cram_subread.out.tuple_sample_ali.set{subread_ali}
            } 
            // index bam
            samtools_index_subread(subread_ali, "alignment/subread")
            samtools_index_subread.out.tuple_sample_ali.set{subread_ali}
            // save aligned reads
            sorted_ali.concat(subread_ali).set{sorted_ali}
            // stat on sorted aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_subread(subread_ali, "fastqc/subread", "subread")
                logs.concat(fastqc_ali_subread.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_subread(subread_ali, reference.collect(), "samtools_stats/subread", "subread")
                logs.concat(samtools_stats_ali_subread.out).set{logs} // save log
            }
        }

        // ---------------- sublong -----------------
        if ( "sublong" in aligner_list ){
            // index
            sublong_index(reference.collect(), "alignment/sublong/indicies")
            // align
            sublong(reads, reference.collect(), sublong_index.out.collect(), "alignment/sublong")
            // merge bam if paired
            samtools_merge_bam_if_paired(sublong.out.tuple_sample_bam)
            samtools_merge_bam_if_paired.out.tuple_sample_bam.set{sublong_ali_tmp} // set name
            // sort and convert to cram
            samtools_sort_sublong(sublong_ali_tmp, reference.collect())
            // index bam/cram files
            samtools_index_sublong(samtools_sort_sublong.out, "alignment/sublong")
            samtools_index_sublong.out.tuple_sample_ali.set{sublong_ali} // set name
            // save aligned reads
            sorted_ali.concat(sublong_ali).set{sorted_ali}
            // stat on sorted aligned reads
            if(params.fastqc && !params.cram){
                fastqc_ali_sublong(sublong_ali, "fastqc/sublong", "sublong")
                logs.concat(fastqc_ali_sublong.out).set{logs} // save log
            }
            if(params.samtools_stats){
                samtools_stats_ali_sublong(sublong_ali, reference.collect(), "samtools_stats/sublong", "sublong")
                logs.concat(samtools_stats_ali_sublong.out).set{logs} // save log
            }
        }
        // ------------------- MULTIQC -----------------
        multiqc(logs.collect(),params.multiqc_config)

        // ------------------- R rendering -----------------
        r_rendering(multiqc.out.multiqc_report_data)

    emit:
        sorted_ali                // channel: [ val(meta), path(alignment), path(index) ]

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
        aligned to a reference using one or more aligners, the ouput is converted in bam in needed, and finaly 
        sorted. If the fastqc option is activated, fastqc is run on raw and aligned reads. A multiqc report is generated containing
        fastqc output and supported aligner output.

        Usage example:
        nextflow run aline.nf --reads /path/to/reads_{1,2}.fastq.gz --reference /path/to/reference.fa --outdir alignment_results --aligner bbmap,bowtie2 --fastqc true

        --help                      prints the help section

    General Parameters
        --reads                     path to the reads file, folder or csv. If a folder is provided, all the files with proper extension in the folder will be used. You can provide remote files (commma separated list).
                                    file extension expected : <.fastq.gz>, <.fq.gz>, <.fastq> or <.fq> 
                                                              for paired reads extra <_R1_001> or <_R2_001> is expected where <R> and <_001> are optional. e.g. <sample_id_1.fastq.gz>, <sample_id_R1.fastq.gz>, <sample_id_R1_001.fastq.gz>)
                                    csv input expects 6 columns: sample, fastq_1, fastq_2, strandedness, read_type and data_type. 
                                    fastq_2 is optional and can be empty. Strandedness, read_type and data_type expect same values as corresponding AliNe parameters; If a value is provided via AliNe paramter, it will override the value in the csv file.
                                    Example of csv file:
                                        sample,fastq_1,fastq_2,strandedness,read_type,data_type
                                        control1,path/to/data1.fastq.gz,,auto,short_single,rna
                                        control2,path/to/data2_R1.fastq.gz,path/to/data2_R2.fastq.gz,auto,short_paired,rna
        --reference                 path to the reference file (fa, fa.gz, fasta or fasta.gz)
        --aligner                   aligner(s) to use among this list (comma or space separated) ${align_tools}
        --outdir                    path to the output directory (default: alignment_results)
        --annotation                [Optional][used by STAR, Tophat2] Absolute path to the annotation file (gtf or gff3)

    Type of input reads
        --data_type                 type of data among this list ${data_type_allowed} (no default)
        --read_type                 type of reads among this list ${read_type_allowed} (no default)
        --strandedness              Set the strandedness of your reads (no default). In auto mode salmon will guess the library type for each sample.
                                    If you know the library type you can set it to one of the following: ${strandedness_allowed}. See https://salmon.readthedocs.io/en/latest/library_type.html for more information.
                                    In such case the sample library type will be used for all the samples.
        --read_length               Length of the reads, if none provided it is automatically deduced. [Optional, used by STAR, Salmon, Kallisto]
        --seqtk_sample_size         Number of reads to sample from the input reads (default: 10000). The subsamples are used to guess strandedness and read_length. 
                                    Note that this step is not needed if the library type and strandedness are provided via AliNe parameters.

    Extra steps 
        --trimming_fastp            run fastp for trimming (default: false)
        --fastqc                    run fastqc on raw and aligned reads (default: false). Note: FastQC will be automatically disabled for alignment files when --cram is enabled.
        --samtools_stats            run samtools stats on aligned reads (default: false)
        --filter_unmapped           filter out unmapped reads from final alignment files (default: false). Filtering is performed during sorting when possible for optimal performance.
        --cram                      output alignment files in sorted CRAM format instead of sorted BAM (default: false). This saves disk space but disables FastQC on alignment files. Conversion is performed during sorting when possible for optimal performance.
        --multiqc_config            path to the multiqc config file (default: config/multiqc_conf.yml)

    Aligner specific options
        --bbmap_options             additional options for bbmap
        --bowtie_options            additional options for bowtie
        --bowtie2_options           additional options for bowtie2
        --bwaaln_options            additional options for bwaaln
        --bwamem_options            additional options for bwamem
        --bwamem2_options           additional options for bwamem2
        --bwasw_options             additional options for bwasw
        --dragmap_options           additional options for dragmap
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
        --subread_options           additional options for subread
        --sublong_options           additional options for sublong

    Other
        --monochrome_logs           set to true to disable color in logs (default: false)
        --relax                     set to true to disable aligner parameter changes made by AliNe (default: false) 

    """
}

def printAlignerOptions(aligner_list) {
    def sentence = ""
    if ("bbmap" in aligner_list){ 
        sentence += """
    bbmap parameters
        bbmap_options           : ${params.bbmap_options}
    """} 
    if ("bowtie" in aligner_list){ 
        sentence += """       
    bowtie parameters
        bowtie_options          : ${params.bowtie_options}
    """} 
    if ("bowtie2" in aligner_list){ 
        sentence += """       
    bowtie2 parameters
        bowtie2_options         : ${params.bowtie2_options}
    """} 
    if ("bwaaln" in aligner_list){
        sentence += """
    bwaaln parameters
        bwa_options             : ${params.bwaaln_options}
    """} 
    if ("bwamem" in aligner_list){
        sentence += """
    bwamem parameters
        bwamem_options          : ${params.bwamem_options}
    """} 
    if ("bwamem2" in aligner_list){
        sentence += """
    bwamem2 parameters
        bwamem2_options         : ${params.bwamem2_options}
    """} 
    if ("bwasw" in aligner_list){
        sentence += """
    bwasw parameters
        bwasw_options           : ${params.bwasw_options}
    """} 
    if ("dragmap" in aligner_list){
        sentence += """
    dragmap parameters
        dragmap_options         : ${params.dragmap_options}
    """} 
    if ("graphmap2" in aligner_list){
        sentence += """
    graphmap2 parameters
        graphmap2_options       : ${params.graphmap2_options}
    """} 
    if ("hisat2" in aligner_list){
        sentence += """
    hisat2 parameters
        hisat2_options          : ${params.hisat2_options}
    """}
    if ("kallisto" in aligner_list){
        sentence += """
    kallisto parameters
        kallisto_options        : ${params.kallisto_options}
        kallisto_index_options  : ${params.kallisto_index_options}
    """}
    if ("minimap2" in aligner_list){
        sentence += """
    minimap2 parameters
        minimap2_options        : ${params.minimap2_options}
        minimap2_index_options  : ${params.minimap2_index_options}
    """} 
    if ("ngmlr" in aligner_list){
        sentence += """
    ngmlr parameters
        ngmlr_options           : ${params.ngmlr_options}
    """} 
    if ("novoalign" in aligner_list){
        sentence += """
    novalign parameters
        novalign_options        : ${params.novoalign_options}
        novoalign_license       : ${params.novoalign_license}
    """} 
    if ("nucmer" in aligner_list){
        sentence += """
    nucmer parameters
        nucmer_options          : ${params.nucmer_options}
    """}
     if ("salmon" in aligner_list){
        sentence += """
    salmon parameters
        salmon_index_options    : ${params.salmon_index_options}
        salmon_options          : ${params.salmon_options}
    """}
    if ("star" in aligner_list){
        sentence += """
    star parameters
        star_index_options      : ${params.star_index_options}
        star_options            : ${params.star_options}      
        star_2pass              : ${params.star_2pass}
    """}
    if ("subread" in aligner_list){
        sentence += """
    subread parameters
        subread_options         : ${params.subread_options}
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
 - dragmap 

How to add an aliger ?
-----------------------
Add a module,
In aline.nf 
    Add the aligner in the aligner_list
    Import module here
    If behaves differently depending sequencing techniology add condition in the PARAMS CHECK sections (module bash.nf)
    If tool need read length guessing add it in the LIBRARY TYPE GUESSING to activate the guessing if read length not provided by user.
    Add in help
    Add in printAlignerOptions
    Add aligner process in the main workflow 
    Think to convert sam to bam if necessary
    Think to sort bam output if necessary 
Add info in README.md
Add tool in multiqc config file (at least for fastqc)
*/
