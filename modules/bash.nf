/*
Here are described all processes related to bash
*/

// A process to compute the mean read length of a FASTQ
process read_length {
    label 'bash'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*", mode: 'copy'

    input:
        tuple val(meta), path(fastq)
        val outpath
        
    output:
        tuple val(meta), env(READLENGTH), emit: tuple_id_readlength
        path "*read_length.txt"

    script:
        // compute made only on first file when is paired
        """
        READLENGTH=\$(cat ${fastq[0]} | awk 'NR % 4 == 2 {sum += length(\$0); count++} END {if (count > 0) print int((sum / count) + 0.5); else print "0"}')
        echo \$READLENGTH > ${fastq[0].baseName.replace('.fastq','')}_read_length.txt
        """
}

/*
This process checks the aligner used according to the read type and inform user if non-optimal choice is made.
*/
process check_aligner{
    label 'bash'
    tag "${meta.id}"

    input:
        tuple val(meta), path(fastq)
        val aligner_list

    output:
        tuple val(meta), path(fastq)

    script:

        // --- bbmap tool ---
        if ( "bbmap" in aligner_list ){
            // short and long reads ok
        }

        // --- bowtie tool ---
        if ( "bowtie" in aligner_list ){
            if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Bowtie aligner is not recommended to align long reads!"
            }
        }

        // --- bowtie2 tool ---
        if ( "bowtie2" in aligner_list ){
            if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Bowtie2 aligner is not recommended to align long reads!"
            }
        }

        // --- bwa aln tool ---
        if ( "bwaaln" in aligner_list ){
            if ( meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Bwaaln aligner is not recommended to align long reads!"
            }
        }
        

        // --- bwa mem tool ---
        if ( "bwamem" in aligner_list ){
            // short and long reads ok
        }

        // --- bwa mem2 tool ---
        if ("bwamem2" in aligner_list ){
            if ( meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Bwamem2 aligner is not recommended to align long reads!"
            }
        }
       

        // --- bwa sw tool ---
        if ( "bwasw" in aligner_list ){
            if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Bwasw aligner is not recommended to align long reads!"
            }
        }

        // --- graphmap2 tool ---
        if ( "graphmap2" in aligner_list ){
            if ( meta.read_type == "short_single" && meta.read_type == "short_paired"){
                log.info "${meta.id} => Graphmap2 aligner is not recommended to align short reads!"
            }
        }
        
        // hisat2
        if ("hisat2" in aligner_list ){
            if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Hisat2 aligner is not recommended to align long reads!"
            }
        }

        // --- kallisto tool ---
        if ( "kallisto" in aligner_list ){
            if ( meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Kallisto aligner is not recommended to align long reads!"
            }
        }

        // --- last tool ---
        if ( "last" in aligner_list ){
            if ( meta.read_type == "short_single" && meta.read_type == "short_paired"){
                log.info "${meta.id} => Last aligner is not recommended to align short reads!"
            }
        }

        // ---- minimap2 tool ---
        // Force -a option to be sure to get sam output
        if ("minimap2" in aligner_list ){
            // short and long reads ok
        }

        // ngmlr tool - check options
        if ("ngmlr" in aligner_list ){
            if ( meta.read_type == "short_single" && meta.read_type == "short_paired"){
                log.info "${meta.id} => Ngmlr aligner is not recommended to align short reads!"
            }
        }

        // novoalign tool - load license into the container
        if ("novoalign" in aligner_list ){
            if ( meta.read_type == "ont"){
                log.info "${meta.id} => Novoalign aligner is not recommended to align long reads!"
            }
        }

        // mummer / nucmer
        if ("nucmer" in aligner_list ){
            if ( meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Nucmer aligner is not recommended to align long reads!"
            }
        }

        // --- salmon tool ---
        if ( "salmon" in aligner_list ){
            if ( meta.read_type == "ont" || meta.read_type == "pacbio"){
                log.info "${meta.id} => Salmon aligner is not recommended to align long reads!"
            }
        }

        // --- star tool ---
        if ( "star" in aligner_list ){
        }
        
        // --- subread tool ---
        if ( "subread" in aligner_list ){
            if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                log.info "${meta.id} => Subread aligner is not recommended to align long reads!"
            }
        }

        // --- sublong tool ---
        if ( "sublong" in aligner_list ){
            if ( meta.read_type == "short_single" || meta.read_type == "short_paired"){
                log.info "${meta.id} => Sublong aligner is not recommended to align short reads!"
            }
        }

        """
        echo "done"
        """
}

/*
This process checks the aligner parameters and adpat them according to the read type and annotation and data type.
*/
process check_aligner_params{
    label 'bash'
    tag "${meta.id}"
    publishDir "${params.outdir}/${outpath}", pattern: "*.txt", mode: 'copy'

    input:
        tuple val(meta), path(fastq)
        val aligner_list
        val annotation
        val outpath

    output:
        tuple val(meta), path(fastq)
        path "*.txt"

    script:

        // --- bbmap tool ---
        if ( "bbmap" in aligner_list ){
            def bbmap_tool = "bbmap.sh"
            def bbmap_options = params.bbmap_options ?: ""
            if ( !params.relax ){
                if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                    bbmap_tool = "mapPacBio.sh"
                    log.info "${meta.id} => Long reads being used, using mapPacBio.sh to align with bbmap!\n" +
                             "    However, if you know what you are doing you can activate the AliNe --relax parameter to use bbmap.sh anyway."
                    // Function to check and set maxlen in params.bbmap_options when long_reads is set
                    // params is supposed to be a immutable. Using params.replace method might not be supported in the future 
                    if ( !bbmap_options.contains("maxlen") ){
                        bbmap_options += " maxlen=5000"
                    }
                }
            }
            // set library type
            if (! bbmap_options.contains(" xs=") && 
                meta.paired && 
                meta.strandedness){ 
                if (meta.strandedness.contains("U") ){ // this test must be before the others
                    bbmap_options += " -xs=us"
                }  
                else if (meta.strandedness.contains("I") ){
                    bbmap_options += " xs=fr"
                } else if (meta.strandedness.contains("O") ){
                    bbmap_options += " xs=ss"
                } 
            }
            // set read orientation
            if (! bbmap_options.contains(" rcs=") &&
                ! bbmap_options.contains(" requirecorrectstrand=") && 
                meta.paired &&
                meta.strandedness &&
                ! meta.strandedness.contains("I") 
            ){ 
                bbmap_options += " rcs=f"
            }

            meta.bbmap_tool = bbmap_tool
            meta.bbmap_options = bbmap_options
        }

        // --- bowtie tool ---
        if ( "bowtie" in aligner_list ){
            def bowtie_options = params.bowtie_options ?: ""
            if (! bowtie_options.contains("--fr ") && 
                ! bowtie_options.contains("--rf ") && 
                ! bowtie_options.contains("--ff ") &&
                meta.paired && meta.strandedness) { 
                if (meta.strandedness.contains("I") ){
                    bowtie_options += " --fr"
                } else if (meta.strandedness.contains("O") ){
                    bowtie_options += " --rf"
                } else if (meta.strandedness.contains("M") ){
                    bowtie_options +=  " --ff"
                }  
            }
            meta.bowtie_options = bowtie_options
        }

        // --- bowtie2 tool ---
        if ( "bowtie2" in aligner_list ){
            def bowtie2_options = params.bowtie2_options ?: ""
            if (! bowtie2_options.contains("--fr ") && 
                ! bowtie2_options.contains("--rf ") && 
                ! bowtie2_options.contains("--ff ") &&
                meta.paired && 
                meta.strandedness){ 
                if (meta.strandedness.contains("I") ){
                    bowtie2_options += " --fr"
                } else if (meta.strandedness.contains("O") ){
                    bowtie2_options += " --rf"
                } else if (meta.strandedness.contains("M") ){
                    bowtie2_options += " --ff"
                }  
            }
            meta.bowtie2_options = bowtie2_options
        }

        // --- bwa aln tool ---
        if ( "bwaaln" in aligner_list ){
            def bwaaln_options = params.bwaaln_options ?: ""
            meta.bwaaln_options = bwaaln_options
        }
        

        // --- bwa mem tool ---
        if ( "bwamem" in aligner_list ){
            def bwamem_options = params.bwamem_options ?: ""
            if( !params.relax ){
                if (meta.read_type == "pacbio"){
                    if ( ! bwamem_options.contains(" pacbio") ){
                        bwamem_options += " -x pacbio"
                        log.info "${meta.id} => Pacbio reads being used, setting -x pacbio to bwamem!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
                if (meta.read_type == "ont"){
                    if ( ! bwamem_options.contains(" ont2d") ){
                        bwamem_options += " -x ont2d"
                        log.info "${meta.id} => Ont reads being used, setting -x ont2d to bwamem!\n" + 
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
            }
            meta.bwamem_options = bwamem_options
        }

        // --- bwa mem2 tool ---
        if ("bwamem2" in aligner_list ){
            def bwamem2_options = params.bwamem2_options ?: ""
            if ( !params.relax ){
                if (meta.read_type == "pacbio"){
                    if ( ! bwamem2_options.contains(" pacbio") ){
                        bwamem2_options += " -x pacbio"
                        log.info "${meta.id} => Pacbio reads being used, setting -x pacbio to bwamem2!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
                if (meta.read_type == "ont"){
                    if ( ! bwamem2_options.contains(" ont2d") ){
                        bwamem2_options += " -x ont2d"
                        log.info "${meta.id} => Ont reads being used, setting -x ont2d to bwamem2!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
            }
            meta.bwamem2_options = bwamem2_options
        }

        // --- bwa sw tool ---
        if ( "bwasw" in aligner_list ){
            def bwasw_options = params.bwasw_options ?: ""
            meta.bwasw_options = bwasw_options
        }

        // --- graphmap2 tool ---
        if ( "graphmap2" in aligner_list ){
            def graphmap2_options = params.graphmap2_options ?: ""
            if (meta.annotation && !graphmap2_options.contains("--gtf ") ){
                graphmap2_options += " --gtf ${annotation}"
            }
            meta.graphmap2_options = graphmap2_options
        }
        
        // hisat2
        if ("hisat2" in aligner_list ){
            def hisat2_options = params.hisat2_options ?: ""
            // deal with library type - default is unstranded.
            if (! hisat2_options.contains("--rna-strandness ") &&
                meta.strandedness && ! meta.strandedness.contains("U")
                ){ // only if -S is not set and if we are not skipping library usage
                if ( meta.paired ) {
                    if (meta.strandedness.contains("SR") ) {
                        hisat2_options += " --rna-strandness RF"
                    } else if (meta.strandedness.contains("SF") ) {
                        hisat2_options += " --rna-strandness FR"
                    }
                } // unpair
                else { // Lib is not unstranded 
                    if ( meta.strandedness.contains("SF") ) {
                        hisat2_options += " --rna-strandness F"
                    } else if ( meta.strandedness.contains("SR") ) {
                        hisat2_options += " --rna-strandness R"
                    } 
                }
            }
            // read_orientation is only for paired-end reads
            if (! hisat2_options.contains("--fr ") && 
                ! hisat2_options.contains("--rf ") && 
                ! hisat2_options.contains("--ff ") && 
                meta.paired && meta.strandedness ){ 
                if ( meta.strandedness.contains("I") ){
                    hisat2_options += " --fr"
                } else if ( meta.strandedness.contains("O") ){
                    hisat2_options += " --rf"
                }
                else if ( meta.strandedness.contains("M") ){
                    hisat2_options += " --ff"
                }  
            }
            meta.hisat2_options = hisat2_options
        }

        // --- kallisto tool ---
        if ( "kallisto" in aligner_list ){
            def kallisto_options = params.kallisto_options ?: ""
            // deal with read_orientation
            if (! kallisto_options.contains("--fr-stranded ") && 
                ! kallisto_options.contains("--rf-stranded ") && 
                meta.strandedness) { 
                
                if (meta.strandedness.contains("I") ){
                    kallisto_options += " --fr-stranded"
                } else if (meta.strandedness.contains("O") ){
                    kallisto_options += " --rf-stranded"
                } 
            }
            // Use read length (-l) and sd (-s) from params? for single-end reads
            if(! meta.paired){
                if ( !kallisto_options.contains("-l ") ){
                    kallisto_options += " -l ${meta.read_length}"
                }
                if ( !kallisto_options.contains("-s ") ){
                    // 10% of read length will be used as Estimated standard deviation of fragment length
                    def tenPercent = (meta.read_length.toInteger() * 10 / 100) as int 
                    kallisto_options += " -s ${tenPercent}"
                }
            }
            meta.kallisto_options = kallisto_options
        }

        // --- last tool ---
        if ( "last" in aligner_list ){
            def last_options = params.last_options ?: ""
            meta.last_options = last_options
        }

        // ---- minimap2 tool ---
        // Force -a option to be sure to get sam output
        if ("minimap2" in aligner_list ){
            def minimap2_options = params.minimap2_options ?: ""
            if ( !params.relax ){
                if (meta.read_type == "short_single" || meta.read_type == "short_paired"){
                    if ( ! minimap2_options.contains("--sr ") ){
                        minimap2_options += " --sr"
                    }
                }
                if (meta.read_type == "pacbio"){
                    if ( ! minimap2_options.contains(" ava-pb") and ! minimap2_options.contains(" splice:hq") and 
                         ! minimap2_options.contains(" map-hifi") and ! minimap2_options.contains(" map-pb") ){
                         log.info("${meta.id} => Warn: <${minimap2_options}> minimap2 options missing or not accepted for pacbio data.\n" +
                                 "    We set the default <map-pb> parameter. If you do not agree, please provide options among this list:\n" +
                                 "    ava-pb, splice:hq, map-hifi, map-pb (see https://github.com/lh3/minimap2).\n" +
                                 "    If you wish to use parameter not intended for pacbio data use the --relax parameter to skip this warning message.")
                        minimap2_options += " -x map-pb"
                    }
                }
                if (meta.read_type == "ont"){
                    if ( ! minimap2_options.contains(" ava-ont") and ! minimap2_options.contains(" splice") and 
                        ! minimap2_options.contains(" lr:hq") and ! minimap2_options.contains(" map-ont") ){
                        log.info("${meta.id} => Warn: <${minimap2_options}> minimap2 options missing or not accepted for ont data.\n" +
                                 "    We set the default <map-ont> option. If you do not agree, please provide options among this list:\n" +
                                 "    ava-ont, splice, lr:hq, map-ont (see https://github.com/lh3/minimap2).\n" +
                                 "    If you wish to use parameter not intended for pacbio data use the --relax parameter to skip this warning message.")
                        minimap2_options += " -x map-ont"
                    }
                }
            }
            // relax or not this option has to be used
            if ( ! minimap2_options.contains("-a ") ){
                minimap2_options += " -a"
            }
            meta.minimap2_options = minimap2_options
        }

        // ngmlr tool - check options
        if ("ngmlr" in aligner_list ){
            def ngmlr_options = params.ngmlr_options ?: ""
            if (!params.relax) {
                // for pacbio reads, set -g 20 and -x 0
                if (meta.read_type == "ont"){
                    if (! ngmlr_options.contains("-x ont")){
                        ngmlr_options += " -x ont"
                        log.info "${meta.id} => Ont reads being used, setting -x ont to ngmlr!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
            }
            meta.ngmlr_options = ngmlr_options
        }
       

        // novoalign tool - load license into the container
        if ("novoalign" in aligner_list ){
            def novoalign_lic = ""
            def novoalign_options = params.novoalign_options ?: ""
            if (!params.relax) {
                // for pacbio reads, set -g 20 and -x 0
                if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                    if (! novoalign_options.contains("-g ")){
                        novoalign_options += " -g 20"
                        log.info "${meta.id} => Long reads being used, setting -g 20 to Novoalign!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                    if (! novoalign_options.contains("-x ")){
                        novoalign_options += " -x 0"
                        log.info "${meta.id} => Long reads being used, setting -x 0 to Novoalign!\n" +
                                 "    However, if you know what you are doing you can activate the AliNe --relax parameter and avoid this behavior."
                    }
                }
            }
            meta.novoalign_options = novoalign_options
        }

        // mummer / nucmer
        if ("nucmer" in aligner_list ){
            def nucmer_options = params.nucmer_options ?: ""
            meta.nucmer_options = nucmer_options
        }

        // --- salmon tool ---
        if ( "salmon" in aligner_list ){
            def salmon_options = params.salmon_options ?: ""
            // deal with library type 
            if (! salmon_options.contains("-l ") && ! salmon_options.contains("--libType ") ){
                if (meta.strandedness){ 
                    salmon_options += " -l ${meta.strandedness}"
                } else {
                    salmon_options += " -l A" // A for automatic
                }
            }
            if ( !meta.paired){
                // Use read length (--fldMean) and sd (--fldSD) from params? only for single-end reads
                if ( !salmon_options.contains("--fldMean ") ){
                    salmon_options += " --fldMean ${meta.read_length}"
                }
                if ( !salmon_options.contains("--fldSD ") ){
                    // 10% of read length will be used as Estimated standard deviation of fragment length
                    def tenPercent = (meta.read_length.toInteger() * 10 / 100) as int 
                    salmon_options += " --fldSD ${tenPercent}"
                }
            }
            meta.salmon_options = salmon_options
        }

        // --- star tool ---
        def star_tool = "STAR"
        def star_options = params.star_options ?: ""
        if ( "star" in aligner_list ){
            if (meta.annotation && ! star_options.contains("--sjdbGTFfile ") ){
                star_options += " --sjdbGTFfile ${annotation}"
            }
            if (!params.relax){
                if (meta.read_type == "pacbio" || meta.read_type == "ont"){
                    star_tool = "STARlong"
                }
            }
            meta.star_tool = star_tool
            meta.star_options = star_options
        }
        
        // --- subread tool ---
        if ( "subread" in aligner_list ){
            def subread_options = params.subread_options ?: ""
            // deal with annotation
            if (meta.annotation && !subread_options.contains("-a ") ){
                subread_options += " -a ${annotation}"
            }
            // -t specifes the type of input sequencing data. Possible values include 0, denoting RNA-seq data, or 1, denoting genomic DNA-seq data. For genomic DNA-seq data, the aligner takes into account both the number of matched bases and the number of mis-matched bases to determine the the best mapping location after applying the ‘seed-and-vote’ approach for read mapping. For RNA-seq data, only the number of mis-matched bases is considered for determining the best mapping location.
            if (meta.data_type.toString().toLowerCase() == "dna" && !subread_options.contains("-t 1") ){
                subread_options += "-t 1"
            }
            if (meta.data_type.toLowerCase() == "rna" && !subread_options.contains("-t 0") ){
                subread_options += " -t 0"
            }
            // deal with library type orientation
            if (! subread_options.contains("-S ") &&
                meta.paired && meta.strandedness) { // only if -S is not set and if we are not skipping library usage
                if (meta.strandedness.contains("M") ){
                    subread_options += " -S ff"
                } else if (meta.strandedness.contains("O") ) {
                    subread_options += " -S rf"
                } else if (meta.strandedness.contains("I") ) {
                    subread_options += " -S fr"
                } 
            }
            meta.subread_options = subread_options
        }

        // --- sublong tool --- dealed apart subread due to different output (sorted or not)
        if ( "sublong" in aligner_list ){
            def sublong_options = params.sublong_options ?: ""
            meta.sublong_options = sublong_options
        }
        
        // ---------------- Display ----------------
        // Create a map of aligner options
        def optionsMap = [:]
        for (tool in aligner_list) {
            def key = "${tool}_options"
            optionsMap[tool] = params.getAt(key) ?: ""
        }

        // Serialize this map as key=value pairs separated by commas (or any delimiter)
        def optionsString = optionsMap.collect { k, v -> "${k}=${v}" }.join(",")


        """
        fileout="${meta.id}.txt"
        meta="${meta}"
        # Remove square brackets
        meta=\$(echo "\$meta" | tr -d '[]')
        # Remove leading and trailing space
        meta=\$(echo "\$meta" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//')

        optionsString='${optionsString}'

        # Parse optionsString into bash associative array - separator is =
        declare -A valueBefore
        IFS=',' read -ra pairs <<< "\$optionsString"
        for pair in "\${pairs[@]}"; do
            key="\${pair%%=*}"
            key=\$(echo "\$key" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//') # remove leading and trailing spaces
            value="\${pair#*=}"
            [ -z "\$value" ] && value="."
            value=\$(echo "\$value" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//') # remove leading and trailing spaces
            valueBefore["\$key"]="\$value"
        done

        # Parse meta into bash associative array - separator is :
        declare -A valueAfter
        IFS=',' read -ra pairs <<< "\$meta"
        for pair in "\${pairs[@]}"; do
            key="\${pair%%:*}"
            key=\$(echo "\$key" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//') # remove leading and trailing spaces
            value="\${pair#*:}"
            [ -z "\$value" ] && value="."
            value=\$(echo "\$value" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//') # remove leading and trailing spaces
            echo "key <\$key>"
            echo "value <\$value>"
            valueAfter["\$key"]="\$value"
        done

        echo -e "# Below are the tool parameters received by AliNe, and the modifications/choices made by Aline according to available information (read type, annotation, etc.)." >> \$fileout
        echo -e "# Some tools have different executable depending on the type of read. Choices made by AliNe are displayed as aligner_tool." >> \$fileout
        echo -e "Parameter\\tBefore\\tAfter" >> \$fileout
        for tool in \${!valueBefore[@]} ; do 
            newKey="\${tool}_options"
            newKey2="\${tool}_tool"
            echo -e "\$tool\\t\${valueBefore[\${tool}]}\\t\${valueAfter[\${newKey}]}" >> \$fileout
            # Print the tool name if it exists in valueAfter
            if [[ -v valueAfter["\$newKey2"] ]]; then
                echo -e "\$newKey2\\t.\\t\${valueAfter[\${newKey2}]}" >> \$fileout
            fi
        done
        """
}