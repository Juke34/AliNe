// dedicated to the utility functions
class AlineUtils {

    // Function to chek if we
    public static Boolean is_url(file) {
        
        // check if the file is a URL
        if (file =~ /^(http|https|ftp|s3|az|gs):\/\/.*/) {
            return true
        } else {
            return false
        }
    }
    // Function to chek if we
    public static Boolean is_fastq(file) {
        
        // check if the file is a URL
        if (file =~ /.*(fq|fastq|fq.gz|fastq.gz)$/) {
            return true
        } else {
            return false
        }
    }

    // Function to extract the basename of a file
    public static String getCleanName(file) {
        def fileClean = file[0].baseName.replaceAll(/\.(gz)$/, '') // remove .gz
        fileClean = fileClean.replaceAll(/\.(fasta|fa)$/, '') // remove .fasta or .fa
        fileClean = fileClean.replaceAll(/\.(fastq|fq)$/, '') // remove .fastq or .fq
        return fileClean
    }

    /*
     * Extract sample ID from filename (nf-core style)
     * Removes file extensions and paired-end suffixes based on read type
     * 
     * @param filename The filename to process (can be String or File object)
     * @param read_type The read type: "short_paired", "short_single", "pacbio", "ont", etc.
     * @return The extracted sample ID
     * 
     * Examples:
     *   Paired-end:
     *     sample_ABC_R1_001.fastq.gz -> sample_ABC
     *     control_1.fq.gz -> control
     *   Single-end:
     *     sample.fastq.gz -> sample
     *     test.bam -> test
     *
     * This function must be synchronized with the counterpart in RAIN!
     */
    public static String get_file_id(filename, read_type = "short_paired") {
        def name = new File(filename.toString()).getName()

        if (read_type == "short_paired") {
            // Pattern for paired-end: _R?[12](_\d+)?(.fastq|.fq)(.gz)?$
            return name
                .replaceAll(/_R?[12](_\d+)?\.fastq\.gz$/, '')  // Remove _R1_001.fastq.gz, _1_001.fastq.gz, etc.
                .replaceAll(/_R?[12](_\d+)?\.fq\.gz$/, '')     // Remove _R1_001.fq.gz, etc.
                .replaceAll(/_R?[12](_\d+)?\.fastq$/, '')      // Remove _R1_001.fastq, etc.
                .replaceAll(/_R?[12](_\d+)?\.fq$/, '')         // Remove _R1_001.fq, etc.
                .replaceAll(/_R?[12](_\d+)?\.bam$/, '')        // Remove _R1.bam, etc.
        } else {
            // Pattern for single-end: (.fastq|.fq)(.gz)?$
            return name
                .replaceAll(/\.fastq\.gz$/, '')  // Remove .fastq.gz
                .replaceAll(/\.fq\.gz$/, '')     // Remove .fq.gz
                .replaceAll(/\.fastq$/, '')      // Remove .fastq
                .replaceAll(/\.fq$/, '')         // Remove .fq
                .replaceAll(/\.bam$/, '')        // Remove .bam
        }
    } 
}