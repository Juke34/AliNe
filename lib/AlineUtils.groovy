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
}