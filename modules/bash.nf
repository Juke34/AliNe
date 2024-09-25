/*
Here are described all processes related to bash
*/

// Deal with annotation file
// annotation is used by different aligner (tophat2, star, etc.). To avoid to duplicate processes according to the presence of the annotation file, a specific process is dedicated to create a fake file is none provided. 
// If process receive a file wich is not the fake one it includes the file in the command. To append the options of aligner we will use the annotation_file variable
// While the processes will be called sending the "annotation" channel created by the prepare_annotation process.
process prepare_annotation {
    label 'seqtk'
    
    input:
        val (annotation)
        
    output:
        path("*"), emit: annotation

    script:

        // set input/output according to short_paired parameter
        if (params.annotation){
            """
            ln -s ${params.annotation}
            """
        } else {
            """
            touch null.gtf
            """
        }
            
}