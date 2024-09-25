AliNe (Alignment in Nextflow) 
=========================================  
<img align="right" src="img/IRD.png" width="200" height="66" /> <img align="right" src="img/MIVEGEC.png" width="100" height="66" />

AliNe is a pipeline written in nextflow that aims to efficiently align reads against a reference genome using the tools of your choice.

Genome + Reads => FastQC -> Alignment -> Sort -> MultiQC

## Table of Contents

   * [Foreword](#foreword)
   * [Installation](#installation)
      * [AliNe](#aline)
      * [Nextflow](#nextflow)
      * [Container platform](#container-platform)
        * [Docker](#docker)
        * [Singularity](#singularity)  
   * [Usage and test](#usage)
   * [Parameters](#parameters)
   * [Uninstall](#uninstall)

## Foreword

AliNe is a pipeline written in nextflow that aims to efficienlty align reads against a reference genome.

A QC with FastQC is made at each step if option activated.
A trimming is feasible before alignment if option activated.
**The pipeline deals automaticallu with all quality encoding ('sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+').** All fastq will be standardised in Phred+33 for downstream alignments by seqkit. 
You can choose to run one or several aligner in parallel.

Here is the list of implemented aligners:

| Tool	| Single End (short reads) | Paired end (short reads) | Pacbio | ONT |
| --- | --- | --- |  --- | --- |
| bbmap | X | X | x | x |
| bowtie2 | X | X | x | x |
| bwaaln | X | X R1 and R2 independently aligned then merged with bwa sampe | X | X |
| bwamem | X | X | x | x |
| bwasw | X | X | x | x |
| graphmap2 | x | x R1 and R2 independently aligned then merged with cat | X | X |
| hisat2 | x | x | x | x |
| minimap2 | x | x | X | X |
| ngmlr | x | x | X | X |
| novoalign | X | X | X | ? |
| nucmer | X | X R1 and R2 are concatenated then aligned | x | x |
| star | X | X | x | x |
| star 2pass mode | X | X | x | x |
| subread | X | X | x | x |
| sublong |  |  |  |  |
| tophat | X | X | na | na |

It is possible to bypass the default authorized read type using the AliNe --relax parameter.

The pipeline deals automatically with the library types. It extract 10 000 reads by default and-d run salmon to guess the library type. 
It is then translated to the correct option in the following aligners:

| Tool	| tool option | Library type by salmon | Comment | 
| --- | --- | --- | --- |
| bbmap | xs=fr / xs=ss / xs=us | ISF ISR / OSF OSR / U | strand information |
| bbmap | - / rcs=f / | ISF ISR IU / OSF OSR OU MSF MSR MU | read orientation |
| bowtie2 | --fr / --rf / --ff |  ISF ISR IU / OSF OSR OU / MSF MSR MU| read orientation |
| bwaaln | na | na | na | |
| bwamem | na | na | na | |
| bwasw | na | na | na | |
| graphmap2 | na | na | na | |
| hisat2 | --rna-strandness [ F / R / FR / RF ] | SF / SR / ISF OSF MSF / ISR OSR MSR | strand information |
| hisat2 | --fr / --rf / --ff | I / O / M | read orientation |
| minimap2 | na | na | |
| ngmlr | na | na | |
| novoalign | na | na | |
| nucmer | na | na | |
| star | na | na | |
| star 2pass mode | na | na | |
| subread | -S fr / -S rf / -S ff | ISF ISR IU / OSF OSR OU / MSF MSR MU | read orientation |
| tophat2 | fr-unstranded / fr-firststrand / fr-secondstrand | U / SR / SF | strand information |

If the skip_libray_usage paramater is set the information provided about the library type provided by the user or guessed by the pipeline via the --library_type parameter is not used.
/!\ If you provide yourself the librairy type via the aligner parameter, it will be used over the information provided or guessed via --library_type.

If you provide an annotation file the pipeline will pass automatically the file to the following aligner:  
| Tool	| accept | 
| --- | --- | 
| bbmap | na | 
| bowtie2 | na |
| bwaaln | na |
| bwamem | na |
| bwasw | na |
| graphmap2 | GTF (--gtf)  |
| hisat2 | na |
| minimap2 | na |
| ngmlr | na |
| novoalign | na |
| nucmer | na |
| star | GTF / GFF ( --sjdbGTFfile + --sjdbGTFtagExonParentTranscript Parent in case of GFF ) |
| star 2pass mode | GTF / GFF (--sjdbGTFfile + --sjdbGTFtagExonParentTranscript Parent in case of GFF ) |
| subread | na |
| tophat | GTF/GFF3 (-G) | 
 

## Installation

The prerequisites to run the pipeline are:  

  * The AliNe repository
  * [Nextflow](https://www.nextflow.io/)  >= 22.04.0
  * [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/)  

### AliNe 

```bash
# clone the workflow repository
git clone https://github.com/Juke34/AliNe.git

# Move in it
cd AliNe
```

### Nextflow 

  * Via conda 

    <details>
      <summary>See here</summary>
      ```
      conda create -n nextflow
      conda activate nextflow
      conda install nextflow
      ```  
    </details>

  * Manually
    <details>
      <summary>See here</summary>
       Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

      ```
      # Make sure 11 or later is installed on your computer by using the command:
      java -version

      # Install Nextflow by entering this command in your terminal(it creates a file nextflow in the current dir):
      curl -s https://get.nextflow.io | bash 

      # Add Nextflow binary to your user's PATH:
      mv nextflow ~/bin/
      # OR system-wide installation:
      # sudo mv nextflow /usr/local/bin
      ```
    </details>

### Container platform

To run the workflow you will need a container platform: docker or singularity.

### Docker

Please follow the instructions at the [Docker website](https://docs.docker.com/desktop/)

### Singularity

Please follow the instructions at the [Singularity website](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)

## Usage

You can first check the available options and parameters by running:
`nextflow run aline.nf --help`

To run the workflow you must select a profile according to the container platform you want to use:   
- `singularity`, a profile using Singularity to run the containers
- `docker`, a profile using Docker to run the containers

The command will look like that: 
```
nextflow run aline.nf -profile docker <rest of paramaters>
```
Another profile is available (/!\\ actually not yet implemented):

- `slurm`, to add if your system has a slurm executor (local by default) 

The use of the `slurm` profile  will give a command like this one: 
```
nextflow run main.nf -profile docker,slurm <rest of paramaters>
```
## Test the workflow

Test data are included in the AliNe repository in the `test` folder.

Test with short single reads:
```
nextflow run -profile docker,test_illumina_single.config aline.nf
```

Test with short paired reads:
```
nextflow run -profile docker,test_illumina_paired.config aline.nf
```

Test with ont reads:
```
nextflow run -profile docker,test_ont.config aline.nf
```

Test with pacbio reads:
```
nextflow run -profile docker,test_pacbio.config aline.nf
```

On success you should get a message looking like this:
```
  AliNe Pipeline execution summary
    --------------------------------------
    Completed at : 2024-03-07T21:40:23.180547+01:00
    UUID         : e2a131e3-3652-4c90-b3ad-78f758c06070
    Duration     : 8.4s
    Success      : true
    Exit Status  : 0
    Error report : -
```

## Parameters

## Uninstall

You can simply remove the `AliNe` directory from your computer, and remove the nextflow conda environment:
```
conda remove -n nextflow
```