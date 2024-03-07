# AliNe (Alignment in Nextflow)
=========================================  
<img align="right" src="img/IRD.png" width="200" height="66" /> <img align="right" src="img/MIVEGEC.png" width="100" height="66" />

AliNe is a pipeline written in nextflow that aims to efficienlty align reads against a reference genome using the tools of your choice.

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
You can choose to run one or several aligner in parallel.

Here is the list of implemented aligners:

| Tool	| Single End (short reads) | Paired end (short reads) | Pacbio | ONT |
| --- | --- | --- |  --- | --- |
| bowtie2 | x | | | |
| bwaaln | x | | | |
| bwamem | x | | | |
| bwasw | x | | | |
| graphmap2 | x | | | |
| hisat2 | x | | | |
| minimap2 | x | | | |
| nucmer | x | | | |
| star | x | | | |

This work is ongoing... adaptation to run with different type of reads is ongoing.

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
cd baargin
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

A typical command to run a test on single end data will look like that:

```
nextflow run -profile docker aline.nf --aligner hisat2,graphmap2,bwamem,nucmer --genome test/hpv16.fa --reads test/U2OS_A1_R1_sub100000.fastq --single_end true --reads_extension .fastq
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