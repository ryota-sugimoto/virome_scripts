# Metagenome assembly pipeline scripts
This directory contains pipeline scripts for the genome assembly and the spacer collection. To execute the pipeline, user must edit the following scripts to specify the paths to the programs; assembly_pipeline.sh, preprocess_pairend.sh, collect_spacers.sh, and crispr_detect.sh. See the comments in the README files for the details.

## assembly_pipeline.sh
The main routine of the protocol. The script requires a file which contains the sample IDs and the ftp URLs for the FASTQ files in tabular format (see example_run_file for example). This script only execute the sample specified in the command argument. User must edit the script to specify the path to the SPAdes and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./assembly_pipeline.sh example_run_file DRR127476 DRR127476_out
```

## example_run_file
An example of sample file. This kind of file can be created from RUN Selector from Sequence Read Archive.

## preprocess_pairend.sh
A script internally used by the assembly_pipeline.sh script for the read preprocessing. User must edit this script to specify the paths to the BBtools and the human reference. Again, the location is commented with `#TODO`.

## collect_spacers.sh
A scrirpt internally used by the assembly_pipeline.sh script for the spacer collection. This script assumes that the input directory is generated from the assembly pipeline. The user must edit the script to specify the path to the BBtools and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./collect_spacers.sh assembly_pipeline_output_directory
```

## crispr_detect.sh
A wrapper script for the CRISPRDetect program. The user must edit the script to specify the paths to the CRISPRDetect program and a tmp directory. Again. the location is commented with `#TODO`.
