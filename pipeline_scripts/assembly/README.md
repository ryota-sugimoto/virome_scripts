# Metagenome assembly pipeline scripts
This directory contains pipeline scripts for the genome assembly and spacer collection. To execute the pipeline, user must edit the following scripts to specify the paths to the programs; assembly_pipeline.sh, preprocess_pairend.sh, collect_spacers.sh, and crispr_detect.sh. See the comments in the README files for the details.

## assembly_pipeline.sh
The main routine of the protocol. The script requires a file which contains the sample IDs and the ftp URLs for the FASTQ files in tabular format (see the example file example_run_file). This script only execute the sample specified in the command argument. User must edit the script to specify the path to the SPAdes and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./assembly_pipeline.sh example_run_file DRR127476  out_dir
```

## example_run_file
An example of run file. This kind of file can be created from RUN Selector from Sequence Read Archive.

## preprocess_pairend.sh
A script internally used by the assembly_pipeline.sh script. User must edit this script to specify the paths to the BBtools and the human reference. Again, the location is commented with `#TODO`.
