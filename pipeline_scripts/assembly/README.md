# Metagenome assembly pipeline scripts
This directory contains pipeline scripts for the genome assembly.

## assembly_pipeline.sh
The main routine of the pipeline for the metagenome assembly. The script requires a file which contains the sample IDs and the ftp URLs in the tabular format (see the example file example_run_file). This script only execute the sample specified in the command argument. The user must edit the script to specify the path to the SPAdes and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./assembly_pipeline.sh example_run_file DRR127476  out_dir
```

## example_run_file
An example of run file. This kind of file can be created from RUN Selector from Sequence Read Archive.

## preprocess_pairend.sh
A script internally used by the assembly_pipeline.sh script. The user must edit this script to specify the paths to the BBtools and the human reference. Again, the location is commented with `#TODO`.
