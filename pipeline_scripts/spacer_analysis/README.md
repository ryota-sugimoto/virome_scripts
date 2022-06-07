# Spacer pipeline scripts
This directory contains pipeline scripts for collecting CRISPR spacers.

## collect_spacers.sh
The main routine of the pipeline for collecting CRISPR spacers. This script assumes that the input directory is generated from the assembly pipeline. The user must edit the script to specify the path to the BBtools and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./collect_spacers.sh assembly_pipeline_output_directory
```

## crispr_detect.sh
A wrapper script for the CRISPRDetect program. The user must edit the script to specify the paths to the CRISPRDetect program and a tmp directory. Again. the location is commented with `#TODO`.
