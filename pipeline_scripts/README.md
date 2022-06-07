# Pipeline scripts for metagenome assembly and spacer collection
This directory contains pipeline scripts for the genome assembly and the spacer collection. To execute the pipeline, user must edit the following scripts to specify the paths to the programs: assembly_pipeline.sh, preprocess_pairend.sh, collect_spacers.sh, and crispr_detect.sh. See the comments in the README file for the details.

## assembly_pipeline.sh
The main routine of the protocol. The script requires a file which contains the sample IDs and the ftp URLs for the FASTQ files in tabular format (see example_run_file for example). This script processes only the sample specified in the command argument. User must edit the script to specify the path to the SPAdes and the prodigal programs, memory limit, and the number of threads. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./assembly_pipeline.sh example_run_file DRR127476 out_dir
```

## example_run_file
An example of sample file. This kind of file can be created from Run Selector from Sequence Read Archive.

## preprocess_pairend.sh
A script internally used by the assembly_pipeline.sh script for the read preprocessing. User must edit this script to specify the paths to the BBtools and the human reference, memory limit. and the number of threads. Again, the location is commented with `#TODO`.

## collect_spacers.sh
A scrirpt internally used by the assembly_pipeline.sh script for the spacer collection. This script assumes that the input directory is generated from the assembly pipeline. The user must edit the script to specify the path to the BBtools and the prodigal programs. The location that requres editing is commented with `#TODO` in the script.

### Example
```
./collect_spacers.sh assembly_pipeline_output_directory
```

## crispr_detect.sh
A wrapper script for the CRISPRDetect program. The user must edit the script to specify the paths to the CRISPRDetect program and a tmp directory. Again. the location is commented with `#TODO`.

## extract_spacers.py
Extract spacer sequences from the masked FASTQ files.
The direct repeat parts of the reads must be masked with "R".
The program does not output short (<20) or long (>50) sequences in default. These values can be changed with -s and -l options.

### Example
```
./extract_spacers.py masked.fastq > spacers.fasta
./extract_spacers.py -s 10 -l 40 -n 'sample_1' masked.fastq > spacers.fasta
