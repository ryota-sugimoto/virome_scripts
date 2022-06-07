# Utility scripts

This directory contains utility scripts for the various purposes.

## circular_contigs.py

Extract terminally redundant (TR) sequences from the input FASTA file.
The script takes a FASTA file as input and outputs only TR sequences to the stdout.

### Example
```
./circular_contigs.py input.fasta > tr.fasta
```

## concatemer.py

Concatenate seuqneces from the input FASTA file in tandem and output to the stdout. The number of repetition can be changed with the `-n` option.

### Example
```
./concatemer.py -n 2 input.fasta > concate.fasta #Concatenate two times
```

## filter_fasta_by_id.py

Filter out the FASTA records with specific IDs from the input FASTA file. The IDs can be specifed by the command or a file with the `-f` option. The ID file should be listing the FASTA records in each line.

### Example
```
./filter_fasta_by_id.py 'record_id_1' input.fasta > filtered.fasta
./filter_fasta_by_id.py -f id_file input.fasta > filtered.fasta
```

## gc_contents.py

Calculate GC contents of the sequence. Each line of the stdout corresponds to a record in the input FASTA file, which is seperated by tabular. The first column is the record ID and the second column is the GC%.

### Example
```
./gc_contents.py input.fasta > gc_contents
```

## inverted_repeat_contigs.py
Extract terminally inverted repeat (TIR) sequences from the input FASTA file. The script takes a FASTA file as input and outputs only TIR sequences to the stdout.

### Example
```
./inverted_repeat_contigs.py input.fasta > tir.fasta
```

## lcc.py
Calculate local complexity over the FASTA sequence. Each line of the stdout corresponds to a record in the input FASTA file, which is seperated by tabular. The first column is the record ID and the second column is the local complexity.

### Example
```
./lcc.py input.fasta > lcc
```

## parse_hhr.py
Parse input HHR file to the tabular format.

### Example
```
./parse_hhr.py input.hhr > parsed_hhr
```

## phase_circular_genome.py
Phase circular genome based on the self-aligned blast result. The script requires two files; fasta and the blast output file genereted with the `-outfmt 6` option.

### Example
```
./phase_circular_genome.py input.fasta blastn_output > phased.fasta
```

## remove_gapped.py
Remove records that contain gapped sequence from the input FASTA file.

### Example
```
./remove_gapped.py input.fasta > without_gapped.fasta
```

## rename_fasta.py
Rename IDs in the input FASTA file according to the tabular formatted text file. The tabular text file contains two columns; the first column is the original ID and the second column is the new ID.

### Example
```
./rename_fasta.py rename_list input.fasta > renamed.fasta
```

## reverse_complement.py
Output revese complement of the all sequnces in the FASTA file.

### Example
```
./reverse_complement.py input.fasta > reverse_complemented.fasta
```

## shuffle_fasta_record_order.py
Shuffle the order of records in the FASTA file.

### Example
```
./shuffle_fasta_record_order.py input.fasta > shuffled_id.fasta
```

## shuffle_sequence.py
Shuffle sequences in the FASTA file randomly.

### Example
```
./shuffle_sequence.py input.fasta > shuffled_sequence.fasta
```

## split_fasta.py
Split records in the input FASTA file to the seperate files. The number of records in a split file can be change with the `-s` option.

### Example
```
./split_fasta.py input.fasta #contain a record in each file
./split_fasta.py -s 10 input.fasta #contain at most 10 records in each file
```

## stockholm2fasta.py
Convert stockhole format MSA file to FASTA formated file.

### Example
```
./stockholm2fasta.py stockholm.sto > msa.fasta
```
