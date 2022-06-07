## extract_spacers.py
Extract spacer sequences from the masked FASTQ files.
The direct repeat parts of the reads must be masked with "R".
The program does not output short (<20) or long (>50) sequences in default. These values can be changed with -s and -l options.

### Example
```
./extract_spacers.py masked.fastq > spacers.fasta
./extract_spacers.py -s 10 -l 40 -n 'sample_1' masked.fastq > spacers.fasta
```
