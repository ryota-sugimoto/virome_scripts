# Utility scripts

This directory contains utility scripts for various purposes.

## circular_contigs.py

Extract terminally redundant (TR) sequences from the input FASTA file.
The script takes a FASTA file as input and outputs only TR sequences to the stdout.
Example: ./circular_contigs input.fasta > tr.fasta

## concatemer.py

Concatenate seuqneces from the input FASTA file in tandem and output to the stdout.

Example: ./concatemer -n 2 input.fasta > concate.fasta #Concatenate two times

## filter_fasta_by_id.py


