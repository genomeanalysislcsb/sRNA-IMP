# Inputs and Outputs

## Required Inputs

Required per sample:
- one raw, non-filtered single-end smallRNA FASTQ file
- host reference FASTA files or prebuilt host Bowtie indexes

## Optional Inputs

Optional per sample:
- assembly FASTA
- assembly GFF

## Main Output Groups

For each sample, the pipeline writes:
- preprocessing FASTQs and QC files
- host analysis outputs
- non-host analysis outputs
- assembly-based outputs
- novel ncRNA outputs
- summary TSVs
- MultiQC report
