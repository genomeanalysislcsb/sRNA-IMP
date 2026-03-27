# Configuration

## Supported Config Styles

`sRNA-IMP` supports three configuration patterns:
- single-sample YAML
- multi-sample YAML with a `samples:` list
- TSV-driven multi-sample mode with `samples_tsv:` in YAML

## Shared Sections

Typical shared sections are:
- `project`
- `settings`
- `tools`
- `workflows`

## Sample-specific Inputs

Per sample, the pipeline needs:
- `name`
- `type`
- `fastq`
- optional `assembly_fasta`
- optional `assembly_gff`

## Important Workflow Settings

### preprocessing
- `min_length`
- `max_length`
- `run_fastqc`
- `host_reference_fastas`
- optional `host_filter_indexes`

### host_specific
- `mirna_index`
- `trna_index`
- `rrna_index`
- `other_ncrna_index`

### nonhost_specific
- `kraken_db`
- `bracken_db`
- `trna_index`
- `mirna_euk_index`
- `rfam_class_index`
- `rfamseq_tsv`

### novel_ncrna
- `min_length`
- `max_length`
- `min_reads`
- `merge_distance`
- `mismatches`
- optional `rfam_cm`
