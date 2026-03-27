# Human Workflow

The host-specific branch performs sequential classification of host reads against:
- mature miRNA
- tRNA
- rRNA
- other ncRNA

The native pipeline now mirrors the updated legacy `human_final_pipeline.sh` behavior more closely:
- intermediate remaining-read FASTQs are retained in the human subtree
- count files use the `.human.*.counts.tsv` naming pattern
- summary generation is delegated to `summarize_human_subtree.py`

Expected summary output:
- `<sample>.human.summary.tsv`
