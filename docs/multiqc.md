# MultiQC

The pipeline generates a per-sample MultiQC report when `settings.run_multiqc` is enabled.

The report is written under:
- `results/<TYPE>/<sample>/multiqc/`

The pipeline uses `scripts/srna_imp_multiqc_config.yaml` as the MultiQC configuration source.
