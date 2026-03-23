# Directory Structure

```text
results/
└── smallRNA/
    └── <TYPE>/
        └── <sample>/
            ├── 01_preprocessing/
            │   ├── <sample>.host.fastq.gz
            │   ├── <sample>.nonhost.fastq.gz
            │   ├── <sample>.read_stats.tsv
            │   └── qc/
            ├── 02_host_specific/
            ├── 03_nonhost_specific/
            ├── reports/
            │   ├── logs/
            │   └── <sample>.summary.tsv
            └── sample.json
```
