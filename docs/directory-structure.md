# Directory Structure

```text
results/<TYPE>/<sample>/
├── Analysis/
│   ├── annotation/
│   ├── taxonomy/kraken/
│   └── <TYPE>/
│       ├── human/
│       ├── nonhuman/
│       └── novel/
├── Assembly/
├── Preprocessing/<TYPE>/
├── Stats/<TYPE>/
├── logs_cluster_jobs/<TYPE>/
└── multiqc/
```

Important files include:
- `sample.json`
- `Preprocessing/<TYPE>/<sample>.preprocessing.summary.tsv`
- `Analysis/<TYPE>/human/<sample>.human.summary.tsv`
- `Analysis/<TYPE>/nonhuman/<sample>.nonhuman.summary.tsv`
- `Analysis/<TYPE>/novel/<sample>.novel.summary.tsv`
- `Stats/<TYPE>/<sample>.pipeline.summary.tsv`
