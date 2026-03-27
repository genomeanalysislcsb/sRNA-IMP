# Quickstart

## Single Sample

```bash
micromamba run -p /scratch/users/pmay/sRNA/sRNA-IMP/.micromamba/envs/srna-imp python /scratch/users/pmay/sRNA/sRNA-IMP/scripts/srna_imp.py   --config /scratch/users/pmay/sRNA/sRNA-IMP/config/BHC01.yaml   --workflow all
```

## Multi-sample YAML

```bash
micromamba run -p /scratch/users/pmay/sRNA/sRNA-IMP/.micromamba/envs/srna-imp python /scratch/users/pmay/sRNA/sRNA-IMP/scripts/srna_imp.py   --config /scratch/users/pmay/sRNA/sRNA-IMP/config/BHC01_BHC02.multisample.yaml   --workflow all
```

## TSV-driven Multi-sample

```bash
micromamba run -p /scratch/users/pmay/sRNA/sRNA-IMP/.micromamba/envs/srna-imp python /scratch/users/pmay/sRNA/sRNA-IMP/scripts/srna_imp.py   --config /scratch/users/pmay/sRNA/sRNA-IMP/config/BHC01_BHC02.tsv_config.yaml   --workflow all
```

## Targeted Re-runs

```bash
python scripts/srna_imp.py --config config/BHC01.yaml --workflow host_specific
python scripts/srna_imp.py --config config/BHC01_BHC02.tsv_config.yaml --workflow all --sample BHC02
```
