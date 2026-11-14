# Figures of the sRNA-IMP paper

Scripts that draw the paper's figures from sRNA-IMP outputs. Environment: `envs/paper.yml`
(`micromamba create -n srna_imp_paper -f envs/paper.yml`); `fig3_tarallo.py` also runs in the pipeline environment.
Every script writes `<OUTPREFIX>.png` and `.pdf` (170 mm wide); the data figures also write `<OUTPREFIX>.tsv` with the
plotted numbers and p-values.

| figure | script | input |
| --- | --- | --- |
| Fig. 1, Supplementary Fig. S1 | `fig1_workflow.py OUT` / `fig1_workflow.py OUT --supp` | none (schematic) |
| Fig. 2 simulation benchmark | `fig2_benchmark.py OUT --summary final=GLOB --summary smetaseq=GLOB` | `srna_sim_eval.py` summaries ([benchmark](../docs/benchmark.md)) |
| Fig. 3 Tarallo et al. 2019 | `fig3_tarallo.py OUTDIR OUT [--compare-outdir OUTDIR2]` | pipeline output folder after `--cohort` ([tutorial](../docs/tutorial_tarallo.md#8-fig-3-of-the-paper)) |
| Supplementary Fig. S2 CRISPR locus | `figS2_crispr_locus.py OUTDIR OUT --samples … [--labels ID=name,…]` | pipeline output folder after the `crispr` module |
| Fig. 4, Fig. 5 Expobiome HC | `fig45_hc_tables.py OUTDIR TABLES --samples …`, then `fig45_hc_figures.py TABLES OUT` | pipeline output folder after `--compare` (with stranded metaT) |

## Fig. 2

Simulated libraries named `SIM_h<host per mille>_s<seed>` (e.g. `SIM_h020_s1`: 2 % host, seed 1), each run with
sRNA-IMP (default config) and the sMETASeq re-implementation and scored with `srna_sim_eval.py`:

```bash
python3 paper/fig2_benchmark.py fig2_benchmark \
    --summary final='bench/sim/*_s[0-9].eval.summary.tsv' \
    --summary smetaseq='bench/sim/*.smetaseq.eval.summary.tsv'
```

## Fig. 4 and 5

```bash
python3 paper/fig45_hc_tables.py /path/to/expobiome/sRNAs hc_tables --type sRNA-Int --samples BHC01,BHC02,BHC03,BHC04,BHC05
python3 paper/fig45_hc_figures.py hc_tables expobiome_hc        # expobiome_hc.fig4.*, expobiome_hc.fig5.*
```

`fig45_hc_tables.py` reads the per-sample outputs (read fate, MAG placement, genus counts, tRNA fragments, CRISPR,
novel loci, the comparison tables) and the metaT strand check
(`Analysis/comparison/metaT_strand/<ID>.<type>.metaT_strand.tsv`, `COMPARE_MT_STRAND`).
