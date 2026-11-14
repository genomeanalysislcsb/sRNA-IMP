#!/usr/bin/env python3
"""fig45_hc_tables.py — gather the Expobiome HC demonstration numbers (sRNA-Int, BHC01-05) from sRNA-IMP outputs.
Usage: fig45_hc_tables.py RUNDIR OUTDIR [--type sRNA-Int] [--samples BHC01,...]"""
import argparse, glob, os, re, sys
import pandas as pd, numpy as np
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("rundir"); ap.add_argument("outdir")
ap.add_argument("--type", default="sRNA-Int")
ap.add_argument("--samples", default="")
ap.add_argument("--mt-strand", default="", help="only without the pipeline's Analysis/comparison/metaT_strand tables: folder with featureCounts outputs <sample>.mt.s1.txt / .s2.txt "
                "(stranded metaT counted per CDS; sample name case-insensitive; s2 = gene strand for TruSeq Stranded)")
ap.add_argument("--mt-own", default="s2", help="which file holds the gene-strand metaT counts [s2]")
a = ap.parse_args()
T = a.type
S = a.samples.split(",") if a.samples else sorted(os.path.basename(p) for p in glob.glob(f"{a.rundir}/*") if os.path.isdir(f"{p}/Analysis/{T}"))
os.makedirs(a.outdir, exist_ok=True)
R = lambda s, *p: os.path.join(a.rundir, s, *p)
rd = lambda p, **k: pd.read_csv(p, sep="\t", **k)

def genus(c):
    m = re.search(r"g__([^;]*)", str(c)); return m.group(1) if m and m.group(1) else None

fate, place, mags, gen, trf, layers, rho = [], [], [], [], [], [], []
for s in S:
    st = rd(R(s, "Stats", T, f"{s}.{T}.read_stats.tsv")).set_index("step")["reads"]
    cat = rd(R(s, "Analysis", T, "assembly", f"{s}.category_counts.tsv")).set_index("category")["reads"]
    ms = rd(R(s, "Analysis", T, "assembly", f"{s}.mapping_summary.tsv")).set_index("category")["value"]
    trimmed, host, nonhost, mapped = st["trimmed"], st["human_total"], st["nonhuman"], float(ms["mapped_reads"])
    g = lambda *k: sum(cat.get(x, 0) for x in k)
    anti = sum(v for k, v in cat.items() if k.startswith("antisense_"))
    other = mapped - g("rRNA", "tRNA", "CDS", "intergenic") - anti
    f = {"sample": s, "trimmed": trimmed, "host": host, "rRNA": g("rRNA"), "tRNA": g("tRNA"),
         "CDS (sense)": g("CDS"), "antisense": anti, "other annotated": other, "intergenic": g("intergenic"),
         "unmapped": nonhost - mapped}
    fate.append(f)
    # placement: share in MAGs per omics (contig basis)
    al = rd(R(s, "Analysis", "comparison", f"{s}.all_sources.bins.tsv"))
    tot = al[[c for c in al if c.endswith("_count")]].sum()
    ub = al[al.bin == "unbinned"].iloc[0]
    hq = al[al.quality.isin(["HQ", "MQ"])]
    for o in [T, "metaT", "metaG"]:
        place.append({"sample": s, "omics": o, "in_MAGs": 1 - ub[f"{o}_count"] / tot[f"{o}_count"],
                      "in_HQ_MQ_MAGs": hq[f"{o}_count"].sum() / tot[f"{o}_count"], "n_MAGs": int((al.bin != "unbinned").sum())})
    # per MAG, rRNA contigs removed
    b = rd(R(s, "Analysis", "comparison", f"{s}.all_sources.bins_no_rrna.tsv"))
    b = b[b.bin != "unbinned"].copy(); b["sample"] = s; b["genus"] = b.classification.map(genus)
    b["log2_sRNA_vs_metaT"] = np.log2((b[f"{T}_CPM"] + 1) / (b["metaT_CPM"] + 1))
    b["log2_sRNA_vs_metaG"] = np.log2((b[f"{T}_CPM"] + 1) / (b["metaG_CPM"] + 1))
    mags.append(b)
    for o in ["metaT", "metaG"]:
        r = spearmanr(b[f"{T}_CPM"], b[f"{o}_CPM"]).correlation
        rho.append({"sample": s, "vs": o, "level": "MAG (no rRNA contigs)", "spearman": r, "n": len(b)})
    cs = rd(R(s, "Analysis", "comparison", f"{s}.compare.summary.tsv"))
    for o in ["metaT", "metaG"]:
        x = cs[(cs.pair == f"{T}_vs_{o}") & (cs.metric.isin(["spearman_cpm", "spearman_cpm_within_type"]))]
        for _, r in x.iterrows():
            rho.append({"sample": s, "vs": o, "level": f"{r.level}:{r.feature_type}", "spearman": float(r.value), "n": np.nan})
    # tRNA fragment classes
    fc = rd(R(s, "Analysis", T, "fragments", f"{s}.fragments.classes.tsv"))
    for _, r in fc[fc.type == "tRNA"].iterrows():
        trf.append({"sample": s, "class": r["class"], "fraction": r.fraction_of_type, "reads": r.weighted_reads})
    # layers
    cr = rd(R(s, "Analysis", T, "crispr", f"{s}.crispr_arrays.tsv"))
    rt = rd(R(s, "Analysis", T, "nonhuman", "rfam", f"{s}.rfam_type.counts.tsv")).set_index("rfam_type")["count"]
    nov = rd(R(s, "Analysis", T, "novel", f"{s}.assembly.summary.tsv")).set_index("category")["reads"]
    mi = rd(R(s, "Analysis", T, "human", "mirna", f"{s}.human.mirna.counts.tsv"))
    hs = rd(R(s, "Analysis", T, "human", f"{s}.human.summary.tsv")).set_index("category")["reads"]
    layers.append({"sample": s,
        "host reads (% of trimmed)": 100 * host / trimmed,
        "host miRNA reads": hs.get("human_miRNA", 0), "host miRNAs >=10 reads": int((mi["count"] >= 10).sum()),
        "antisense reads (% of mapped)": 100 * anti / mapped,
        "antisense CDS reads": cat.get("antisense_CDS", 0),
        "CRISPR arrays": len(cr), "CRISPR arrays >=10 reads": int((cr.reads >= 10).sum()),
        "spacers with reads": int(cr.spacers_with_reads.sum()), "spacers in arrays": int(cr.spacers.sum()),
        "crRNA reads (spacer-containing)": cr.reads_with_spacer.sum(),
        "riboswitch reads (Rfam)": rt.get("Cis-reg;riboswitch", 0), "leader reads (Rfam)": rt.get("Cis-reg;leader", 0),
        "novel intergenic loci": int(nov.get("novel_intergenic_candidates", 0)),
        "novel antisense loci": int(nov.get("novel_antisense_candidates", 0))})

# antisense on CDS vs metaT of the gene (needs assembly/<s>.feature_counts.tsv and the sRNA vs metaT features table);
# with --mt-strand also metaT reads on the gene strand and on the opposite (= antisense sRNA) strand
def fcounts(path, name):
    f = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    return f.iloc[:, [0, -1]].set_axis(["feature_id", name], axis=1)

def strand_class(o, p):
    return np.where((o >= 10) & (p >= 10), "both strands", np.where(p >= 10, "sRNA strand only",
                    np.where(o >= 10, "gene strand only", "<10 metaT reads")))

anti = []
for s in S:
    fc = R(s, "Analysis", T, "assembly", f"{s}.feature_counts.tsv")
    mt = R(s, "Analysis", "comparison", f"{T}_vs_metaT", f"{s}.{T}_vs_metaT.features.tsv")
    if not (os.path.exists(fc) and os.path.exists(mt)):
        continue
    c = rd(fc, low_memory=False); c = c[c.type == "CDS"].copy(); c["len"] = c.end - c.start
    m = rd(mt, usecols=["feature", "metaT_count"])
    c = c.merge(m, left_on="feature_id", right_on="feature", how="left"); c["metaT_count"] = c.metaT_count.fillna(0)
    c["sample"] = s
    c["metaT_class"] = pd.cut(c.metaT_count, [-1, 0, 9.999, np.inf], labels=["no metaT reads", "1-9 metaT reads", ">=10 metaT reads"])
    c["hypothetical"] = c["product"].fillna("").str.contains("hypothetical")
    cols = ["sample", "feature_id", "contig", "start", "end", "strand", "len", "product", "hypothetical", "sense_reads",
            "antisense_reads", "metaT_count", "metaT_class"]
    pm = R(s, "Analysis", "comparison", "metaT_strand", f"{s}.{T}.metaT_strand.tsv")
    if os.path.exists(pm):   # pipeline output (--compare, COMPARE_MT_STRAND)
        x = rd(pm)[["feature_id", "metaT_gene_strand", "metaT_opposite_strand", "metaT_strand_class"]]
        c = c.merge(x, on="feature_id", how="left")
        c["metaT_strand_class"] = c.metaT_strand_class.fillna("low metaT")
        c[["metaT_gene_strand", "metaT_opposite_strand"]] = c[["metaT_gene_strand", "metaT_opposite_strand"]].fillna(0)
        cols += ["metaT_gene_strand", "metaT_opposite_strand", "metaT_strand_class"]
    elif a.mt_strand:
        opp = "s1" if a.mt_own == "s2" else "s2"
        fs = {k: [f for f in glob.glob(f"{a.mt_strand}/*.mt.{k}.txt") if os.path.basename(f).lower().startswith(s.lower() + ".")]
              for k in (a.mt_own, opp)}
        if all(fs.values()):
            c = c.merge(fcounts(fs[a.mt_own][0], "metaT_gene_strand"), on="feature_id", how="left")
            c = c.merge(fcounts(fs[opp][0], "metaT_opposite_strand"), on="feature_id", how="left")
            c[["metaT_gene_strand", "metaT_opposite_strand"]] = c[["metaT_gene_strand", "metaT_opposite_strand"]].fillna(0)
            c["metaT_strand_class"] = strand_class(c.metaT_gene_strand, c.metaT_opposite_strand)
            cols += ["metaT_gene_strand", "metaT_opposite_strand", "metaT_strand_class"]
            tot = c.metaT_gene_strand.sum() + c.metaT_opposite_strand.sum()
            print(f"{s}: metaT on gene strand {c.metaT_gene_strand.sum() / tot:.3f} of CDS reads", file=sys.stderr)
    anti.append(c[c.antisense_reads >= 10][cols])
anti = pd.concat(anti) if anti else pd.DataFrame()
anti.to_csv(f"{a.outdir}/hc_antisense_cds.tsv", sep="\t", index=False)
fate = pd.DataFrame(fate); place = pd.DataFrame(place); mags = pd.concat(mags); trf = pd.DataFrame(trf)
layers = pd.DataFrame(layers); rho = pd.DataFrame(rho)
# genus recurrence: sum CPM per genus per sample
gg = mags.dropna(subset=["genus"]).groupby(["genus", "sample"])[[f"{T}_CPM", "metaT_CPM", "metaG_CPM"]].sum().reset_index()
gg["log2_sRNA_vs_metaT"] = np.log2((gg[f"{T}_CPM"] + 1) / (gg["metaT_CPM"] + 1))
for n, d in [("fate", fate), ("placement", place), ("mags", mags), ("genus", gg), ("trna_fragments", trf),
             ("layers", layers), ("spearman", rho)]:
    d.to_csv(f"{a.outdir}/hc_{n}.tsv", sep="\t", index=False)
print(fate.to_string()); print(place.to_string()); print(layers.T.to_string())
print(rho[rho.level.str.startswith(("MAG", "feature:all", "contig", "bin:all"))].to_string())
