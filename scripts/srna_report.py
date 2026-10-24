#!/usr/bin/env python3
"""
srna_report.py — per-sample PDF report of the sRNA pipeline (srna_pipeline.sh --report)

Reads only the summary tables the pipeline already writes under <outdir>/<sample>/ and puts the key
numbers and plots of all sRNA types of the sample (e.g. sRNA-Ext and sRNA-Int) side by side:

  1  overview        inputs, read fate from raw reads to classes, key numbers
  2  QC              trimming, adapter dimers, contaminants, overrepresented sequences
  3  human           human share, assembly rescue, length distribution, classes, top miRNAs
  4  non-human       classes, Bracken genus / family, Rfam types, tRNA amino acids
  5  assembly        reads per feature category, tRNA / rRNA fragment classes and lengths
     integrity       RNA integrity indicators (fragment end definition, rRNA spread, read lengths)
  6  novel           candidate loci, read support, Rfam families of candidates, top candidates
  7  bins / omics    top bins, Spearman correlations (all pairs, bases), scatter plots
                     (only with --compare results)

Missing inputs (modules not run) are skipped with a note. Optionally writes MultiQC custom-content
tables (*_mqc.tsv) so that MultiQC can later summarise all samples.

Usage
  srna_report.py --outdir <pipeline outdir> --sample BXXX1 [--types sRNA-Ext,sRNA-Int]
                 [--pdf <file>] [--html <file>] [--mqc-dir <dir>]

--html writes the same pages as one self-contained HTML file (PNG images embedded, page menu) for
quick checks in a browser.
"""
import argparse
import urllib.parse
import base64
import collections
import gzip
import textwrap
import csv
import html as htmlmod
import io
import datetime
import glob
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

PAGE = (11.69, 8.27)                       # A4 landscape
TYPE_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
SOURCE_COLORS = {"sRNA-Ext": "#0072B2", "sRNA-Int": "#D55E00", "metaT": "#009E73", "metaG": "#CC79A7"}
CLASS_COLORS = {"rRNA": "#999999", "tRNA": "#0072B2", "Rfam": "#009E73", "miRNA-like": "#CC79A7",
                "unknown": "#E69F00", "human": "#D55E00", "<20 nt": "#F0E442", "too short / dimer": "#cccccc"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False})


# ----------------------------------------------------------------------------- helpers
def read_tsv(path):
    """list of dicts; [] if the file is missing/empty; '#' lines are skipped"""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def kv(path, key="category", val="count"):
    """two-column style table -> {key: float}; also handles 3-column sample/category/count"""
    out = {}
    for r in read_tsv(path):
        k = r.get(key) or r.get("metric") or r.get("step") or list(r.values())[1 if len(r) > 2 else 0]
        v = r.get(val) or r.get("value") or r.get("reads") or list(r.values())[-1]
        try:
            out[k] = float(v)
        except (TypeError, ValueError):
            out[k] = v
    return out


def num(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def fmt(x, nd=0):
    if x is None or x == "":
        return "–"
    if isinstance(x, str):
        return x
    if nd == 0:
        return "{:,.0f}".format(x)
    return "{:,.{}f}".format(x, nd)


def pct(a, b, nd=1):
    return "–" if not b else "{:.{}f}%".format(100.0 * a / b, nd)


def short(s, n=38):
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "…"


def tax_tail(classification):
    """GTDB string -> deepest named rank"""
    parts = [p for p in str(classification).split(";") if len(p) > 3]
    return parts[-1][3:] if parts else "unclassified"


def table(ax, rows, col_labels=None, col_widths=None, fontsize=7, title=None, header_color="#e8eef5", left=(0,), rowscale=1.25):
    ax.axis("off")
    if title:
        ax.set_title(title, loc="left", fontweight="bold")
    if not rows:
        ax.text(0, 0.5, "no data", fontsize=8, color="#888")
        return
    t = ax.table(cellText=[[str(c) for c in r] for r in rows], colLabels=col_labels, loc="upper left",
                 cellLoc="right", colLoc="right", colWidths=col_widths)
    t.auto_set_font_size(False)
    t.set_fontsize(fontsize)
    t.scale(1, rowscale)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if c in left:
            cell.set_text_props(ha="left")
            cell._loc = "left"
        if r == 0 and col_labels:
            cell.set_facecolor(header_color)
            cell.set_text_props(fontweight="bold")


def note(ax, text):
    ax.axis("off")
    ax.text(0.0, 0.5, text, fontsize=8, color="#888", va="center")


EXPLAIN = {
    "overview":
        "sRNA-Int is the cell fraction (bacterial, human and food-derived cells), sRNA-Ext the cell-free fraction "
        "(RNA released by secretion, membrane vesicles or lysis). In stool most small RNAs are processed or degraded "
        "bacterial rRNA and tRNA; human RNA is a small minority. The read fate shows where reads are removed and "
        "how the rest is assigned; shares are relative (compositional), not absolute amounts.",
    "qc":
        "Reads are adapter- and quality-trimmed and kept at 18–MAXLEN nt. A large 'too short' share means heavily "
        "degraded RNA or adapter/primer dimers (typical of low-input cell-free RNA). The most frequent sequences in gut "
        "sRNA libraries are usually tRNA halves (30–35 nt) and 5′ ends of 16S/23S rRNA; reads of the full read length "
        "carry no adapter, i.e. the RNA fragment was longer than the read.",
    "human":
        "Human reads come mainly from cells shed from the gut lining and, if present, immune cells or blood. Marker "
        "miRNAs indicate the origin: intestinal/epithelial (miR-192/194/215, miR-200 family, miR-141, miR-375) = shed "
        "colonocytes; blood (miR-451a, miR-486, miR-144, miR-16) = bleeding; immune (miR-142, miR-150, miR-223, "
        "miR-155) = immune cells / inflammation; muscle (miR-1, miR-133, miR-206) = most likely dietary meat. Intact "
        "miRNAs are protected in cells and Argonaute complexes, whereas cell-free human RNA is degraded quickly, so Ext "
        "usually holds fewer human reads and hardly any miRNAs.",
    "nonhuman":
        "Taxonomy of rRNA fragments (Kraken2/Bracken, GTDB) shows which taxa contribute rRNA-derived small RNAs. This "
        "mixes abundance, ribosome content (growth) and RNA release, so it differs from metagenomic abundance; species "
        "calls from short fragments are uncertain – genus and family are more robust. Rfam types list the other "
        "non-coding RNAs (sRNAs, riboswitches, tmRNA, RNase P, ribozymes); tRNA amino-acid profiles reflect isoacceptor "
        "abundance plus differential cleavage and stability.",
    "species":
        "Species-level view of the rRNA taxonomy, shown with the evidence behind each call. Most rRNA fragments are "
        "20–50 nt from conserved regions, so Kraken2 usually stops at genus or family; Bracken then redistributes "
        "these reads to species using the k-mer profiles of the database, i.e. largely model-based. Species resting on "
        "few distinct minimizers (one short conserved stretch) or few directly assigned reads are flagged 'weak'. "
        "For species-level questions the genome bins (metaG, GTDB) are the more reliable reference.",
    "asm_cmc":
        "All non-human reads are mapped once to the sample's own contigs (bowtie, uncollapsed, strand kept); this BAM is "
        "used by all assembly-based analyses. Reads with several equally good positions count 1/hits per position. Each "
        "read is assigned to the annotated feature it overlaps best on its own strand (sense); only if there is none, to "
        "a feature on the opposite strand (antisense); otherwise it is intergenic. Reads on protein-coding genes are "
        "mostly mRNA degradation fragments; antisense reads point to antisense transcription or cis-encoded sRNAs.",
    "ncrna":
        "Non-coding RNA genes (ncRNA) and cis-regulatory elements (riboswitches, leaders, thermometers, ...) annotated "
        "by Bakta with Rfam models, with the reads (bowtie mapping, weighted 1/hits) on their own strand (sense) or the "
        "opposite strand (antisense). Housekeeping ncRNAs (RNase P, SRP, 6S, tmRNA) are abundant and stable; trans-"
        "acting sRNAs are regulators of bacterial gene expression. Riboswitches and leaders are parts of mRNA 5′ UTRs "
        "and usually get few reads unless premature termination releases them as short RNAs.",
    "integrity":
        "sRNA libraries are size-selected, so RIN-like 5′→3′ coverage skew does not apply. Instead, integrity is read "
        "from how defined the fragment ends are: processed RNAs (mature tRNAs, tRNA halves, crRNAs) have sharp ends, "
        "random chemical or nuclease degradation gives diffuse ends (i-tRFs, loss of the 3′-CCA, high end entropy) and "
        "an even spread of fragments along 16S/23S. tRNA halves are a stress response, not degradation. Compare samples "
        "within one fraction: cell-free sRNA-Ext is expected to be more fragmented than sRNA-Int for biological "
        "reasons. Values are relative indicators, best checked against RIN/DV200 of the same extractions.",
    "qcverdict":
        "Rule-based QC of each library (thresholds in the QC rules file, manual decisions in QC_EXCLUDE). FAIL "
        "libraries are left out of the cohort statistics; WARN libraries stay in but are flagged. The FastQC grid "
        "shows the module verdicts of raw and trimmed reads: for sRNA libraries 'per base sequence content', "
        "'sequence length distribution', 'duplication' and 'overrepresented sequences' are usually WARN/FAIL by "
        "nature (short, highly repeated molecules) and are not used as rules by default.",
    "crispr":
        "CRISPR arrays are transcribed into a long pre-crRNA that is cut into mature crRNAs (one spacer plus part of a "
        "repeat), which guide Cas proteins to matching phage or plasmid DNA. sRNA reads on an array that pile up on "
        "single spacers, come from one strand and have defined ends are processed crRNAs, i.e. an active CRISPR-Cas "
        "system; type II crRNAs typically start inside the spacer and end ~20 nt into the next repeat. Highly "
        "expressed spacers point to the phages or plasmids the taxon is currently defending against. Type II (and "
        "V-B) systems also need a tracrRNA, often encoded next to cas9: a read cluster outside the array whose "
        "sequence pairs with the crRNA repeat (anti-repeat) is marked green as tracrRNA candidate.",
    "assembly":
        "Non-human reads mapped back to the sample's own contigs. tRNA halves (cleavage in the anticodon loop, "
        "≈30–40 nt) and 5′/3′ tRNA fragments (tRFs) are produced by stress-induced ribonucleases and are major "
        "components of extracellular RNA. rRNA reads are mostly internal fragments. 'Intergenic' covers unannotated "
        "regions: partial rRNA operons missed by the annotation, bacterial sRNAs, UTRs and leaders.",
    "novel":
        "Read clusters outside annotated features (intergenic) or antisense to a protein-coding gene. Most clusters with "
        "an Rfam hit are rRNA/tRNA pieces the (partial) contig annotation missed; they are summarised here. Clusters "
        "with another Rfam family are known ncRNAs, and clusters without any hit, especially with a stable predicted "
        "structure (low MFE), are candidates for so far unknown bacterial sRNAs, UTR-derived or antisense RNAs – to be "
        "checked for conservation, promoters/terminators and expression in metaT.",
    "bins":
        "Per-bin comparison of sRNA with genome abundance (metaG) and transcription (metaT). A high sRNA–metaG "
        "correlation means sRNA output scales with how abundant a taxon is; bins far from the diagonal release more or "
        "less sRNA than their abundance predicts. Most sRNA reads lie on unbinned contigs (multi-copy rRNA/tRNA), so "
        "bins represent only part of the signal.",
}


PIPELINE_NAME = "sRNA-IMP pipeline"      # shown in page footers, HTML header and PDF metadata
PIPELINE_URL = "github.com/genomeanalysislcsb/sRNA-IMP"


def page(pdf, sample, title, subtitle="", explain=""):
    fig = plt.figure(figsize=PAGE)
    fig.text(0.03, 0.965, "%s — %s" % (sample, title), fontsize=14, fontweight="bold", va="top")
    if subtitle:
        fig.text(0.03, 0.927, subtitle, fontsize=7.5, color="#555", va="top")
    fig._explain = EXPLAIN.get(explain, explain)
    fig.text(0.97, 0.012, "%s report · %s · %s" % (PIPELINE_NAME, PIPELINE_URL, datetime.date.today().isoformat()),
             fontsize=6, color="#999", ha="right")
    return fig


def finalize(fig):
    """put the explanation box under the title and shrink all axes so that nothing (titles, tables) overlaps it"""
    text = getattr(fig, "_explain", "")
    if not text:
        return
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    # wrap so that the box stays inside the page: measure the rendered width and narrow if needed
    width = 190
    while True:
        t = fig.text(0.03, 0.905, "\n".join(textwrap.wrap(text, width)), fontsize=7.2, color="#333", va="top",
                     linespacing=1.35, bbox=dict(boxstyle="round,pad=0.4", fc="#f3f6fa", ec="#d5dde8"))
        bb = t.get_window_extent(renderer=renderer).transformed(inv)
        if bb.x1 <= 0.965 or width <= 80:
            break
        t.remove()
        width = int(width * 0.97 * 0.945 / max(bb.x1 - 0.03, 0.01)) if bb.x1 > 1.2 else width - 6
    # bottom of the box incl. its padding
    t.update_bbox_position_size(renderer)
    patch = t.get_bbox_patch()
    box_bottom = patch.get_window_extent(renderer).transformed(inv).y0 if patch else bb.y0 - 0.01
    # highest point of all axes incl. titles / tables, then scale the axes (and their gaps) into the space below
    top = 0.0
    for ax in fig.axes:
        try:
            tb = ax.get_tightbbox(renderer)
            top = max(top, tb.transformed(inv).y1 if tb is not None else ax.get_position().y1)
        except Exception:
            top = max(top, ax.get_position().y1)
    limit = box_bottom - 0.012
    if top > limit and top > 0:
        k = limit / top
        for ax in fig.axes:
            p = ax.get_position()
            ax.set_position([p.x0, p.y0 * k, p.width, p.height * k])


def grouped_barh(ax, labels, series, colors, title, xlabel="", log=False, pct_axis=False):
    """series: list of (name, [values aligned with labels])"""
    n = len(series)
    h = 0.8 / max(n, 1)
    y = list(range(len(labels)))
    for i, (name, vals) in enumerate(series):
        ax.barh([v + (i - (n - 1) / 2) * h for v in y], vals, height=h, label=name, color=colors[i % len(colors)])
    ax.set_yticks(y)
    ax.set_yticklabels([short(l, 32) for l in labels])
    ax.invert_yaxis()
    ax.set_title(title, loc="left")
    if xlabel:
        ax.set_xlabel(xlabel)
    if log:
        ax.set_xscale("log")
    if n > 1:
        ax.legend(frameon=False, loc="lower right")


def top_union(dicts, k):
    tot = {}
    for d in dicts:
        s = sum(d.values()) or 1
        for key, v in d.items():
            tot[key] = max(tot.get(key, 0), v / s)
    return [x for x, _ in sorted(tot.items(), key=lambda kv_: -kv_[1])[:k]]


def frac(d):
    s = sum(d.values()) or 1
    return {k: v / s for k, v in d.items()}


# ----------------------------------------------------------------------------- data
class TypeData:
    def __init__(self, base, sample, t):
        self.t = t
        st = os.path.join(base, "Stats", t)
        an = os.path.join(base, "Analysis", t)
        self.inputs = {r[0]: r[1] for r in self._rows(os.path.join(base, "log_cluster_jobs", t, "inputs.tsv"))}
        self.read_stats = kv(os.path.join(st, "%s.%s.read_stats.tsv" % (sample, t)), "step", "reads")
        self.qc = {r["metric"]: num(r["value"]) for r in read_tsv(os.path.join(st, "%s.%s.qc_summary.tsv" % (sample, t)))}
        self.over = read_tsv(os.path.join(st, "%s.%s.overrepresented.tsv" % (sample, t)))
        self.hlen = read_tsv(os.path.join(st, "%s.%s.human_length_dist.tsv" % (sample, t)))
        self.human = kv(os.path.join(an, "human", "%s.human.summary.tsv" % sample), "category", "count")
        self.mirna = read_tsv(os.path.join(an, "human", "mirna", "%s.human.mirna.counts.tsv" % sample))
        self.nonhuman = kv(os.path.join(an, "nonhuman", "%s.nonhuman.summary.tsv" % sample), "category", "count")
        self.bracken = {lv: read_tsv(os.path.join(an, "nonhuman", "rrna_tax", "%s.rRNA.bracken.%s.tsv" % (sample, lv)))
                        for lv in "DPCOFGS"}
        self.kraken = self._kraken(os.path.join(an, "nonhuman", "rrna_tax", "%s.rRNA.kraken.report" % sample))
        self.rfam_type = {r["rfam_type"]: num(r["count"]) for r in
                          read_tsv(os.path.join(an, "nonhuman", "rfam", "%s.rfam_type.counts.tsv" % sample))}
        self.trna_aa = {r["aa"]: num(r["count"]) for r in
                        read_tsv(os.path.join(an, "nonhuman", "trna", "%s.trna_aa.counts.tsv" % sample))
                        if r.get("aa") not in (None, "", "aa")}
        self.category = {r["category"]: num(r["reads"]) for r in
                         read_tsv(os.path.join(an, "bins", "%s.category_counts.tsv" % sample))}
        self.frag_classes = read_tsv(os.path.join(an, "fragments", "%s.fragments.classes.tsv" % sample))
        self.frag_lengths = read_tsv(os.path.join(an, "fragments", "%s.fragments.lengths.tsv" % sample))
        self.integrity = read_tsv(os.path.join(an, "fragments", "%s.integrity.tsv" % sample))
        self.integrity_profile = read_tsv(os.path.join(an, "fragments", "%s.integrity.profile.tsv" % sample))
        asm = os.path.join(an, "assembly") if os.path.isdir(os.path.join(an, "assembly")) else an
        ap_ = os.path.join(asm, sample)
        self.asm_summary = kv(ap_ + ".mapping_summary.tsv", "category", "value")
        self.asm_category = {r["category"]: num(r["reads"]) for r in read_tsv(ap_ + ".category_counts.tsv")}
        self.asm_category_nr = {r["category"]: num(r["reads_non_rRNA"]) for r in read_tsv(ap_ + ".category_counts.tsv")
                                if r.get("reads_non_rRNA") not in (None, "")}
        self.asm_len = read_tsv(ap_ + ".length_dist.tsv")
        self.asm_functions = read_tsv(ap_ + ".functions.tsv")
        cp_ = os.path.join(an, "crispr", sample)
        self.cr_arrays = read_tsv(cp_ + ".crispr_arrays.tsv")
        self.cr_units = read_tsv(cp_ + ".crispr_units.tsv")
        self.cr_ends = read_tsv(cp_ + ".crispr_ends.tsv")
        self.cr_lengths = read_tsv(cp_ + ".crispr_lengths.tsv")
        self.cr_genes = read_tsv(cp_ + ".crispr_genes.tsv")
        self.cr_cov = read_tsv(cp_ + ".crispr_coverage.tsv")
        self.cr_tracr = read_tsv(cp_ + ".crispr_tracr.tsv")
        self.feat_counts = ap_ + ".feature_counts.tsv"
        if not os.path.exists(self.feat_counts):
            self.feat_counts = os.path.join(an, "bins", "%s.feature_counts.tsv" % sample)
        self.novel_sum = kv(os.path.join(an, "novel", "%s.assembly.summary.tsv" % sample), "category", "reads")
        self.cand = {b: self._bed(os.path.join(an, "novel", "%s.candidates.%s.bed" % (sample, b)))
                     for b in ("intergenic", "antisense_to_CDS")}
        self.rfam_tbl = self._tbl(os.path.join(an, "novel", "%s.candidates.rfam.tbl" % sample))
        self.mfe = self._rnafold(os.path.join(an, "novel", "%s.candidates.no_rfam.RNAfold.txt" % sample))

    @staticmethod
    def _rows(path):
        if not os.path.exists(path):
            return []
        with open(path) as fh:
            return [ln.rstrip("\n").split("\t", 1) for ln in fh if "\t" in ln]

    @staticmethod
    def _bed(path):
        out = []
        if os.path.exists(path):
            with open(path) as fh:
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) >= 7:
                        out.append((f[0], int(f[1]), int(f[2]), f[5], num(f[6]), num(f[7]) if len(f) > 7 else 0))
        return out

    @staticmethod
    def _two(path, header=False):
        """two-column table (key, number) -> {key: float}"""
        out = {}
        if os.path.exists(path):
            with open(path) as fh:
                for i, ln in enumerate(fh):
                    f = ln.split()
                    if len(f) < 2 or (header and i == 0):
                        continue
                    try:
                        out[" ".join(f[:-1])] = out.get(" ".join(f[:-1]), 0) + float(f[-1])
                    except ValueError:
                        pass
        return out

    @staticmethod
    def _kraken(path):
        """Kraken2 report -> list of dicts (clade, taxon, distinct minimizers or None, rank, taxid, name)"""
        out = []
        if os.path.exists(path):
            with open(path) as fh:
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) >= 8:        # --report-minimizer-data
                        out.append({"clade": num(f[1]), "taxon": num(f[2]), "distinct": num(f[4]),
                                    "rank": f[5].strip(), "taxid": f[6].strip(), "name": f[7].strip()})
                    elif len(f) >= 6:
                        out.append({"clade": num(f[1]), "taxon": num(f[2]), "distinct": None,
                                    "rank": f[3].strip(), "taxid": f[4].strip(), "name": f[5].strip()})
        return out

    @staticmethod
    def _rnafold(path):
        """RNAfold output: >name / sequence / structure ( -12.30) -> {name: MFE}"""
        out, name = {}, None
        if os.path.exists(path):
            with open(path) as fh:
                for ln in fh:
                    ln = ln.strip()
                    if ln.startswith(">"):
                        name = ln[1:].split()[0]
                    elif name and ln.endswith(")") and "(" in ln:
                        m = re.search(r"\(\s*(-?[0-9.]+)\)\s*$", ln)
                        if m:
                            out[name] = float(m.group(1))
                            name = None
        return out

    @staticmethod
    def _tbl(path):
        """cmscan --tblout: best hit per query -> (family, accession, description, evalue)"""
        best = {}
        if os.path.exists(path):
            with open(path) as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    f = ln.split(None, 17)
                    if len(f) < 17:
                        continue
                    q, ev = f[2], num(f[15], 1.0)
                    if q not in best or ev < best[q][3]:
                        best[q] = (f[0], f[1], f[17].strip() if len(f) > 17 else "", ev)
        return best

    # convenience
    def rs(self, k):
        return num(self.read_stats.get(k))


# ----------------------------------------------------------------------------- pages
def page_overview(pdf, sample, td, outdir):
    fig = page(pdf, sample, "overview", "sRNA types: %s · output: %s" % (", ".join(d.t for d in td), outdir),
               explain="overview")
    # key numbers table
    rows = []
    def line(label, f, nd=0, rel=None):
        r = [label]
        for d in td:
            v = f(d)
            r.append(fmt(v, nd) + ("" if rel is None or v is None else "  (%s)" % pct(v, rel(d))))
        rows.append(r)
    raw = lambda d: d.rs("raw") or d.qc.get("raw_reads")
    line("raw reads", raw)
    line("after trimming", lambda d: d.rs("trimmed"), rel=raw)
    line("adapter dimers (raw)", lambda d: d.qc.get("adapter_dimers_raw"), rel=raw)
    line("< 20 nt (never human)", lambda d: d.rs("below_host_min_len_20nt"), rel=lambda d: d.rs("trimmed"))
    line("human", lambda d: d.rs("human_total"), rel=lambda d: d.rs("trimmed"))
    line("  rescued by assembly", lambda d: d.rs("host_rescued_by_assembly"))
    line("non-human", lambda d: d.rs("nonhuman"), rel=lambda d: d.rs("trimmed"))
    nh = lambda k: (lambda d: d.nonhuman.get(k))
    line("  rRNA (RiboDetector)", nh("rRNA_reads"), rel=lambda d: d.rs("nonhuman"))
    line("  tRNA", nh("tRNA"), rel=lambda d: d.rs("nonhuman"))
    line("  Rfam", nh("Rfam_classified"), rel=lambda d: d.rs("nonhuman"))
    line("  miRNA-like", nh("miRNA_like"))
    line("  unknown", nh("unknown"), rel=lambda d: d.rs("nonhuman"))
    line("mapped to own assembly", lambda d: d.novel_sum.get("mapped_reads"), rel=lambda d: d.novel_sum.get("input_reads"))
    line("novel intergenic / antisense loci",
         lambda d: "%s / %s" % (fmt(d.novel_sum.get("novel_intergenic_candidates")),
                                fmt(d.novel_sum.get("novel_antisense_candidates"))))
    line("human miRNA reads", lambda d: d.human.get("human_miRNA"))
    ax = fig.add_axes([0.03, 0.42, 0.55, 0.44])
    table(ax, rows, ["metric"] + [d.t for d in td], col_widths=[0.42] + [0.29] * len(td), fontsize=7.5,
          title="key numbers", rowscale=1.15)

    # read fate stacked bars (fraction of raw reads)
    ax = fig.add_axes([0.66, 0.52, 0.31, 0.34])
    # read fate = exclusive partition of the raw reads (sums to 1); reads < HOST_MIN_LEN are part of the
    # non-human classes (they are not tested against human), so they are not a segment of their own
    cats = ["too short / dimer", "human", "rRNA", "tRNA", "Rfam", "miRNA-like", "unknown"]
    for i, d in enumerate(td):
        R = raw(d) or 1
        vals = {
            "too short / dimer": max(R - d.rs("trimmed"), 0),
            "human": d.rs("human_total"),
            "rRNA": num(d.nonhuman.get("rRNA_reads")),
            "tRNA": num(d.nonhuman.get("tRNA")),
            "Rfam": num(d.nonhuman.get("Rfam_classified")),
            "miRNA-like": num(d.nonhuman.get("miRNA_like")),
            "unknown": num(d.nonhuman.get("unknown")),
        }
        left = 0
        for c in cats:
            w = vals[c] / R
            ax.barh(i, w, left=left, color=CLASS_COLORS.get(c, "#888"), label=c if i == 0 else None)
            left += w
    ax.set_yticks(range(len(td)))
    ax.set_yticklabels([d.t for d in td])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of raw reads")
    ax.set_title("read fate (fate of every raw read)", loc="left")
    ax.legend(frameon=False, ncol=4, bbox_to_anchor=(0, -0.18), loc="upper left")

    # inputs
    ax = fig.add_axes([0.03, 0.04, 0.94, 0.30])
    rows = []
    for d in td:
        if d.inputs.get("raw_fastq"):
            rows.append(["raw_fastq (%s)" % d.t, d.inputs["raw_fastq"], d.inputs.get("date", "")])
    for k in ("assembly", "gff", "contig2bin", "gtdbtk", "metaT_features", "metaT_contigs",
              "metaG_features", "metaG_contigs", "config"):
        vals = sorted({d.inputs.get(k, "") for d in td} - {""})
        if vals:
            rows.append([k, " | ".join(vals), ""])
    table(ax, rows, ["input", "path", "run"], col_widths=[0.14, 0.74, 0.12], fontsize=6,
          title="inputs (log_cluster_jobs/<type>/inputs.tsv)", left=(0, 1), rowscale=1.0)
    pdf.savefig(fig)
    plt.close(fig)


QC_STATE = {"rules": [], "exclude": {}, "outdir": "", "result": {}}
FASTQC_MODULES = ["basic_statistics", "per_base_sequence_quality", "per_tile_sequence_quality",
                  "per_sequence_quality_scores", "per_base_sequence_content", "per_sequence_gc_content",
                  "per_base_n_content", "sequence_length_distribution", "sequence_duplication_levels",
                  "overrepresented_sequences", "adapter_content"]


def qc_verdicts(sample, td):
    """evaluate the QC rules for all types of the sample (cached)"""
    import srna_qc as Q
    if sample not in QC_STATE["result"]:
        res = {}
        for d in td:
            m = Q.library_metrics(d, QC_STATE["outdir"], sample, d.t)
            res[d.t] = (m,) + Q.evaluate(m, sample, d.t, QC_STATE["rules"], QC_STATE["exclude"])
        QC_STATE["result"][sample] = res
    return QC_STATE["result"][sample]


def page_qc_verdict(pdf, sample, td):
    import srna_qc as Q
    res = qc_verdicts(sample, td)
    fig = page(pdf, sample, "QC verdict", "rules: %s%s" % (QC_STATE.get("rules_path") or "scripts/qc_rules.tsv",
               " · manual decisions: %s" % QC_STATE["exclude_path"] if QC_STATE.get("exclude_path") else ""),
               explain="qcverdict")
    n = len(td)
    for i, d in enumerate(td):
        m, lv, hits = res[d.t]
        x0 = 0.03 + i * (0.94 / n)
        w = 0.94 / n - 0.02
        fig.text(x0, 0.795, "%s:" % d.t, fontsize=11, fontweight="bold", va="top")
        fig.text(x0 + 0.09, 0.795, " %s " % lv, fontsize=11, fontweight="bold", va="top", color="white",
                 bbox=dict(boxstyle="round,pad=0.3", fc=Q.LEVEL_COLOR[lv], ec="none"))
        ax = fig.add_axes([x0, 0.52, w, 0.24])
        rows = [[h[1], h[0], R_short(h[5], 60), "" if not h[3] else "%s %s %s" % (_fmtv(h[2]), h[3], h[4])]
                for h in hits]
        table(ax, rows or [["PASS", "-", "no rule triggered", ""]], ["level", "rule", "reason", "value"],
              fontsize=6.5, left=(0, 1, 2), col_widths=[0.11, 0.17, 0.50, 0.22], title="triggered rules")
        for (r, c), cell in ax.tables[0].get_celld().items() if ax.tables else []:
            if r > 0 and c == 0:
                cell.set_facecolor(Q.LEVEL_COLOR.get(cell.get_text().get_text(), "#ffffff") + "55")
    # FastQC module grid
    ax = fig.add_axes([0.22, 0.06, 0.55, 0.40])
    cols = [(d.t, st) for d in td for st in ("raw", "trimmed")]
    any_ = False
    for j, (t, st) in enumerate(cols):
        m = res[t][0]
        for i_, mod in enumerate(FASTQC_MODULES):
            v = m.get("fastqc_%s_%s" % (st, mod), "")
            any_ = any_ or bool(v)
            ax.add_patch(plt.Rectangle((j, i_), 0.95, 0.9, color=Q.LEVEL_COLOR.get(v, "#eeeeee")))
            ax.text(j + 0.47, i_ + 0.45, v or "–", ha="center", va="center", fontsize=6.5,
                    color="white" if v else "#999")
    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(FASTQC_MODULES), 0)
    ax.set_xticks([j + 0.47 for j in range(len(cols))])
    ax.set_xticklabels(["%s\n%s" % (t, st) for t, st in cols], fontsize=7)
    ax.set_yticks([i + 0.45 for i in range(len(FASTQC_MODULES))])
    ax.set_yticklabels([x.replace("_", " ") for x in FASTQC_MODULES], fontsize=7)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("FastQC modules (raw / trimmed reads)%s" % ("" if any_ else " — no FastQC zip files found"),
                 loc="left", fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


def _fmtv(v):
    try:
        v = float(v)
        return "%.3g" % v if abs(v) < 1e4 else "{:,.0f}".format(v)
    except (TypeError, ValueError):
        return str(v)


def R_short(s, n):
    return short(s, n)


def write_qc_flags(path, sample, td):
    res = qc_verdicts(sample, td)
    with open(path, "w") as out:
        out.write("sample\ttype\tqc\trule\tlevel\tvalue\tthreshold\treason\n")
        for d in td:
            m, lv, hits = res[d.t]
            if not hits:
                out.write("%s\t%s\t%s\t-\t-\t-\t-\t-\n" % (sample, d.t, lv))
            for h in hits:
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s %s\t%s\n" % (sample, d.t, lv, h[0], h[1], _fmtv(h[2]), h[3],
                                                                    h[4], h[5]))


def page_qc(pdf, sample, td):
    fig = page(pdf, sample, "QC and trimming", explain="qc")
    rows = []
    keys = [("raw_reads", "raw reads"), ("cutadapt_reads_with_adapter", "reads with adapter"),
            ("adapter_dimers_raw", "adapter dimers"), ("cutadapt_too_short", "too short (< MINLEN)"),
            ("cutadapt_too_long", "too long (> MAXLEN)"), ("trimmed_reads", "trimmed reads"),
            ("distinct_sequences", "distinct sequences"), ("singleton_sequences", "singletons (of distinct)"),
            ("top10_sequences_reads", "top-10 sequences")]
    keys += sorted({(k, k.replace("contaminant_", "contaminant ")) for d in td for k in d.qc if k.startswith("contaminant_")})
    for k, lab in keys:
        r = [lab]
        for d in td:
            v = d.qc.get(k)
            base = d.qc.get("distinct_sequences") if k == "singleton_sequences" else \
                (d.qc.get("trimmed_reads") if k.startswith("contaminant") or k.startswith("top10") or k == "distinct_sequences"
                 else d.qc.get("raw_reads"))
            r.append("–" if v is None else "%s  (%s)" % (fmt(v), pct(v, base, 2)))
        rows.append(r)
    ax = fig.add_axes([0.03, 0.45, 0.45, 0.42])
    table(ax, rows, ["metric"] + [d.t for d in td], col_widths=[0.36] + [0.32] * len(td), title="library QC")

    # trimming fate bars
    ax = fig.add_axes([0.60, 0.60, 0.36, 0.25])
    labs = ["trimmed", "too short", "adapter dimers"]
    series = []
    for d in td:
        R = d.qc.get("raw_reads") or 1
        series.append((d.t, [num(d.qc.get("trimmed_reads")) / R, num(d.qc.get("cutadapt_too_short")) / R,
                             num(d.qc.get("adapter_dimers_raw")) / R]))
    grouped_barh(ax, labs, series, TYPE_COLORS, "fraction of raw reads")
    ax.set_xlim(0, 1)

    # overrepresented sequences
    n = len(td)
    for i, d in enumerate(td):
        ax = fig.add_axes([0.03 + i * (0.94 / n), 0.04, 0.94 / n - 0.02, 0.36])
        rows = [[short(r["sequence"], 52), r["length"], fmt(num(r["count"])), "%.2f%%" % (100 * num(r["fraction"]))]
                for r in d.over[:12]]
        table(ax, rows, ["sequence", "nt", "reads", "%"], col_widths=[0.66, 0.08, 0.15, 0.11], fontsize=6,
              title="%s — most frequent trimmed sequences" % d.t)
    pdf.savefig(fig)
    plt.close(fig)


# built-in default; normally read from mirna_origin_markers.tsv (next to this script or --mirna-markers)
MIRNA_ORIGIN_DEFAULT = [
    ("intestinal / epithelial", "#009E73",
     ["miR-192-5p", "miR-194-5p", "miR-215-5p", "miR-200a-3p", "miR-200b-3p", "miR-200c-3p", "miR-141-3p",
      "miR-429", "miR-375-3p", "miR-148a-3p", "miR-196b-5p"], []),
    ("blood", "#D55E00", ["miR-451a", "miR-486-5p", "miR-144-3p", "miR-144-5p", "miR-16-5p", "miR-92a-3p"], []),
    ("immune", "#0072B2", ["miR-142-3p", "miR-142-5p", "miR-150-5p", "miR-223-3p", "miR-155-5p", "miR-146a-5p"], []),
    ("muscle (diet?)", "#CC79A7", ["miR-1-3p", "miR-206", "miR-133a-3p", "miR-133b"], []),
]
MIRNA_ORIGIN = []           # (group, color, [miRNAs], [references]); filled by load_markers()
MIRNA_GENERAL_REFS = []     # references of "_..." rows
MIRNA_MARKER_SOURCE = "built-in default"


def load_markers(path):
    """read the marker TSV (group, color, mirnas, references); fall back to the built-in list"""
    global MIRNA_ORIGIN, MIRNA_GENERAL_REFS, MIRNA_MARKER_SOURCE
    MIRNA_ORIGIN, MIRNA_GENERAL_REFS = [], []
    if path and os.path.exists(path):
        palette = ["#009E73", "#D55E00", "#0072B2", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"]
        for r in read_tsv(path):
            g = (r.get("group") or "").strip()
            if not g:
                continue
            refs = [x.strip() for x in (r.get("references") or "").split("|") if x.strip()]
            if g.startswith("_"):
                MIRNA_GENERAL_REFS.append((g.lstrip("_").strip(), refs))
                continue
            mirnas = [re.sub(r"^hsa-", "", x.strip()) for x in (r.get("mirnas") or "").split(",") if x.strip()]
            color = (r.get("color") or "").strip() or palette[len(MIRNA_ORIGIN) % len(palette)]
            MIRNA_ORIGIN.append((g, color, mirnas, refs))
        MIRNA_MARKER_SOURCE = path
    if not MIRNA_ORIGIN:
        MIRNA_ORIGIN = [tuple(x) for x in MIRNA_ORIGIN_DEFAULT]
        MIRNA_MARKER_SOURCE = "built-in default (no marker file)"


def all_references():
    """numbered, de-duplicated list of all marker references -> (list, {ref: number})"""
    refs = []
    for _g, _c, _m, rr in MIRNA_ORIGIN:
        refs += [r for r in rr if r not in refs]
    for _g, rr in MIRNA_GENERAL_REFS:
        refs += [r for r in rr if r not in refs]
    return refs, {r: i + 1 for i, r in enumerate(refs)}


def mirna_origin(rows):
    """human miRNA counts -> {group: reads}, plus 'other'"""
    cnt = {}
    for r in rows:
        name = re.sub(r"^hsa-", "", r.get("feature", ""))
        cnt[name] = cnt.get(name, 0) + num(r.get("count"))
    out, used = {}, set()
    for g, _c, markers, _r in MIRNA_ORIGIN:
        out[g] = sum(cnt.get(m, 0) for m in markers)
        used |= set(markers)
    out["other"] = sum(v for k, v in cnt.items() if k not in used)
    return out


def origin_of(feature):
    name = re.sub(r"^hsa-", "", feature)
    for g, _c, markers, _r in MIRNA_ORIGIN:
        if name in markers:
            return g.split(" ")[0]
    return ""


def page_human(pdf, sample, td):
    fig = page(pdf, sample, "human reads",
               "human = exact host-genome/transcript match (≥ 20 nt) not explained better by the sample's own assembly",
               explain="human")
    rows = []
    for lab, f in [("trimmed reads", lambda d: d.rs("trimmed")),
                   ("host genome hits", lambda d: d.rs("host_genome_mapped")),
                   ("host transcript hits", lambda d: d.rs("host_tx_mapped")),
                   ("rescued by assembly", lambda d: d.rs("host_rescued_by_assembly")),
                   ("human", lambda d: d.rs("human_total")),
                   ("  miRNA", lambda d: d.human.get("human_miRNA")),
                   ("  tRNA", lambda d: d.human.get("human_tRNA")),
                   ("  rRNA", lambda d: d.human.get("human_rRNA")),
                   ("  other ncRNA", lambda d: d.human.get("human_otherRNA")),
                   ("  unassigned", lambda d: d.human.get("unknown_after_all"))]:
        rows.append([lab] + ["%s  (%s)" % (fmt(f(d)), pct(num(f(d)), d.rs("human_total"))) if lab.startswith("  ")
                             else fmt(f(d)) for d in td])
    rows.append(["human per million trimmed"] + [fmt(1e6 * d.rs("human_total") / (d.rs("trimmed") or 1)) for d in td])
    ax = fig.add_axes([0.03, 0.50, 0.38, 0.38])
    table(ax, rows, ["human reads"] + [d.t for d in td], col_widths=[0.40] + [0.30] * len(td))

    ax = fig.add_axes([0.50, 0.55, 0.46, 0.31])
    for i, d in enumerate(td):
        xs = [int(r["length"]) for r in d.hlen]
        tot = sum(num(r["human_reads"]) for r in d.hlen) or 1
        ys = [100 * num(r["human_reads"]) / tot for r in d.hlen]
        if xs:
            ax.plot(xs, ys, marker="o", ms=2, lw=1.2, color=TYPE_COLORS[i], label=d.t)
    ax.set_xlabel("read length (nt)")
    ax.set_ylabel("% of human reads")
    ax.set_title("human read length (grey: miRNA size 21–23 nt)", loc="left")
    ax.axvspan(20.5, 23.5, color="#eee", zorder=0)
    ax.legend(frameon=False)

    # miRNA origin: stacked bars (share of human miRNA reads) + reads per group
    ax = fig.add_axes([0.06, 0.14, 0.24, 0.24])
    groups = [g for g, _c, _m, _r in MIRNA_ORIGIN] + ["other"]
    colors = {g: c for g, c, _m, _r in MIRNA_ORIGIN}
    colors["other"] = "#bbbbbb"
    orig = [mirna_origin(d.mirna) for d in td]
    for i, (d, o) in enumerate(zip(td, orig)):
        tot = sum(o.values())
        left = 0
        for g in groups:
            w = o[g] / tot if tot else 0
            ax.barh(i, w, left=left, color=colors[g], label=g if i == 0 else None)
            left += w
        ax.text(1.01, i, "n = %s" % fmt(tot), va="center", fontsize=7)
    ax.set_yticks(range(len(td)))
    ax.set_yticklabels([d.t for d in td])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of human miRNA reads")
    ax.set_title("human miRNA origin (marker miRNAs)", loc="left")
    _refs, refno = all_references()
    handles, labels = ax.get_legend_handles_labels()
    lab2 = []
    for l in labels:
        rr = [refno[r] for g, _c, _m, refs_ in MIRNA_ORIGIN if g == l for r in refs_]
        lab2.append("%s [%s]" % (l, ",".join(map(str, rr))) if rr else l)
    ax.legend(handles, lab2, frameon=False, fontsize=6.3, ncol=2, bbox_to_anchor=(0, -0.20), loc="upper left")
    fig.text(0.06, 0.012, "marker sets: %s · [n] = reference on the last page" % short(MIRNA_MARKER_SOURCE, 90),
             fontsize=6, color="#777")

    n = len(td)
    w = 0.60 / n
    for i, d in enumerate(td):
        ax = fig.add_axes([0.37 + i * w, 0.03, w - 0.02, 0.38])
        top = sorted(d.mirna, key=lambda r: -num(r.get("count")))[:13]
        tot = sum(num(r.get("count")) for r in d.mirna) or 1
        rows = [[r["feature"], fmt(num(r["count"]), 1), "%.1f%%" % (100 * num(r["count"]) / tot), origin_of(r["feature"])]
                for r in top]
        table(ax, rows, ["miRNA", "reads", "%", "origin"], col_widths=[0.42, 0.18, 0.16, 0.24], fontsize=6.3,
              title="%s — top human miRNAs" % d.t, left=(0, 3), rowscale=1.15)
    pdf.savefig(fig)
    plt.close(fig)


def page_nonhuman(pdf, sample, td):
    fig = page(pdf, sample, "non-human reads",
               "rRNA: RiboDetector → Kraken2 (GTDB) + Bracken; non-rRNA: tRNA → miRNA-like → Rfam; "
               "Bracken fractions are of rRNA reads assigned at that rank", explain="nonhuman")
    # Bracken genus / family
    for j, lv in enumerate(("G", "F")):
        dicts = [{r["name"]: num(r["fraction_total_reads"]) for r in d.bracken[lv]} for d in td]
        if not any(dicts):
            note(fig.add_axes([0.06 + j * 0.49, 0.52, 0.40, 0.33]), "no Bracken %s output" % lv)
            continue
        labels = top_union(dicts, 15)
        ax = fig.add_axes([0.17 + j * 0.49, 0.52, 0.30, 0.33])
        grouped_barh(ax, labels, [(d.t, [dd.get(l, 0) for l in labels]) for d, dd in zip(td, dicts)],
                     TYPE_COLORS, "Bracken %s (fraction)" % {"G": "genus", "F": "family"}[lv])
    # Rfam types
    dicts = [frac(d.rfam_type) for d in td]
    if any(dicts):
        labels = top_union(dicts, 10)
        ax = fig.add_axes([0.17, 0.06, 0.30, 0.36])
        grouped_barh(ax, labels, [(d.t, [dd.get(l, 0) for l in labels]) for d, dd in zip(td, dicts)],
                     TYPE_COLORS, "Rfam types (fraction of Rfam reads)", log=True)
    else:
        note(fig.add_axes([0.06, 0.06, 0.40, 0.36]), "no Rfam type table (RFAM_FAMILY_TXT not set?)")
    # tRNA amino acids
    dicts = [frac(d.trna_aa) for d in td]
    if any(dicts):
        labels = sorted({k for dd in dicts for k in dd if k not in ("unknown", "")},
                        key=lambda k: -max(dd.get(k, 0) for dd in dicts))[:21]
        ax = fig.add_axes([0.58, 0.06, 0.39, 0.36])
        n = len(td)
        w = 0.8 / n
        for i, (d, dd) in enumerate(zip(td, dicts)):
            ax.bar([x + (i - (n - 1) / 2) * w for x in range(len(labels))], [dd.get(l, 0) for l in labels],
                   width=w, color=TYPE_COLORS[i], label=d.t)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90)
        ax.set_title("tRNA reads by amino acid (fraction)", loc="left")
        ax.legend(frameon=False)
    pdf.savefig(fig)
    plt.close(fig)


RFAM_FAMILY = {}       # RFxxxxx -> (id, description, type); from --rfam-family


def load_rfam_family(path):
    if not path or not os.path.exists(path):
        return
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 19 and c[0].startswith("RF"):
                parts = [x.strip() for x in c[18].split(";") if x.strip()]
                if len(parts) > 1 and parts[0] == "Gene":
                    parts = parts[1:]
                RFAM_FAMILY[c[0]] = (c[1], c[3], ";".join(parts) or "unknown")


def element_type(gff_type, product, rfam):
    if rfam in RFAM_FAMILY:
        t = RFAM_FAMILY[rfam][2]
        return "Gene (other ncRNA)" if t == "Gene" else t
    p = product.lower()
    for key, t in (("riboswitch", "Cis-reg;riboswitch"), ("leader", "Cis-reg;leader"),
                   ("thermometer", "Cis-reg;thermoregulator"), ("rnase p", "ribozyme"), ("srp", "SRP"),
                   ("tmrna", "tmRNA")):
        if key in p:
            return t
    return "Cis-reg" if gff_type == "regulatory_region" else gff_type


def read_elements(bed, types=("ncRNA", "regulatory_region", "tmRNA")):
    """annotation BED (gff2bed) or GFF -> {feature_id: (gff_type, product, rfam_acc)} for ncRNA-like features"""
    out = {}
    if not bed or not os.path.exists(bed):
        return out
    op = gzip.open if bed.endswith(".gz") else open
    with op(bed, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 10 and f[7] in types:            # gff2bed
                ftype, attr = f[7], f[9]
            elif len(f) >= 9 and f[2] in types:           # GFF
                ftype, attr = f[2], f[8]
            else:
                continue
            a = dict(x.split("=", 1) for x in attr.split(";") if "=" in x)
            fid = a.get("ID", "")
            m = re.search(r"RFAM:(RF\d+)", a.get("Dbxref", ""))
            prod = urllib.parse.unquote(a.get("product") or a.get("Name") or "unknown")
            prod = re.sub(r"^\((?:5'|3'|5' and 3'|5'/3')\s*truncated\)\s*", "", prod)   # merge partial copies
            out[fid] = (ftype, prod, m.group(1) if m else "")
    return out


def page_ncrna(pdf, sample, td, base):
    bed = os.path.join(base, "Analysis", "annotation", "%s.bed" % sample)
    elems = read_elements(bed)
    if not elems or not any(os.path.exists(d.feat_counts) for d in td):
        return False
    fig = page(pdf, sample, "annotated ncRNAs and regulatory elements",
               "Bakta ncRNA / regulatory_region / tmRNA features (Rfam); reads from the assembly mapping "
               "(weighted 1/hits, sense / antisense); type from Rfam family.txt", explain="ncrna")
    # per type and per product
    n_annot = collections.Counter(element_type(t, p, r) for t, p, r in elems.values())
    bytype, byprod = [], []
    for d in td:
        bt = collections.defaultdict(lambda: [0.0, 0.0, 0])
        bp = collections.defaultdict(lambda: [0.0, 0.0, 0, ""])
        for r in read_tsv(d.feat_counts):
            e = elems.get(r["feature_id"])
            if not e:
                continue
            t = element_type(*e)
            se, an_ = num(r["sense_reads"]), num(r["antisense_reads"])
            bt[t][0] += se
            bt[t][1] += an_
            bt[t][2] += 1
            bp[e[1]][0] += se
            bp[e[1]][1] += an_
            bp[e[1]][2] += 1
            bp[e[1]][3] = t
        bytype.append(bt)
        byprod.append(bp)
    types = sorted(set(n_annot) | {t for bt in bytype for t in bt},
                   key=lambda t: -max([bt.get(t, [0])[0] for bt in bytype] + [0]))
    rows = []
    for t in types[:14]:
        rows.append([t, fmt(n_annot.get(t, 0))] +
                    ["%s / %s" % (fmt(bt[t][0]), fmt(bt[t][1])) if t in bt else "0" for bt in bytype])
    rows.append(["all", fmt(sum(n_annot.values()))] +
                ["%s / %s" % (fmt(sum(v[0] for v in bt.values())), fmt(sum(v[1] for v in bt.values()))) for bt in bytype])
    ax = fig.add_axes([0.03, 0.47, 0.50, 0.40])
    table(ax, rows, ["Rfam type", "annotated"] + ["%s sense / anti" % d.t for d in td],
          col_widths=[0.34, 0.14] + [0.52 / len(td)] * len(td), fontsize=6.6, title="reads per element type")
    # bars: sense reads per type (log)
    ax = fig.add_axes([0.70, 0.50, 0.27, 0.34])
    lab = types[:10]
    grouped_barh(ax, lab, [(d.t, [bt.get(t, [0])[0] + 0.1 for t in lab]) for d, bt in zip(td, bytype)],
                 TYPE_COLORS, "sense reads per type (log)", log=True)
    # top products
    allp = collections.defaultdict(float)
    for bp in byprod:
        tot = sum(v[0] + v[1] for v in bp.values()) or 1
        for k, v in bp.items():
            allp[k] = max(allp[k], (v[0] + v[1]) / tot)
    top = sorted(allp, key=lambda k: -allp[k])[:13]
    rows = []
    for k in top:
        t = next((bp[k][3] for bp in byprod if k in bp), "")
        r = [short(k, 44), short(t, 22)]
        for bp in byprod:
            v = bp.get(k, [0, 0, 0, ""])
            r += [fmt(v[0]), fmt(v[1])]
        rows.append(r)
    ax = fig.add_axes([0.03, 0.03, 0.94, 0.38])
    cols = ["element (product)", "type"]
    for d in td:
        cols += ["%s sense" % d.t, "%s anti" % d.t]
    table(ax, rows, cols, col_widths=[0.32, 0.18] + [0.50 / (2 * len(td))] * (2 * len(td)), fontsize=6.4,
          title="top elements (features summed by product)", left=(0, 1))
    pdf.savefig(fig)
    plt.close(fig)


def page_asm_cmc(pdf, sample, td):
    if not any(d.asm_summary or d.asm_category for d in td):
        return False
    fig = page(pdf, sample, "assembly mapping",
               "module 'assembly': uncollapsed non-human reads → bowtie on the contigs (-v ASSEMBLY_MISMATCHES, "
               "≤ ASSEMBLY_MAX_HITS best hits, NH tag) → features of the annotation (sense before antisense)",
               explain="asm_cmc")
    rows = []
    for lab, k in (("input reads", "input_reads"), ("mapped reads", "mapped_reads"), ("not mapped", "unmapped_reads"),
                   ("alignments", "alignments")):
        rows.append([lab] + ["%s%s" % (fmt(num(d.asm_summary.get(k))),
                                       "  (%s)" % pct(num(d.asm_summary.get(k)), num(d.asm_summary.get("input_reads")))
                                       if k in ("mapped_reads", "unmapped_reads") else "") for d in td])
    cats = top_union([frac(d.asm_category) for d in td], 12)
    for c in cats:
        rows.append(["  " + c] + ["%s  (%s)" % (fmt(d.asm_category.get(c, 0)),
                                               pct(d.asm_category.get(c, 0), sum(d.asm_category.values())))
                                  for d in td])
    def as_share(cat):
        a, b = cat.get("antisense_CDS", 0), cat.get("CDS", 0)
        return pct(a, a + b) if a + b else "–"
    if any(d.asm_category_nr for d in td):
        rows.append(["antisense on CDS: all / no rRNA"] + ["%s / %s" % (as_share(d.asm_category), as_share(
            d.asm_category_nr)) if d.asm_category_nr else as_share(d.asm_category) for d in td])
    ax = fig.add_axes([0.03, 0.42, 0.42, 0.45])
    table(ax, rows, ["reads"] + [d.t for d in td], col_widths=[0.40] + [0.30] * len(td), fontsize=6.6,
          title="mapping and reads per feature category (weighted)", rowscale=1.12)
    # length distributions
    ax = fig.add_axes([0.53, 0.55, 0.43, 0.29])
    for i, d in enumerate(td):
        for col, ls, lab in (("mapped_weighted", "-", "mapped"), ("unmapped", "--", "not mapped")):
            tot = sum(num(r[col]) for r in d.asm_len)
            if tot:
                xs = [int(r["length"]) for r in d.asm_len]
                ax.plot(xs, [100 * num(r[col]) / tot for r in d.asm_len], ls, color=TYPE_COLORS[i], lw=1.2,
                        label="%s %s" % (d.t, lab))
    ax.set_xlabel("read length (nt)")
    ax.set_ylabel("% of reads")
    ax.set_title("read length: mapped (solid) vs not mapped (dashed)", loc="left")
    ax.legend(frameon=False, ncol=2, fontsize=6.5)
    # top protein-coding genes (by product / gene)
    per = []
    for d in td:
        agg = collections.defaultdict(lambda: [0.0, 0.0])
        for r in read_tsv(d.feat_counts):
            if r.get("type") != "CDS" or r.get("flag"):     # flag: spurious ORF opposite a rescued rRNA gene
                continue
            k = r.get("gene") or r.get("product") or r["feature_id"]
            if k.lower() == "hypothetical protein":
                continue
            agg[k][0] += num(r["sense_reads"])
            agg[k][1] += num(r["antisense_reads"])
        per.append(agg)
    keys = top_union([{k: v[0] + v[1] for k, v in pf.items()} for pf in per], 12)
    rows = []
    for k in keys:
        r = [short(k, 30)]
        for pf in per:
            a_, b_ = pf.get(k, (0, 0))
            r += [fmt(a_), fmt(b_)]
        rows.append(r)
    cols = ["protein-coding gene"]
    for d in td:
        cols += ["%s sense" % d.t.replace("sRNA-", ""), "%s anti" % d.t.replace("sRNA-", "")]
    ax = fig.add_axes([0.03, 0.03, 0.46, 0.32])
    table(ax, rows, cols, col_widths=[0.36] + [0.64 / (2 * len(td))] * (2 * len(td)), fontsize=6.2,
          title="top protein-coding genes (hypothetical proteins left out)", rowscale=1.1)
    # top KEGG pathways
    dicts = []
    for d in td:
        dicts.append({r["term"]: num(r["sense_reads"]) for r in d.asm_functions if r.get("db") == "kegg_pathway"})
    ax = fig.add_axes([0.53, 0.03, 0.44, 0.32])
    if any(dicts):
        keys = top_union(dicts, 12)
        rows = [[k] + ["%s (%s)" % (fmt(dd.get(k, 0)), pct(dd.get(k, 0), sum(dd.values()))) for dd in dicts]
                for k in keys]
        table(ax, rows, ["KEGG pathway"] + ["%s sense" % d.t for d in td],
              col_widths=[0.30] + [0.70 / len(td)] * len(td), fontsize=6.2,
              title="top KEGG pathways of reads on CDS (a read counts for every pathway of its gene)", rowscale=1.1)
    else:
        note(ax, "no functional table (functions.tsv)")
    pdf.savefig(fig)
    plt.close(fig)


def page_crispr(pdf, sample, td):
    if not any(d.cr_arrays for d in td):
        return False
    fig = page(pdf, sample, "CRISPR arrays and crRNAs",
               "module 'crispr': reads on Bakta CRISPR arrays (assembly BAM, 1/NH); transcribed strand from the "
               "reads; Cas type from cas genes within 10 kb", explain="crispr")
    # union of top arrays
    score = collections.defaultdict(float)
    info = {}
    for d in td:
        for r in d.cr_arrays:
            score[r["array_id"]] = max(score[r["array_id"]], num(r["reads"]))
            info.setdefault(r["array_id"], r)
    top = [k for k in sorted(score, key=lambda k: -score[k]) if score[k] > 0][:10]
    # arrays with a tracrRNA candidate are always listed
    top += [k for k in sorted({x["array_id"] for d in td for x in getattr(d, "cr_tracr", [])
                               if x.get("tracr_candidate") == "yes"}, key=lambda k: -score.get(k, 0))
            if k not in top and k in info][:4]
    rows = []
    for k in top:
        r = info[k]
        row = [k.split("_")[-1], short(tax_tail(r.get("classification", "")) if r.get("classification") else r["bin"], 24),
               r["cas_type"], r["repeats"]]
        for d in td:
            dr = next((x for x in d.cr_arrays if x["array_id"] == k), None)
            row.append("%s (%s%s)" % (fmt(num(dr["reads"])), dr["transcribed_strand"],
                                       "%.0f%%" % (100 * num(dr["strand_fraction"]))) if dr else "–")
        # tracrRNA candidate (anti-repeat) in any fraction: anti-repeat length and reads per fraction
        tc = [(d.t, x) for d in td for x in getattr(d, "cr_tracr", [])
              if x["array_id"] == k and x.get("tracr_candidate") == "yes"]
        row.append("%s nt; %s" % (max(int(x["aligned_nt"]) for _, x in tc),
                                  " / ".join(fmt(max([num(x["reads"]) for t_, x in tc if t_ == d.t] or [0]))
                                             for d in td)) if tc else "–")
        rows.append(row)
    tot_rows = [["all arrays", "", "", ""] + [fmt(sum(num(x["reads"]) for x in d.cr_arrays)) for d in td] + [
        "%d arrays" % len({x["array_id"] for d in td for x in getattr(d, "cr_tracr", [])
                           if x.get("tracr_candidate") == "yes"})]]
    ax = fig.add_axes([0.03, 0.49, 0.53, 0.38])
    table(ax, rows + tot_rows, ["array", "taxon (bin)", "Cas type", "rep."] +
          ["%s reads" % d.t.replace("sRNA-", "") for d in td] + ["tracrRNA"],
          col_widths=[0.08, 0.26, 0.12, 0.05] + [0.30 / len(td)] * len(td) + [0.19], fontsize=6.0, left=(0, 1, 2),
          title="top CRISPR arrays — reads (strand, share); tracrRNA: anti-repeat length; reads %s" % "/".join(
              d.t.replace("sRNA-", "") for d in td), rowscale=1.12)
    # spacer profile of the top array
    if top:
        a0 = top[0]
        ax = fig.add_axes([0.60, 0.56, 0.37, 0.28])
        n = len(td)
        w = 0.8 / n
        nsp = 0
        for i, d in enumerate(td):
            us = [u for u in d.cr_units if u["array_id"] == a0]
            sp = [u for u in us if u["kind"] == "spacer"]
            nsp = max(nsp, len(sp))
            ax.bar([j + (i - (n - 1) / 2) * w for j in range(len(sp))], [num(u["reads"]) for u in sp], width=w,
                   color=TYPE_COLORS[i], label=d.t)
        ax.set_yscale("symlog", linthresh=10)
        ax.set_xlabel("spacer (transcription order, 1 = leader side)")
        ax.set_ylabel("reads")
        ax.set_xticks(range(0, nsp, max(1, nsp // 15)))
        ax.set_xticklabels([str(x + 1) for x in range(0, nsp, max(1, nsp // 15))])
        r0 = info[a0]
        ax.set_title("reads per spacer — array %s (%s, type %s)" % (a0.split("_")[-1],
                     short(tax_tail(r0.get("classification", "")) if r0.get("classification") else r0["bin"], 26),
                     r0["cas_type"]), loc="left", fontsize=8)
        ax.legend(frameon=False)
        # read ends of the top array
        ax = fig.add_axes([0.06, 0.08, 0.40, 0.30])
        for i, d in enumerate(td):
            for end, ls in (("5p", "-"), ("3p", "--")):
                h = {int(r["offset"]): num(r["reads"]) for r in d.cr_ends if r["array_id"] == a0 and r["end"] == end}
                tot = sum(h.values())
                if tot:
                    xs = [x for x in range(-30, 51)]
                    ax.plot(xs, [100 * h.get(x, 0) / tot for x in xs], ls, color=TYPE_COLORS[i], lw=1.2,
                            label="%s %s end" % (d.t, end.replace("p", "′")))
        ax.axvline(0, color="#bbb", lw=0.8)
        ax.set_xlabel("offset (nt): 5′ end vs spacer start, 3′ end vs spacer end (+ = downstream)")
        ax.set_ylabel("% of reads")
        ax.set_title("read ends on spacers — array %s" % a0.split("_")[-1], loc="left")
        ax.legend(frameon=False, fontsize=6.5)
    # read lengths on arrays
    ax = fig.add_axes([0.55, 0.08, 0.42, 0.30])
    for i, d in enumerate(td):
        for cl, ls in (("with_spacer", "-"), ("repeat_only", "--")):
            h = {int(r["length"]): num(r["reads"]) for r in d.cr_lengths if r["class"] == cl}
            tot = sum(num(r["reads"]) for r in d.cr_lengths)
            if h:
                xs = sorted(h)
                ax.plot(xs, [100 * h[x] / tot for x in xs], ls, color=TYPE_COLORS[i], lw=1.2, marker="o", ms=2,
                        label="%s %s" % (d.t, cl.replace("_", " ")))
    ax.set_xlabel("read length (nt)")
    ax.set_ylabel("% of array reads")
    ax.set_title("length of reads on CRISPR arrays", loc="left")
    ax.legend(frameon=False, fontsize=6.5)
    pdf.savefig(fig)
    plt.close(fig)


def _steps(rows, aid):
    """run-length coverage rows of one array -> xs, plus, minus (step function, value from pos onwards)"""
    rr = [r for r in rows if r["array_id"] == aid]
    return [int(r["pos"]) for r in rr], [num(r["plus"]) for r in rr], [num(r["minus"]) for r in rr]


def draw_locus(fig, box, aid, td, title):
    """gene/array track + one mirrored coverage track (+ up, - down) per sRNA type"""
    x0, y0, w, h = box
    ref = next((d for d in td if any(r["array_id"] == aid for r in d.cr_cov)), None)
    if ref is None:
        return
    xs_all = [x for d in td for x in _steps(d.cr_cov, aid)[0]]
    lo, hi = min(xs_all), max(xs_all)
    n = len(td)
    th = h * 0.22                                   # annotation track height
    ch = (h - th) / n                               # coverage track height
    # annotation track
    ax = fig.add_axes([x0, y0 + h - th, w, th])
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1.2, 1.2)
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=8)
    ax.plot([lo, hi], [0, 0], color="#999", lw=0.6)
    for g in ref.cr_genes:
        if g["array_id"] != aid:
            continue
        gs, ge = int(g["start"]), int(g["end"])
        y = 0.45 if g["strand"] == "+" else -0.45
        dx = (ge - gs) if g["strand"] == "+" else -(ge - gs)
        xs = gs if g["strand"] == "+" else ge
        ax.add_patch(matplotlib.patches.FancyArrow(xs, y, dx, 0, width=0.45, head_width=0.7,
                                                   head_length=min(250, 0.25 * (ge - gs)), length_includes_head=True,
                                                   color="#E69F00", lw=0))
        ax.text((gs + ge) / 2, y + (0.55 if y > 0 else -0.85), short(g["name"].replace("CRISPR-associated ", ""), 18),
                ha="center", fontsize=5.5)
    for u in ref.cr_units:
        if u["array_id"] != aid:
            continue
        us, ue = int(u["start"]), int(u["end"])
        if u["kind"] == "repeat":
            ax.add_patch(matplotlib.patches.Rectangle((us, -0.35), ue - us, 0.7, color="#333333", lw=0))
        else:
            ax.add_patch(matplotlib.patches.Rectangle((us, -0.25), ue - us, 0.5, color="#56B4E9", lw=0))
    # tracrRNA candidates (anti-repeat, any fraction): green box on the strand of the candidate
    seen = set()
    for d in td:
        for r in getattr(d, "cr_tracr", []):
            if r["array_id"] != aid or r.get("tracr_candidate") != "yes":
                continue
            ts, te = int(r["start"]), int(r["end"])
            key = (ts // 20, r["strand"])
            if key in seen or te < lo or ts > hi:
                continue
            seen.add(key)
            y = 0.45 if r["strand"] == "+" else -0.45
            ax.add_patch(matplotlib.patches.Rectangle((ts, y - 0.3), max(te - ts, (hi - lo) * 0.004), 0.6,
                                                      color="#009E73", lw=0, zorder=5))
            ax.text((ts + te) / 2, y + (0.55 if y > 0 else -0.85),
                    "tracrRNA? (%s nt anti-repeat, p %s)" % (r["aligned_nt"], "%.2g" % num(r["p_shuffle"])),
                    ha="center", fontsize=5.5, color="#007A5A")
    # coverage tracks
    for i, d in enumerate(td):
        xs, pl, mi = _steps(d.cr_cov, aid)
        ax = fig.add_axes([x0, y0 + h - th - (i + 1) * ch, w, ch * 0.92])
        ax.set_xlim(lo, hi)
        if xs:
            ax.fill_between(xs, pl, step="post", color=TYPE_COLORS[i], lw=0)
            ax.fill_between(xs, [-v for v in mi], step="post", color=TYPE_COLORS[i], alpha=0.45, lw=0)
            m = max(max(pl), max(mi), 1)
            ax.set_ylim(-m * 1.05, m * 1.05)
            ax.set_yticks([-m, 0, m])
            ax.set_yticklabels(["−%s" % fmt(m), "0", "+%s" % fmt(m)], fontsize=6)
        ax.axhline(0, color="#999", lw=0.5)
        ax.set_ylabel(d.t.replace("sRNA-", ""), fontsize=7, rotation=0, ha="right", va="center")
        ax.tick_params(axis="x", labelsize=6)
        if i < n - 1:
            ax.set_xticklabels([])
        else:
            gl = [g for g in ref.cr_genes if g["array_id"] == aid]
            ax.set_xlabel("%s (bp; + strand up, − strand down; reads weighted 1/NH)" % gl[0]["contig"]
                          if gl else "position (bp)", fontsize=6.5)


def page_crispr_loci(pdf, sample, td):
    score = collections.defaultdict(float)
    info = {}
    for d in td:
        covered = {r["array_id"] for r in d.cr_cov}
        for r in d.cr_arrays:
            if r["array_id"] in covered:
                score[r["array_id"]] = max(score[r["array_id"]], num(r["reads"]))
                info.setdefault(r["array_id"], r)
    top = sorted(score, key=lambda k: -score[k])[:2]
    if not top:
        return False
    fig = page(pdf, sample, "CRISPR loci",
               "Cas genes (orange, above: + strand), array (dark: repeats, blue: spacers), tracrRNA candidates "
               "(green) and sRNA coverage per fraction; full locus left, array only right", explain="crispr")
    for k, aid in enumerate(top):
        r = info[aid]
        tl = "%s — %s, type %s, %s repeats" % (aid.split("_")[-1], tax_tail(r.get("classification", ""))
                                                if r.get("classification") else r["bin"], r["cas_type"], r["repeats"])
        yb = 0.50 - k * 0.46
        draw_locus(fig, (0.08, yb, 0.52, 0.34), aid, td, tl)
        # zoom on the array: reuse the locus drawer with coverage cut to the array
        sub = []
        for d in td:
            dd = type("D", (), {})()
            dd.t = d.t
            a0 = int(r["start"]) - 60
            a1 = int(r["end"]) + 60
            dd.cr_cov = [x for x in d.cr_cov if x["array_id"] == aid and a0 <= int(x["pos"]) <= a1]
            dd.cr_genes = []
            dd.cr_units = d.cr_units
            sub.append(dd)
        draw_locus(fig, (0.67, yb, 0.30, 0.34), aid, sub, "array %s only" % aid.split("_")[-1])
    pdf.savefig(fig)
    plt.close(fig)


SPECIES_MIN_DISTINCT = 100     # flag species with fewer distinct minimizers
SPECIES_MIN_DIRECT = 50        # ... or fewer reads assigned directly by Kraken2


def assignment_depth(kr):
    """share of rRNA reads by the rank Kraken2 assigned them to (taxon-level reads)"""
    cats = {"species": 0.0, "genus": 0.0, "family": 0.0, "order and above": 0.0, "unclassified": 0.0}
    for r in kr:
        rk = r["rank"]
        if rk == "U":
            cats["unclassified"] += r["taxon"]
        elif rk.startswith("S"):
            cats["species"] += r["taxon"]
        elif rk.startswith("G"):
            cats["genus"] += r["taxon"]
        elif rk.startswith("F"):
            cats["family"] += r["taxon"]
        else:
            cats["order and above"] += r["taxon"]
    return cats


def page_species(pdf, sample, td):
    if not any(d.bracken["S"] or d.kraken for d in td):
        return False
    fig = page(pdf, sample, "rRNA taxonomy at species level",
               "Kraken2 (GTDB) on RiboDetector rRNA reads → Bracken -l S on the filtered report; "
               "weak = < %d distinct minimizers or < %d reads assigned directly by Kraken2"
               % (SPECIES_MIN_DISTINCT, SPECIES_MIN_DIRECT), explain="species")
    depth_cols = ["species", "genus", "family", "order and above", "unclassified"]
    depth_colors = {"species": "#0072B2", "genus": "#56B4E9", "family": "#009E73", "order and above": "#E69F00",
                    "unclassified": "#cccccc"}
    # assignment depth: stacked bars + table
    ax = fig.add_axes([0.10, 0.62, 0.36, 0.22])
    deps = []
    for i, d in enumerate(td):
        dep = assignment_depth(d.kraken)
        deps.append(dep)
        tot = sum(dep.values())
        left = 0
        for c in depth_cols:
            w = dep[c] / tot if tot else 0
            ax.barh(i, w, left=left, color=depth_colors[c], label=c if i == 0 else None)
            left += w
    rows = [["rRNA reads"] + [fmt(sum(dp.values())) for dp in deps]]
    rows += [[c] + ["%s  (%s)" % (fmt(dp[c]), pct(dp[c], sum(dp.values()))) for dp in deps] for c in depth_cols]
    rows += [["Bracken reads at species level"] + [fmt(sum(num(r.get("new_est_reads")) for r in d.bracken["S"]))
                                                   for d in td]]
    ax.set_yticks(range(len(td)))
    ax.set_yticklabels([d.t for d in td])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of rRNA reads")
    ax.set_title("deepest rank assigned by Kraken2", loc="left")
    ax.legend(frameon=False, ncol=5, fontsize=6.5, bbox_to_anchor=(0, -0.28), loc="upper left")
    ax = fig.add_axes([0.56, 0.58, 0.41, 0.26])
    table(ax, rows, ["deepest rank"] + [d.t for d in td], col_widths=[0.40] + [0.30] * len(td), fontsize=6.8,
          title="reads per deepest Kraken2 rank")

    # top species per type with evidence
    n = len(td)
    for i, d in enumerate(td):
        ax = fig.add_axes([0.03 + i * (0.94 / n), 0.03, 0.94 / n - 0.02, 0.44])
        ev = {r["taxid"]: r for r in d.kraken if r["rank"].startswith("S")}
        top = sorted(d.bracken["S"], key=lambda r: -num(r.get("fraction_total_reads")))[:12]
        rows = []
        for r in top:
            k = ev.get(str(r.get("taxonomy_id", "")).strip(), {})
            dist = k.get("distinct")
            direct = num(r.get("kraken_assigned_reads"))
            weak = (dist is not None and dist < SPECIES_MIN_DISTINCT) or direct < SPECIES_MIN_DIRECT
            rows.append([short(r["name"], 34), "%.3f" % num(r["fraction_total_reads"]), fmt(num(r["new_est_reads"])),
                         fmt(direct), fmt(dist) if dist is not None else "–", "weak" if weak else ""])
        sp_reads = sum(num(r.get("new_est_reads")) for r in d.bracken["S"])
        table(ax, rows, ["species (GTDB)", "fraction", "Bracken reads", "Kraken direct", "distinct min.", "flag"],
              col_widths=[0.38, 0.11, 0.14, 0.13, 0.13, 0.11], fontsize=6.2, left=(0, 5),
              title="%s — top species (Bracken: %s reads at species level)" % (d.t, fmt(sp_reads)))
        if not rows:
            note(ax, "no Bracken species output")
    pdf.savefig(fig)
    plt.close(fig)


def page_assembly(pdf, sample, td):
    fig = page(pdf, sample, "fragment classes of tRNAs and rRNAs",
               "non-human reads mapped to the sample's own contigs; multi-mapping reads weighted 1/hits",
               explain="assembly")
    # categories
    dicts = [frac(d.category) for d in td]
    if any(dicts):
        labels = top_union(dicts, 10)
        ax = fig.add_axes([0.16, 0.54, 0.32, 0.31])
        grouped_barh(ax, labels, [(d.t, [dd.get(l, 0) for l in labels]) for d, dd in zip(td, dicts)],
                     TYPE_COLORS, "reads per feature category (fraction)")
    # fragment classes tRNA / rRNA
    for j, ft in enumerate(("tRNA", "rRNA")):
        dicts = [{r["class"]: num(r["fraction_of_type"]) for r in d.frag_classes if r["type"] == ft} for d in td]
        ax = fig.add_axes([0.62, 0.54 if j == 0 else 0.06, 0.35, 0.31 if j == 0 else 0.22])
        if not any(dicts):
            note(ax, "no %s fragment classes" % ft)
            continue
        order = ["full_length", "5p_long", "5p_half", "5p_tRF", "i_tRF", "3p_tRF", "3p_half", "3p_long",
                 "5p_end", "internal", "3p_end"]
        labels = [c for c in order if any(c in dd for dd in dicts)] + \
                 sorted({c for dd in dicts for c in dd} - set(order))
        grouped_barh(ax, labels, [(d.t, [dd.get(l, 0) for l in labels]) for d, dd in zip(td, dicts)],
                     TYPE_COLORS, "%s fragment classes (fraction of %s reads)" % (ft, ft))
    # length profiles
    ax = fig.add_axes([0.06, 0.06, 0.45, 0.36])
    styles = {"tRNA": "-", "rRNA": "--"}
    plotted = False
    for i, d in enumerate(td):
        for ft in ("tRNA", "rRNA"):
            h = {}
            for r in d.frag_lengths:
                if r["type"] == ft:
                    L = int(r["read_length"])
                    h[L] = h.get(L, 0) + num(r["weighted_reads"])
            tot = sum(h.values())
            if tot:
                xs = sorted(h)
                ax.plot(xs, [100 * h[x] / tot for x in xs], styles[ft], color=TYPE_COLORS[i], lw=1.2,
                        label="%s %s" % (d.t, ft))
                plotted = True
    if plotted:
        ax.set_xlabel("read length (nt)")
        ax.set_ylabel("% of reads of the type")
        ax.set_title("fragment length (tRNA solid, rRNA dashed)", loc="left")
        ax.legend(frameon=False, ncol=2)
    else:
        note(ax, "no fragment length table")
    pdf.savefig(fig)
    plt.close(fig)


INTEGRITY_ROWS = [
    ("reads_le20nt_fraction", "reads ≤20 nt", "↑ fragmented"),
    ("reads_ge30nt_fraction", "reads ≥30 nt", ""),
    ("median_read_length", "median read length (nt)", ""),
    ("assembly_mapped_fraction", "mapped to own contigs", ""),
    ("tRNA_mature_end_fraction", "tRNA reads with a mature end", "↑ intact"),
    ("tRNA_internal_fraction", "tRNA internal fragments (i-tRF)", "↑ degraded"),
    ("tRNA_3p_CCA_fraction", "3′-anchored tRNA reads with CCA", "↑ intact"),
    ("tRNA_halves_fraction", "tRNA halves", "stress, not degradation"),
    ("tRNA_5p_top3_share", "tRNA 5′ ends at main positions", "↑ sharp"),
    ("tRNA_3p_top3_share", "tRNA 3′ ends at main positions", "↑ sharp"),
    ("tRNA_5p_end_entropy", "tRNA 5′-end entropy", "↑ diffuse"),
    ("rRNA_5p_end_entropy", "rRNA 5′-end entropy", "↑ diffuse"),
    ("rRNA_bin_evenness", "rRNA spread along 16S/23S", "↑ random cuts"),
    ("rRNA_16S_5p_3p_log2ratio", "16S 5′/3′ reads (log2)", ""),
    ("rRNA_23S_5p_3p_log2ratio", "23S 5′/3′ reads (log2)", ""),
]


RNA_QC = {}          # (sample, type) -> {column: value}; from --rna-qc
RNA_QC_KEYS = ("sample", "sample_id", "type")


def load_rna_qc(path):
    """wet-lab QC table: sample, [sample_id], type, then any columns (RIN, DV200, Qubit, A260/280, ...)"""
    if not path:
        return
    if not os.path.exists(path):
        sys.stderr.write("[report] RNA QC table %s not found — ignored\n" % path)
        return
    for r in read_tsv(path):
        smp, t = r.get("sample", ""), r.get("type", "")
        if smp and t:
            RNA_QC[(smp, t)] = r
            if r.get("sample_id"):
                RNA_QC.setdefault((r["sample_id"], t), r)


def rna_qc_for(sample, t):
    return RNA_QC.get((sample, t)) or RNA_QC.get((sample.split(".")[0], t)) or {}


def page_integrity(pdf, sample, td):
    qcs = [rna_qc_for(sample, d.t) for d in td]
    if not any(d.integrity for d in td) and not any(qcs):
        return False
    fig = page(pdf, sample, "RNA integrity",
               "module 'fragments': end definition of tRNA/rRNA fragments on the assembly, read lengths; "
               "relative indicators, compare within a fraction", explain="integrity")
    vals = [{r["metric"]: r["value"] for r in d.integrity} for d in td]
    rows = []
    for key, label, direction in INTEGRITY_ROWS:
        if not any(key in v for v in vals):
            continue
        cells = []
        for v in vals:
            x = v.get(key, "")
            if x == "":
                cells.append("–")
            elif key == "median_read_length" or "log2" in key:
                cells.append("%.1f" % num(x) if "log2" in key else "%d" % num(x))
            else:
                cells.append("%.3f" % num(x))
        rows.append([label] + cells + [direction])
    qcols = []
    for q in qcs:
        for k in q:
            if k not in RNA_QC_KEYS and k not in qcols:
                qcols.append(k)
    qrows = []
    for k in qcols:
        cells = [(q.get(k) or "–").replace("NA", "–") for q in qcs]
        if all(c == "–" for c in cells):
            continue
        qrows.append(["wet lab: " + k.replace("_", " ")] + [short(c, 18) for c in cells] +
                     ["extraction QC" if "rin" not in k.lower() else "↑ intact (long RNA)"])
    rows = qrows + rows
    ax = fig.add_axes([0.03, 0.06, 0.46, 0.80])
    table(ax, rows, ["indicator"] + [d.t for d in td] + ["direction"],
          col_widths=[0.50] + [0.15] * len(td) + [0.26], left=(0, len(td) + 1), rowscale=1.45,
          title="integrity indicators")
    for j, k in enumerate(("16S", "23S")):
        ax = fig.add_axes([0.58, 0.50 - j * 0.42, 0.39, 0.32])
        plotted = False
        for i, d in enumerate(td):
            pr = [r for r in d.integrity_profile if r["rRNA"] == k]
            if pr:
                xs = [int(r["bin"]) for r in pr]
                ax.plot(xs, [100 * num(r["fraction"]) for r in pr], "-o", ms=3, color=TYPE_COLORS[i], lw=1.2,
                        label=d.t)
                plotted = True
        if not plotted:
            note(ax, "no complete %s genes with enough reads" % k)
            continue
        ax.axhline(10, color="#999", lw=0.8, ls=":")
        ax.set_xticks(range(1, 11))
        ax.set_xticklabels(["5′"] + [""] * 8 + ["3′"])
        ax.set_ylabel("% of fragments")
        ax.set_title("fragment midpoints along complete %s genes (dotted: even)" % k, loc="left")
        ax.legend(frameon=False)
    pdf.savefig(fig)
    plt.close(fig)


def rfam_category(fam):
    f = fam.lower()
    if "rrna" in f:
        return "rRNA"
    if f.startswith("trna") or f.startswith("mt-trna"):
        return "tRNA"
    if f.startswith("tmrna") or f == "tmrna":
        return "tmRNA"
    return "other Rfam"


def page_novel(pdf, sample, td):
    fig = page(pdf, sample, "novel ncRNA candidate loci",
               "read clusters on the contigs (18–250 nt, ≥ NOVEL_MINREADS weighted reads) outside annotated "
               "features (intergenic) or antisense to a CDS; best cmscan (Rfam) hit per locus; MFE from RNAfold",
               explain="novel")
    cats = ["rRNA", "tRNA", "tmRNA", "other Rfam", "no Rfam hit"]
    loci = []
    for d in td:
        seen = {}
        for b in ("intergenic", "antisense_to_CDS"):
            for c in d.cand[b]:
                key = "%s:%d-%d(%s)" % (c[0], c[1], c[2], c[3])
                if key in seen:
                    seen[key]["cls"] = "intergenic+antisense"
                    continue
                hit = d.rfam_tbl.get(key)
                seen[key] = {"key": key, "c": c, "cls": "intergenic" if b == "intergenic" else "antisense",
                             "hit": hit, "cat": rfam_category(hit[0]) if hit else "no Rfam hit",
                             "mfe": d.mfe.get(key)}
        loci.append(list(seen.values()))

    # summary per category
    rows = []
    for cat in cats:
        r = [cat]
        for L in loci:
            sel = [x for x in L if x["cat"] == cat]
            r.append("%s  /  %s" % (fmt(len(sel)), fmt(sum(x["c"][4] for x in sel))))
        rows.append(r)
    rows.append(["all loci"] + ["%s  /  %s" % (fmt(len(L)), fmt(sum(x["c"][4] for x in L))) for L in loci])
    rows.append(["  of which antisense to CDS"] + [fmt(sum(1 for x in L if "antisense" in x["cls"])) for L in loci])
    ax = fig.add_axes([0.03, 0.56, 0.40, 0.32])
    table(ax, rows, ["Rfam class of locus"] + ["%s: loci / reads" % d.t for d in td],
          col_widths=[0.36] + [0.32] * len(td), title="candidate loci by Rfam class (reads = weighted)")

    # read support per category
    ax = fig.add_axes([0.50, 0.60, 0.20, 0.25])
    bins = [10 ** (x / 4) for x in range(4, 26)]
    for i, (d, L) in enumerate(zip(td, loci)):
        v1 = [x["c"][4] for x in L if x["cat"] in ("rRNA", "tRNA", "tmRNA") and x["c"][4] > 0]
        v2 = [x["c"][4] for x in L if x["cat"] in ("other Rfam", "no Rfam hit") and x["c"][4] > 0]
        if v1:
            ax.hist(v1, bins=bins, histtype="step", lw=1.0, ls="--", color=TYPE_COLORS[i], label="%s rRNA/tRNA" % d.t)
        if v2:
            ax.hist(v2, bins=bins, histtype="step", lw=1.4, color=TYPE_COLORS[i], label="%s other" % d.t)
    ax.set_xscale("log")
    ax.set_xlabel("weighted reads per locus")
    ax.set_ylabel("loci")
    ax.set_title("read support", loc="left")
    ax.legend(frameon=False, fontsize=6)

    # other Rfam families
    dicts = []
    for L in loci:
        h = {}
        for x in L:
            if x["cat"] == "other Rfam":
                h[x["hit"][0]] = h.get(x["hit"][0], 0) + 1
        dicts.append(h)
    ax = fig.add_axes([0.85, 0.60, 0.12, 0.25])
    if any(dicts):
        labels = top_union(dicts, 10)
        grouped_barh(ax, labels, [(d.t, [dd.get(l, 0) for l in labels]) for d, dd in zip(td, dicts)],
                     TYPE_COLORS, "other Rfam families (loci)")
    else:
        note(ax, "no other Rfam families")

    # top other candidates (no rRNA/tRNA/tmRNA)
    n = len(td)
    for i, (d, L) in enumerate(zip(td, loci)):
        ax = fig.add_axes([0.03 + i * (0.94 / n), 0.03, 0.94 / n - 0.02, 0.47])
        other = sorted([x for x in L if x["cat"] in ("other Rfam", "no Rfam hit")], key=lambda x: -x["c"][4])
        rows = []
        for x in other[:17]:
            c = x["c"]
            rows.append([short(c[0].replace("_bridges_", "_b_").replace("_merged_", "_m_"), 24),
                         "%d-%d%s" % (c[1], c[2], c[3]), c[2] - c[1], x["cls"].replace("intergenic+antisense", "both"),
                         fmt(c[4]), short(x["hit"][0], 16) if x["hit"] else "—",
                         "%.1f" % x["mfe"] if x["mfe"] is not None else "–"])
        table(ax, rows, ["contig", "pos", "nt", "class", "reads", "Rfam", "MFE"],
              col_widths=[0.28, 0.16, 0.06, 0.13, 0.11, 0.16, 0.10], fontsize=6,
              title="%s — top candidates other than rRNA / tRNA / tmRNA" % d.t, left=(0, 3, 5))
    pdf.savefig(fig)
    plt.close(fig)


def page_references(pdf, sample):
    fig = page(pdf, sample, "marker sets and references")
    y = 0.88
    fig.text(0.03, y, "Human miRNA origin marker sets (source: %s)" % MIRNA_MARKER_SOURCE, fontsize=10,
             fontweight="bold", va="top")
    y -= 0.04
    refs, refno = all_references()
    for g, c, mirnas, rr in MIRNA_ORIGIN:
        fig.patches.append(matplotlib.patches.Rectangle((0.03, y - 0.017), 0.012, 0.016, transform=fig.transFigure,
                                                        color=c, figure=fig))
        txt = "%s%s: %s" % (g, " [%s]" % ",".join(str(refno[r]) for r in rr) if rr else "", ", ".join(mirnas))
        lines = textwrap.wrap(txt, 170)
        fig.text(0.05, y, "\n".join(lines), fontsize=8, va="top", linespacing=1.3)
        y -= 0.028 * len(lines) + 0.012
    for g, rr in MIRNA_GENERAL_REFS:
        fig.text(0.05, y, "%s: [%s]" % (g, ",".join(str(refno[r]) for r in rr)), fontsize=8, va="top", color="#444")
        y -= 0.04
    y -= 0.02
    fig.text(0.03, y, "References", fontsize=10, fontweight="bold", va="top")
    y -= 0.04
    for i, r in enumerate(refs):
        lines = textwrap.wrap("[%d]  %s" % (i + 1, r), 175)
        fig.text(0.03, y, "\n".join(lines), fontsize=8, va="top", linespacing=1.3)
        y -= 0.026 * len(lines) + 0.01
    y -= 0.02
    fig.text(0.03, y, "Notes", fontsize=10, fontweight="bold", va="top")
    notes = ("Matching is exact by mature miRNA name (miRBase, without 'hsa-'); each read counts for one group only, shares "
             "are relative to all human miRNA reads of the fraction. Markers are indicative, not absolute: the miR-200 family "
             "and miR-141 mark epithelia in general, miR-148a and miR-375 are also high in other tissues, miR-16 and "
             "miR-92a are abundant in blood but widely expressed, miR-146a and miR-155 also respond to inflammation outside "
             "immune cells. Low read numbers (tens to hundreds of miRNA reads) allow patterns, not quantitative comparisons.")
    fig.text(0.03, y - 0.035, "\n".join(textwrap.wrap(notes, 175)), fontsize=8, va="top", linespacing=1.35, color="#333")
    pdf.savefig(fig)
    plt.close(fig)


def page_bins(pdf, sample, base):
    cdir = os.path.join(base, "Analysis", "comparison")
    summ = read_tsv(os.path.join(cdir, "%s.compare.summary.tsv" % sample))
    allb = read_tsv(os.path.join(cdir, "%s.all_sources.bins.tsv" % sample))
    if not summ or not allb:
        return False
    srcs = [c[:-6] for c in allb[0] if c.endswith("_count")]
    fig = page(pdf, sample, "bins and multi-omics comparison",
               "CPM among binned reads; Spearman over HQ+MQ bins (≥ COMPARE_MIN_COUNT reads in both sources); "
               "bases: contig = all reads, cds = metaT/metaG CDS only, no_rrna = rRNA contigs left out",
               explain="bins")
    binned = [r for r in allb if r["bin"] != "unbinned"]
    tot = {s: sum(num(r[s + "_count"]) for r in binned) or 1 for s in srcs}
    allt = {s: sum(num(r[s + "_count"]) for r in allb) or 1 for s in srcs}
    unb = [r for r in allb if r["bin"] == "unbinned"]
    # Spearman table
    sp = {}
    for r in summ:
        if r.get("feature_type") == "HQ+MQ" and r.get("metric") == "spearman_cpm":
            sp.setdefault(r["pair"], {})[r["level"]] = r["value"]
    levels = ["bin", "bin_cds", "bin_no_rrna"]
    feat = {}
    for r in summ:
        if r.get("level") == "feature" and r.get("metric") == "spearman_cpm_within_type" and \
                r.get("feature_type") in ("CDS", "tRNA", "rRNA"):
            feat.setdefault(r["pair"], {})[r["feature_type"]] = r["value"]
    rows = []
    for pair in sorted(sp, key=lambda p: [s in p for s in ("metaG", "metaT")]):
        rows.append([pair.replace("_vs_", " vs ")] + [sp[pair].get(l, "–") for l in levels] +
                    [feat.get(pair, {}).get(t, "–") for t in ("CDS", "tRNA", "rRNA")])
    ax = fig.add_axes([0.03, 0.62, 0.50, 0.24])
    table(ax, rows, ["pair", "bins", "bins cds", "bins no_rRNA", "CDS feat.", "tRNA feat.", "rRNA feat."],
          col_widths=[0.30, 0.11, 0.11, 0.13, 0.11, 0.12, 0.12], fontsize=6.5,
          title="Spearman ρ — bins (HQ+MQ) and features within type")
    # unbinned share
    ax = fig.add_axes([0.66, 0.66, 0.30, 0.19])
    ax.barh(range(len(srcs)), [num(unb[0][s + "_count"]) / allt[s] if unb else 0 for s in srcs],
            color=[SOURCE_COLORS.get(s, "#888") for s in srcs])
    ax.set_yticks(range(len(srcs)))
    ax.set_yticklabels(srcs)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_title("share of reads on unbinned contigs", loc="left")
    # scatter plots
    good = [r for r in binned if r.get("quality") in ("HQ", "MQ")] or binned

    def cpm(r, s):
        return 1e6 * num(r[s + "_count"]) / tot[s]

    pairs = [p for p in (("sRNA-Ext", "sRNA-Int"), ("metaG", "sRNA-Ext"), ("metaG", "sRNA-Int"), ("metaG", "metaT"))
             if p[0] in srcs and p[1] in srcs]
    for k, (x, y) in enumerate(pairs[:4]):
        ax = fig.add_axes([0.05 + k * 0.245, 0.36, 0.19, 0.20])
        xs = [cpm(r, x) + 1 for r in good]
        ys = [cpm(r, y) + 1 for r in good]
        ax.scatter(xs, ys, s=8, alpha=0.7, color=SOURCE_COLORS.get(y, "#444"),
                   edgecolors="none")
        ax.set_xscale("log")
        ax.set_yscale("log")
        lo = min(xs + ys)
        hi = max(xs + ys)
        ax.plot([lo, hi], [lo, hi], color="#bbb", lw=0.8, zorder=0)
        ax.set_xlabel("%s CPM + 1" % x)
        ax.set_ylabel("%s CPM + 1" % y)
        ax.set_title("%s vs %s (ρ %s)" % (y, x, sp.get("%s_vs_%s" % tuple(sorted((x, y), key=_src_order)), {}).get("bin", "–")),
                     loc="left", fontsize=8)
    # top bins table
    ref = "sRNA-Ext" if "sRNA-Ext" in srcs else srcs[0]
    top = sorted(binned, key=lambda r: -num(r[ref + "_count"]))[:7]
    rows = [[r["bin"], short(tax_tail(r.get("classification", "")), 30), r.get("quality", "")] +
            [fmt(cpm(r, s)) for s in srcs] for r in top]
    ax = fig.add_axes([0.03, 0.03, 0.94, 0.22])
    table(ax, rows, ["bin", "taxon (GTDB)", "Q"] + ["%s CPM" % s for s in srcs],
          col_widths=[0.07, 0.26, 0.04] + [0.63 / len(srcs)] * len(srcs), fontsize=6,
          title="top bins by %s (CPM among binned reads)" % ref, left=(0, 1))
    pdf.savefig(fig)
    plt.close(fig)
    return True


def _src_order(s):
    return {"sRNA-Ext": 0, "sRNA-Int": 1, "metaT": 3, "metaG": 4}.get(s, 2)



# ------------------------------------------------------------------ Markdown pages (automated Methods)
SUP = str.maketrans("0123456789,–-", "⁰¹²³⁴⁵⁶⁷⁸⁹˒⁻⁻")


def md_inline_plain(t):
    """Markdown inline -> plain text for matplotlib (superscripts as Unicode, emphasis markers removed)"""
    t = re.sub(r"<sup>([^<]*)</sup>", lambda m: m.group(1).translate(SUP), t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"\1", t)
    t = re.sub(r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])", r"\1", t)
    t = re.sub(r"`([^`]*)`", r"\1", t)
    return re.sub(r"<!--.*?-->", "", t)


def md_inline_html(t):
    t = htmlmod.escape(t, quote=False)
    t = re.sub(r"&lt;sup&gt;(.*?)&lt;/sup&gt;", r"<sup>\1</sup>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])", r"<i>\1</i>", t)
    t = re.sub(r"`([^`]*)`", r"<code>\1</code>", t)
    t = re.sub(r"(https?://[^\s<)]+[^\s<).,;])", r'<a href="\1">\1</a>', t)
    return t


def md_blocks(text):
    """-> [(kind, payload)]: h1/h2/h3 (text), p (text), li (indent, marker, text), table (rows)"""
    out, para, rows = [], [], []

    def flush():
        if para:
            out.append(("p", " ".join(para)))
            para.clear()
        if rows:
            out.append(("table", list(rows)))
            rows.clear()
    for ln in text.splitlines():
        if ln.strip().startswith("<!--"):
            continue
        if not ln.strip():
            flush()
            continue
        m = re.match(r"(#{1,3})\s+(.*)", ln)
        if m:
            flush()
            out.append(("h%d" % len(m.group(1)), m.group(2)))
            continue
        if ln.startswith("|"):
            if para:
                flush()
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if not all(re.match(r"^-+$", c) for c in cells if c):
                rows.append(cells)
            continue
        m = re.match(r"(\s*)(\d+\.|[a-z]\.|[-*])\s+(.*)", ln)
        if m:
            flush()
            out.append(("li", (len(m.group(1)) // 2, m.group(2), m.group(3))))
            continue
        para.append(ln.strip())
    flush()
    return out


def md_to_html(text):
    h, inlist = [], False
    for kind, x in md_blocks(text):
        if kind == "li":
            lvl, mark, t = x
            h.append('<p class="li l%d"><span class="mk">%s</span> %s</p>' % (lvl, htmlmod.escape(
                mark if mark not in "-*" else "•"), md_inline_html(t)))
        elif kind == "p":
            h.append("<p>%s</p>" % md_inline_html(x))
        elif kind == "table":
            body = "".join("<tr>%s</tr>" % "".join(("<th>%s</th>" if i == 0 else "<td>%s</td>") % md_inline_html(c)
                                                    for c in r) for i, r in enumerate(x))
            h.append("<table>%s</table>" % body)
        else:
            h.append("<%s>%s</%s>" % (kind, md_inline_html(x), kind))
    return "\n".join(h)


def pages_markdown(pdf, sample, md_path, title="methods"):
    """render a Markdown file (headings, paragraphs, numbered steps, tables) onto as many A4 pages as needed,
    two columns; adds the text as real HTML to the HTML report. Returns False if the file is missing."""
    if not md_path or not os.path.exists(md_path):
        return False
    with open(md_path) as fh:
        text = fh.read()
    blocks = md_blocks(text)
    fs, lh = 7.0, 0.0182                     # font size, line height (figure fraction)
    cols = [(0.03, 0.485), (0.515, 0.97)]
    width_chars = 92                          # ~ characters per column line at 7 pt (DejaVu Sans)
    top, bottom = 0.89, 0.045
    lines = []                                # (text, style, indent) ; style: h1 h2 h3 p li th td gap
    for kind, x in blocks:
        if kind in ("h1", "h2", "h3"):
            lines += [("", "gap", 0), (md_inline_plain(x), kind, 0)]
        elif kind == "p":
            lines += [(w, "p", 0) for w in textwrap.wrap(md_inline_plain(x), width_chars)] + [("", "gap", 0)]
        elif kind == "li":
            lvl, mark, t = x
            ind = 3 + 3 * lvl
            ws = textwrap.wrap(md_inline_plain(t), width_chars - 2 * ind)
            mark = "•" if mark in "-*" else mark
            for i, w in enumerate(ws):
                lines.append(((mark + " " if i == 0 else "  ") + w, "li", lvl))
        elif kind == "table":
            ncol = max(len(r) for r in x)
            wcol = [40, 22, 36] if ncol == 3 else [width_chars // ncol] * ncol
            for ri, r in enumerate(x):
                cells = [textwrap.wrap(md_inline_plain(c), wcol[i]) or [""] for i, c in enumerate(r)]
                for k in range(max(len(c) for c in cells)):
                    lines.append(("\t".join((c[k] if k < len(c) else "") for c in cells),
                                  "th" if ri == 0 else ("tb" if r[0].startswith("**") else "td"), wcol))
            lines.append(("", "gap", 0))
    n_page, i = 0, 0
    per_col = int((top - bottom) / lh)
    html_done = False
    while i < len(lines):
        n_page += 1
        fig = page(pdf, sample, title + (" (%d)" % n_page if n_page > 1 else ""),
                   os.path.basename(md_path) + " — automated, check every [CHECK: …] before use")
        for x0, x1 in cols:
            y = top
            used = 0
            while i < len(lines) and used < per_col:
                t, st, ind = lines[i]
                i += 1
                if st == "gap":
                    if used:
                        y -= lh * 0.5
                        used += 0.5
                    continue
                if st in ("th", "td", "tb"):
                    xx = x0
                    tot = sum(ind)
                    for c, wc in zip(t.split("\t"), ind):
                        fig.text(xx, y, c, fontsize=fs - 0.6, va="top",
                                 fontweight="bold" if st in ("th", "tb") else "normal")
                        xx += (x1 - x0) * wc / tot
                else:
                    kw = dict(fontsize={"h1": 11, "h2": 9, "h3": 8}.get(st, fs), va="top",
                              fontweight="bold" if st.startswith("h") else "normal")
                    fig.text(x0 + (0.012 * (1 + ind) if st == "li" else 0), y, t, **kw)
                f = {"h1": 1.7, "h2": 1.5, "h3": 1.3}.get(st, 1.0)
                y -= lh * f
                used += f
        fig._no_html = True
        pdf.savefig(fig)
        plt.close(fig)
        if getattr(pdf, "html", False) and not html_done:
            pdf.pages.append(("%s — %s" % (sample, title), None, md_to_html(text)))
            html_done = True
    return True


class Sink:
    """receives the finished pages: PDF page and, optionally, a PNG for the HTML report"""

    def __init__(self, pdf, html=False, dpi=110):
        self.pdf, self.html, self.dpi, self.pages = pdf, html, dpi, []

    def savefig(self, fig):
        finalize(fig)
        self.pdf.savefig(fig)
        if self.html and not getattr(fig, "_no_html", False):
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=self.dpi)
            title = fig.texts[0].get_text() if fig.texts else "page %d" % (len(self.pages) + 1)
            self.pages.append((title, base64.b64encode(buf.getvalue()).decode("ascii")))


def write_html(path, sample, types, pages, pdf_name):
    esc = htmlmod.escape
    nav = " · ".join('<a href="#p%d">%s</a>' % (i, esc(p[0].split("—", 1)[-1].strip()))
                     for i, p in enumerate(pages))

    def sec(i, p):
        if len(p) > 2 and p[1] is None:            # text page (Markdown rendered to HTML)
            return '<section id="p%d"><h2>%s</h2><div class="md">%s</div></section>' % (i, esc(p[0]), p[2])
        return '<section id="p%d"><h2>%s</h2><img alt="%s" src="data:image/png;base64,%s"></section>' % (
            i, esc(p[0]), esc(p[0]), p[1])
    body = "\n".join(sec(i, p) for i, p in enumerate(pages))
    with open(path, "w") as out:
        out.write("""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s report %s</title>
<style>
 body{font-family:system-ui,-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f4f5f7;color:#222}
 header{position:sticky;top:0;background:#fff;border-bottom:1px solid #ddd;padding:10px 16px;z-index:1}
 header h1{font-size:18px;margin:0 0 4px} header .meta{font-size:12px;color:#666} nav{font-size:13px;margin-top:6px}
 nav a{color:#0072B2;text-decoration:none} nav a:hover{text-decoration:underline}
 section{max-width:1200px;margin:16px auto;padding:0 16px} section h2{font-size:14px;color:#555;margin:8px 0}
 img{width:100%%;height:auto;background:#fff;border:1px solid #ddd;border-radius:4px}
 .md{background:#fff;border:1px solid #ddd;border-radius:4px;padding:8px 24px;font-size:14px;line-height:1.5}
 .md h1{font-size:20px} .md h2{font-size:17px;color:#222} .md h3{font-size:15px}
 .md table{border-collapse:collapse;font-size:12.5px} .md td,.md th{border:1px solid #ddd;padding:3px 6px;text-align:left;vertical-align:top}
 .md .li{margin:2px 0 2px 1.2em} .md .l1{margin-left:2.6em} .md .mk{color:#555} .md code{background:#f3f3f3;padding:0 3px}
</style></head><body>
<header><h1>%s report — %s</h1>
<div class="meta">types: %s · generated %s · PDF: %s · <a href="https://%s">%s</a></div><nav>%s</nav></header>
%s
</body></html>
""" % (esc(PIPELINE_NAME), esc(sample), esc(PIPELINE_NAME), esc(sample), esc(", ".join(types)), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            esc(pdf_name), PIPELINE_URL, PIPELINE_URL, nav, body))


# ----------------------------------------------------------------------------- MultiQC
def _qc(sample, td, t):
    try:
        return qc_verdicts(sample, td)[t][1]
    except Exception:
        return ""


def _rin(sample, t):
    q = rna_qc_for(sample, t)
    for k in q:
        if "rin" in k.lower() and q[k] not in ("", "NA"):
            return q[k]
    return ""


def write_mqc(mqc_dir, sample, td):
    os.makedirs(mqc_dir, exist_ok=True)
    gs = os.path.join(mqc_dir, "%s_srna_general_stats_mqc.tsv" % sample)
    with open(gs, "w") as out:
        # one id for all samples: otherwise MultiQC makes a separate column set per file (179 x columns)
        out.write("# id: 'srna_general_stats'\n# section_name: 'sRNA-IMP'\n# plot_type: 'generalstats'\n# pconfig:\n")
        out.write("Sample\tsrna_raw\tsrna_trimmed_pct\tsrna_human_pct\tsrna_nonhuman_rRNA_pct\tsrna_tRNA_pct\t"
                  "srna_mapped_assembly_pct\tsrna_novel_loci\tsrna_lRNA_RIN\tsrna_qc\n")
        for d in td:
            raw = d.rs("raw") or 1
            nh = d.rs("nonhuman") or 1
            out.write("%s_%s\t%d\t%.2f\t%.3f\t%.2f\t%.2f\t%.2f\t%d" % (
                sample, d.t, raw, 100 * d.rs("trimmed") / raw, 100 * d.rs("human_total") / (d.rs("trimmed") or 1),
                100 * num(d.nonhuman.get("rRNA_reads")) / nh, 100 * num(d.nonhuman.get("tRNA")) / nh,
                100 * num(d.novel_sum.get("mapped_reads")) / (num(d.novel_sum.get("input_reads")) or 1),
                len(d.cand["intergenic"]) + len(d.cand["antisense_to_CDS"])) +
                      "\t%s\t%s\n" % (_rin(sample, d.t), _qc(sample, td, d.t)))
    bg = os.path.join(mqc_dir, "%s_srna_read_fate_mqc.tsv" % sample)
    with open(bg, "w") as out:
        out.write("# id: 'srna_read_fate'\n# section_name: 'sRNA read fate'\n# plot_type: 'bargraph'\n"
                  "# pconfig:\n#    id: 'srna_read_fate_plot'\n#    title: 'sRNA read fate'\n")
        cols = ["removed_in_trimming", "human", "rRNA", "tRNA", "Rfam", "miRNA_like", "unknown"]
        out.write("Sample\t" + "\t".join(cols) + "\n")
        for d in td:
            v = [max(d.rs("raw") - d.rs("trimmed"), 0), d.rs("human_total"),
                 num(d.nonhuman.get("rRNA_reads")), num(d.nonhuman.get("tRNA")), num(d.nonhuman.get("Rfam_classified")),
                 num(d.nonhuman.get("miRNA_like")), num(d.nonhuman.get("unknown"))]
            out.write("%s_%s\t" % (sample, d.t) + "\t".join("%.0f" % x for x in v) + "\n")
    bg = os.path.join(mqc_dir, "%s_srna_bracken_genus_mqc.tsv" % sample)
    with open(bg, "w") as out:
        out.write("# id: 'srna_bracken_genus'\n# section_name: 'rRNA genera (Bracken)'\n# plot_type: 'bargraph'\n"
                  "# pconfig:\n#    id: 'srna_bracken_genus_plot'\n#    title: 'Bracken genus (rRNA reads)'\n"
                  "#    cpswitch_c_active: false\n")
        dicts = [{r["name"]: num(r["new_est_reads"]) for r in d.bracken["G"]} for d in td]
        labels = top_union(dicts, 20)
        out.write("Sample\t" + "\t".join(labels) + "\tother\n")
        for d, dd in zip(td, dicts):
            other = sum(dd.values()) - sum(dd.get(l, 0) for l in labels)
            out.write("%s_%s\t" % (sample, d.t) + "\t".join("%.0f" % dd.get(l, 0) for l in labels) +
                      "\t%.0f\n" % other)
    ig = os.path.join(mqc_dir, "%s_srna_integrity_mqc.tsv" % sample)
    keys = [k for k, _, _ in INTEGRITY_ROWS]
    with open(ig, "w") as out:
        out.write("# id: 'srna_integrity'\n# section_name: 'sRNA integrity indicators'\n# plot_type: 'table'\n")
        qk = []
        for d in td:
            for k in rna_qc_for(sample, d.t):
                if k not in RNA_QC_KEYS and k not in qk and not "note" in k.lower():
                    qk.append(k)
        out.write("Sample\t" + "\t".join(qk + keys) + "\n")
        for d in td:
            v = {r["metric"]: r["value"] for r in d.integrity}
            q = rna_qc_for(sample, d.t)
            # NA instead of empty cells: MultiQC drops trailing empty fields ("inconsistent number of columns")
            out.write("%s_%s\t" % (sample, d.t) + "\t".join([q.get(k, "") or "NA" for k in qk] +
                                                          [v.get(k, "") or "NA" for k in keys]) + "\n")
    return mqc_dir


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="pipeline output directory (-o of srna_pipeline.sh)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--types", default="", help="comma-separated sRNA types (default: all found)")
    ap.add_argument("--pdf", default="", help="output PDF (default <outdir>/<sample>/Report/<sample>.report.pdf)")
    ap.add_argument("--rfam-family", default="", help="Rfam family.txt(.gz) for element types (ncRNA page)")
    ap.add_argument("--mirna-markers", default="",
                    help="TSV of human miRNA origin markers (default: mirna_origin_markers.tsv next to this script)")
    ap.add_argument("--rna-qc", default="",
                    help="wet-lab RNA QC TSV (sample, sample_id, type, RIN, ...), shown on the integrity page")
    ap.add_argument("--qc-rules", default="", help="QC rules TSV (default: qc_rules.tsv next to this script)")
    ap.add_argument("--qc-exclude", default="", help="manual QC decisions TSV (sample, type, level, reason)")
    ap.add_argument("--html", default="", help="also write a self-contained HTML report (same pages)")
    ap.add_argument("--methods", default="", help="Markdown file with the automated Methods (srna_methods.py), "
                    "added as the last pages")
    ap.add_argument("--mqc-dir", default="", help="also write MultiQC custom-content tables here")
    a = ap.parse_args()

    base = os.path.join(a.outdir, a.sample)
    if not os.path.isdir(base):
        sys.exit("[report] no such sample directory: %s" % base)
    if a.types:
        types = [t for t in a.types.split(",") if t]
    else:
        types = sorted({os.path.basename(p) for p in glob.glob(os.path.join(base, "Stats", "*")) if os.path.isdir(p)} |
                       {os.path.basename(p) for p in glob.glob(os.path.join(base, "Analysis", "*"))
                        if os.path.isdir(p) and os.path.basename(p) not in ("annotation", "comparison")},
                       key=_src_order)
    if not types:
        sys.exit("[report] %s: no sRNA types found" % a.sample)
    load_markers(a.mirna_markers or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "mirna_origin_markers.tsv"))
    load_rfam_family(a.rfam_family)
    load_rna_qc(a.rna_qc)
    import srna_qc as Q
    QC_STATE.update(rules=Q.load_rules(a.qc_rules), exclude=Q.load_exclude(a.qc_exclude), outdir=a.outdir,
                    rules_path=a.qc_rules, exclude_path=a.qc_exclude)
    td = [TypeData(base, a.sample, t) for t in types]
    pdf_path = a.pdf or os.path.join(base, "Report", "%s.report.pdf" % a.sample)
    os.makedirs(os.path.dirname(os.path.abspath(pdf_path)), exist_ok=True)

    done = []
    with PdfPages(pdf_path) as pdf_file:
        pdf = Sink(pdf_file, html=bool(a.html))
        for name, fn in [("overview", lambda: page_overview(pdf, a.sample, td, a.outdir)),
                         ("qc_verdict", lambda: page_qc_verdict(pdf, a.sample, td)),
                         ("qc", lambda: page_qc(pdf, a.sample, td)),
                         ("human", lambda: page_human(pdf, a.sample, td)),
                         ("nonhuman", lambda: page_nonhuman(pdf, a.sample, td)),
                         ("species", lambda: page_species(pdf, a.sample, td)),
                         ("assembly", lambda: page_assembly(pdf, a.sample, td)),
                         ("integrity", lambda: page_integrity(pdf, a.sample, td)),
                         ("assembly_cmc", lambda: page_asm_cmc(pdf, a.sample, td)),
                         ("ncrna", lambda: page_ncrna(pdf, a.sample, td, base)),
                         ("crispr", lambda: page_crispr(pdf, a.sample, td)),
                         ("crispr_loci", lambda: page_crispr_loci(pdf, a.sample, td)),
                         ("novel", lambda: page_novel(pdf, a.sample, td)),
                         ("bins", lambda: page_bins(pdf, a.sample, base)),
                         ("references", lambda: page_references(pdf, a.sample)),
                         ("methods", lambda: pages_markdown(pdf, a.sample, a.methods))]:
            try:
                r = fn()
                if r is not False:
                    done.append(name)
            except Exception as e:                       # one broken table must not kill the report
                sys.stderr.write("[report] %s: page %s skipped (%s: %s)\n" % (a.sample, name, type(e).__name__, e))
                plt.close("all")
        meta = pdf_file.infodict()
        meta["Title"] = "%s report — %s" % (PIPELINE_NAME, a.sample)
        meta["Subject"] = ", ".join(types)
    sys.stderr.write("[report] %s: %s (pages: %s)\n" % (a.sample, pdf_path, ", ".join(done)))
    if a.html:
        os.makedirs(os.path.dirname(os.path.abspath(a.html)), exist_ok=True)
        write_html(a.html, a.sample, types, pdf.pages, os.path.basename(pdf_path))
        sys.stderr.write("[report] %s: %s\n" % (a.sample, a.html))
    try:
        write_qc_flags(os.path.join(os.path.dirname(os.path.abspath(pdf_path)), "%s.qc_flags.tsv" % a.sample),
                       a.sample, td)
    except Exception as e:
        sys.stderr.write("[report] %s: qc_flags.tsv not written (%s)\n" % (a.sample, e))
    if a.mqc_dir:
        write_mqc(a.mqc_dir, a.sample, td)
        sys.stderr.write("[report] %s: MultiQC tables in %s\n" % (a.sample, a.mqc_dir))


if __name__ == "__main__":
    main()
