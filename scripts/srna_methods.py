#!/usr/bin/env python3
"""srna_methods.py — automated Methods text of the sRNA-IMP pipeline (per-sample processing and cohort analysis)
with tools, versions, parameters and references.

Styles
  nature  Nature-style Methods: run-in subsection headings, numbered references (superscript, order of first
          citation), 'Data availability' / 'Code availability' placeholders.
  star    STAR Protocols (Cell Press): Key resources table, step-by-step method details (numbered), quantification
          and statistical analysis, limitations, resource availability; numbered Cell Press references.
          ('cell' is accepted as an alias.)

Inputs
  --params FILE      NAME<TAB>VALUE of the resolved pipeline configuration (written by srna_pipeline.sh)
  --envs DIRS        comma-separated conda/micromamba environments: tool versions are read from conda-meta/
  --resources FILE   databases / upstream workflow (name, version, identifier, references); default
                     scripts/methods_resources.tsv — edit versions there once for your installation
  --outdir / --cohort-dir   pipeline output and cohort directory (numbers of samples, libraries, groups)
  --mode sample|cohort|all  sample = processing only (per-sample report); cohort / all = processing + cohort

Output: Markdown (--out); with --docx and pandoc in PATH also a .docx. Values that could not be determined are
written as [CHECK: …] so they are easy to find before submission.
"""
import argparse
import collections
import csv
import datetime
import glob
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srna_methods_refs import REFS  # noqa: E402

PIPELINE = "sRNA-IMP pipeline"
REPO_URL = "https://github.com/genomeanalysislcsb/sRNA-IMP"
DOCS_URL = "https://genomeanalysislcsb.github.io/sRNA-IMP/"
HERE = os.path.dirname(os.path.abspath(__file__))

# tool -> (display name, conda package, reference keys, URL)
TOOLS = collections.OrderedDict([
    ("cutadapt", ("Cutadapt", "cutadapt", "cutadapt", "https://github.com/marcelm/cutadapt")),
    ("fastqc", ("FastQC", "fastqc", "fastqc", "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/")),
    ("multiqc", ("MultiQC", "multiqc", "multiqc", "https://multiqc.info")),
    ("bowtie", ("Bowtie", "bowtie", "bowtie", "https://bowtie-bio.sourceforge.net")),
    ("samtools", ("SAMtools", "samtools", "samtools", "https://www.htslib.org")),
    ("bedtools", ("BEDTools", "bedtools", "bedtools", "https://bedtools.readthedocs.io")),
    ("bedops", ("BEDOPS", "bedops", "bedops", "https://bedops.readthedocs.io")),
    ("seqkit", ("SeqKit", "seqkit", "seqkit", "https://bioinf.shenwei.me/seqkit/")),
    ("ribodetector", ("RiboDetector", "ribodetector", "ribodetector", "https://github.com/hzi-bifo/RiboDetector")),
    ("kraken2", ("Kraken 2", "kraken2", "kraken2", "https://github.com/DerrickWood/kraken2")),
    ("bracken", ("Bracken", "bracken", "bracken", "https://github.com/jenniferlu717/Bracken")),
    ("infernal", ("Infernal (cmscan)", "infernal", "infernal", "http://eddylab.org/infernal/")),
    ("aragorn", ("ARAGORN", "aragorn", "aragorn", "http://www.ansikte.se/ARAGORN/")),
    ("trnascan", ("tRNAscan-SE", "trnascan-se", "trnascan", "http://lowelab.ucsc.edu/tRNAscan-SE/")),
    ("viennarna", ("ViennaRNA (RNAfold)", "viennarna", "viennarna", "https://www.tbi.univie.ac.at/RNA/")),
    ("python", ("Python", "python", "", "https://www.python.org")),
    ("numpy", ("NumPy", "numpy", "numpy", "https://numpy.org")),
    ("matplotlib", ("Matplotlib", "matplotlib", "matplotlib", "https://matplotlib.org")),
])


# ----------------------------------------------------------------------------------------------- inputs
def read_params(path):
    p = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for ln in fh:
                if "\t" in ln:
                    k, v = ln.rstrip("\n").split("\t", 1)
                    p[k] = v
    return p


def conda_versions(envs):
    """package -> version from conda-meta/*.json of the given environments (first hit wins)"""
    out = {}
    for e in envs:
        for f in sorted(glob.glob(os.path.join(e, "conda-meta", "*.json"))):
            m = re.match(r"(.+)-([^-]+)-([^-]+)\.json$", os.path.basename(f))
            if m and m.group(1) not in out:
                out[m.group(1)] = m.group(2)
    return out


def yml_versions():
    """fallback: pinned versions in envs/*.yml next to the scripts directory"""
    out = {}
    for f in glob.glob(os.path.join(HERE, "..", "envs", "*.yml")):
        with open(f) as fh:
            for ln in fh:
                m = re.match(r"\s*-\s*([A-Za-z0-9_.\-]+)\s*=\s*([0-9][^=\s#]*)", ln)
                if m and m.group(1) not in out:
                    out[m.group(1)] = m.group(2)
    return out


def read_resources(path):
    res = collections.OrderedDict()
    if path and os.path.exists(path):
        with open(path) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"):
                    continue
                c = ln.rstrip("\n").split("\t") + [""] * 6
                res[c[0]] = {"name": c[1], "version": c[2], "source": c[3], "id": c[4],
                             "refs": [x for x in c[5].split(",") if x]}
    return res


def git_version():
    try:
        v = subprocess.run(["git", "-C", HERE, "describe", "--always", "--dirty", "--tags"], capture_output=True,
                           text=True, timeout=10).stdout.strip()
        return v or None
    except Exception:
        return None


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ----------------------------------------------------------------------------------------------- citations
class Cites:
    def __init__(self):
        self.order = []

    def __call__(self, *keys):
        nums = []
        for k in keys:
            if not k or k not in REFS:
                continue
            if k not in self.order:
                self.order.append(k)
            nums.append(self.order.index(k) + 1)
        if not nums:
            return ""
        return "<sup>" + compress(sorted(set(nums))) + "</sup>"


def compress(nums):
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if j == i else ("%d,%d" % (nums[i], nums[j]) if j == i + 1 else
                                               "%d–%d" % (nums[i], nums[j])))
        i = j + 1
    return ",".join(out)


def fmt_ref(k, style):
    au, year, title, jour, vol, pages, doi = REFS[k]
    etal = au and au[-1] == "et al."
    au = [a for a in au if a != "et al."]
    pages = pages.replace("-", "–")
    url = doi if doi.startswith("http") else ("https://doi.org/" + doi if doi else "")
    if style == "nature":
        if etal or len(au) > 5:
            a = au[0] + " et al."
        elif len(au) == 1:
            a = au[0]
        else:
            a = ", ".join(au[:-1]) + " & " + au[-1]
        if not jour:
            return "%s %s. %s (%d)." % (a, title, url, year)
        loc = " ".join(x for x in ("*%s*" % jour, "**%s**," % vol if vol else "", pages) if x)
        return "%s %s. %s (%d)." % (a, title, loc.rstrip(","), year) + (" %s" % url if url else "")
    # Cell Press (STAR Protocols): up to 10 authors, then 'et al.'
    if etal or len(au) > 10:
        a = ", ".join(au[:10]) + ", et al."
    elif len(au) == 1:
        a = au[0]
    else:
        a = ", ".join(au[:-1]) + ", and " + au[-1]
    if not jour:
        return "%s (%d). %s. %s." % (a, year, title, url)
    loc = "*%s*" % jour + (" *%s*" % vol if vol else "") + (", %s" % pages if pages else "")
    return "%s (%d). %s. %s." % (a, year, title, loc) + (" %s." % url if url else "")


# ----------------------------------------------------------------------------------------------- context
class Ctx:
    def __init__(self, a):
        self.a = a
        self.p = read_params(a.params)
        envs = [e for e in (a.envs or "").split(",") if e]
        self.ver = yml_versions()
        self.ver.update(conda_versions(envs))           # installed versions override pinned ones
        self.res = read_resources(a.resources)
        self.c = Cites()
        self.git = git_version()
        self.used_tools, self.used_res = [], []

    def P(self, k, default=""):
        v = self.p.get(k, "")
        return v if v != "" else default

    def on(self, k, default="1"):
        return self.P(k, default) not in ("0", "", "no", "false")

    def v(self, tool):
        pkg = TOOLS[tool][1]
        return self.ver.get(pkg) or self.ver.get(pkg.lower()) or "[CHECK: version]"

    def tool(self, tool, with_version=True):
        self.used_tools.append(tool) if tool not in self.used_tools else None
        name, _, refs, _ = TOOLS[tool]
        return "%s%s%s" % (name, " v%s" % self.v(tool) if with_version else "",
                           self.c(*refs.split(",")) if refs else "")

    def r(self, key, what="name"):
        x = self.res.get(key)
        if not x:
            return "[CHECK: %s]" % key
        val = x.get(what, "")
        return val if val and not val.startswith("set ") else "[CHECK: %s %s]" % (key, what)

    def rc(self, key):
        """resource name (+ version) with citation"""
        self.used_res.append(key) if key not in self.used_res else None
        x = self.res.get(key)
        if not x:
            return "[CHECK: %s]" % key
        ver = x["version"] if x["version"] and not x["version"].startswith("set ") else "[CHECK: version]"
        name = x["name"] if not x["name"].startswith("set ") else "[CHECK: %s name]" % key
        return "%s%s%s" % (name, "" if ver in name else " (%s)" % ver, self.c(*x["refs"]))


# ----------------------------------------------------------------------------------------------- content
# Each section: (title, [items]). Nature joins the items of a section into one paragraph; STAR Protocols numbers
# them as steps. Items are complete sentences.
def qc_rules_text(X):
    path = X.P("QC_RULES") or os.path.join(HERE, "qc_rules.tsv")
    fail, warn = [], []
    if os.path.exists(path):
        with open(path) as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
        for r in rows[1:]:
            if len(r) >= 7:
                d = r[6] + ("" if r[5] in ("*", "") else " (%s only)" % r[5])
                (fail if r[4] == "FAIL" else warn).append(d)
    return fail, warn


def processing_sections(X, types):
    c, P = X.c, X.P
    S = []
    norc = X.on("DB_NORC")
    strand = (" Because the libraries are strand-specific, only the forward strand of the RNA reference "
              "databases was searched (bowtie --norc)." if norc else "")
    ver = " (%s)" % X.git if X.git else ""
    S.append(("Small RNA data processing overview", [
        "Small RNA sequencing libraries (%s) were processed per sample with the %s%s, a Bash workflow with Python "
        "helper scripts that runs locally or on a SLURM cluster." % (", ".join(types) if types else
                                                                      "[CHECK: fractions]", PIPELINE, ver),
        "Reads were trimmed, separated into human and non-human reads, classified into RNA classes, and mapped to "
        "the co-assembly of the metagenomic and metatranscriptomic reads of the same sample, reconstructed with "
        "the integrated multi-omic workflow %s, so that small RNAs could be placed on the sample's own contigs, genes, genomic bins and functional "
        "annotation." % X.rc("imp"),
    ]))
    adapter = P("ADAPTER_3P", "[CHECK: adapter]")
    tr = P("ADAPTER_3P_TRUNC", "0")
    items = [
        "Read quality was assessed before and after trimming with %s, and all QC outputs were summarised with "
        "%s." % (X.tool("fastqc"), X.tool("multiqc")),
        "The 3′ adapter (%s) was removed with %s%s, together with 3′ bases below Phred quality %s (-q %s)%s; "
        "reads of %s–%s nt were retained (-m %s -M %s)." % (
            adapter, X.tool("cutadapt"),
            (" (a second copy of the adapter lacking its first %s bases was trimmed as well, to remove adapter "
             "dimers whose reads start inside the adapter)" % tr) if tr not in ("", "0") else "",
            P("TRIM_QUALITY", "15"), P("TRIM_QUALITY", "15"),
            (" and with the options %s" % P("CUTADAPT_EXTRA")) if P("CUTADAPT_EXTRA") else "",
            P("MINLEN", "18"), P("MAXLEN", "100"), P("MINLEN", "18"), P("MAXLEN", "100")),
        "Adapter dimers, library complexity (distinct and singleton sequences) and the most frequent sequences were "
        "recorded per library.",
    ]
    if P("CONTAM_FA"):
        items.append("Technical contaminants were quantified by aligning the trimmed reads with %s (-v %s -k 1) to "
                     "%s." % (X.tool("bowtie"), P("CONTAM_MISMATCHES", "1"), X.rc("contaminants")))
    fail, warn = qc_rules_text(X)
    if fail or warn:
        items.append("Each library received a rule-based QC verdict: FAIL if %s; WARN if %s." % (
            "; ".join(fail) or "none", "; ".join(warn) or "none"))
    S.append(("Read preprocessing and quality control", items))

    S.append(("Separation of human and microbial reads", [
        "Reads shorter than %s nt were never assigned to the host." % P("HOST_MIN_LEN", "20"),
        "All other reads were aligned end-to-end with %s (-k 1 --best) to the human genome and, if unmapped, to "
        "the human transcriptome (exon junctions) of %s, %s." % (
            X.tool("bowtie"), X.rc("host_genome"),
            ("allowing %s mismatch(es)" % P("HOST_MISMATCHES", "0") if P("HOST_MISMATCHES", "0") != "0"
             else "requiring exact matches")
            if P("HOST_SHORT_MISMATCHES", "0") == P("HOST_MISMATCHES", "0") else
            "allowing %s mismatches for reads shorter than %s nt and %s for longer reads" % (
                P("HOST_SHORT_MISMATCHES", "0"), P("HOST_SHORT_LEN", "22"), P("HOST_MISMATCHES", "0"))),
        ("Reads with a human hit that matched the sample's own contigs with no more mismatches than the human "
         "reference were reassigned to the microbial fraction, so that conserved short sequences were not "
         "attributed to the host." if X.on("HOST_ASSEMBLY_RESCUE") else ""),
        "Human reads were de-duplicated by read name with %s." % X.tool("seqkit"),
    ]))
    S.append(("Human small RNA classes", [
        "Human reads were classified hierarchically, each read being assigned to the first matching class, with "
        "%s (-v 1 -a --best): mature miRNAs (%s), nuclear and mitochondrial tRNAs (%s), rRNAs (%s), other "
        "ncRNAs (%s), and unknown." % (X.tool("bowtie", False), X.rc("hsa_mirna"), X.rc("hsa_trna"),
                                        X.rc("hsa_rrna"), X.rc("hsa_other")),
        "A read with *n* equally good hits contributed 1/*n* to each reference sequence.",
    ]))
    S.append(("Microbial small RNA classes and rRNA-based taxonomy", [
        "Ribosomal RNA reads were identified among the non-human reads with %s (--len set to the mean read "
        "length)." % X.tool("ribodetector"),
        "The rRNA reads were classified taxonomically with %s (--confidence %s, --minimum-hit-groups %s%s) "
        "against %s; report entries with fewer than %s clade reads%s were removed, and abundances were "
        "re-estimated with %s (-r %s -t %s) at the ranks %s." % (
            X.tool("kraken2"), P("KRAKEN_CONF", "0.5"), P("KRAKEN_MINHITGROUPS", "3"),
            ", --report-minimizer-data" if X.on("KRAKEN_MINIMIZER_DATA") else "", X.rc("kraken_db"),
            P("KRAKEN_REPORT_MINREADS", "10"),
            (" and species with fewer than %s distinct minimizers" % P("KRAKEN_MIN_DISTINCT_MINIMIZERS"))
            if P("KRAKEN_MIN_DISTINCT_MINIMIZERS", "0") != "0" else "",
            X.tool("bracken"), P("BRACKEN_READ_LENGTH", "50"), P("BRACKEN_MIN_READS", "10"),
            ", ".join(P("BRACKEN_LEVELS", "D P C O F G S").split())),
        "The remaining reads were collapsed into unique sequences and classified hierarchically with %s: tRNAs "
        "(-v 1 -k 1 --best; %s), then sequences of %s–%s nt as miRNA-like (exact matches, -v 0, to %s), then the "
        "rest against Rfam sequences grouped by RNA class (-v 1 -k 1; %s); unmatched sequences were classed as "
        "unknown.%s" % (X.tool("bowtie", False), X.rc("trna_all"), P("MIRNA_LIKE_MINLEN", "18"),
                        P("MIRNA_LIKE_MAXLEN", "26"), X.rc("mirna_euk"), X.rc("rfam"), strand),
    ]))
    items = [
        "Contigs, gene and ncRNA annotation (%s) and functional annotation (%s; %s) of each sample were taken "
        "from the integrated multi-omic analysis; the GFF annotation was converted to BED with %s." % (
            X.rc("annotation"), X.rc("functional"), "KEGG%s, Pfam%s and COG%s" % (c("kegg"), c("pfam"), c("cog")),
            X.tool("bedops")),
        "All non-human reads (not collapsed) were aligned to the contigs with %s (-v %s --best --strata, %s) and "
        "processed with %s; the number of hits *n* of a read was stored as NH tag and every alignment was "
        "weighted 1/*n* in all downstream counts." % (
            X.tool("bowtie", False), P("ASSEMBLY_MISMATCHES", "1"),
            ("at most %s equally good hits, -k %s" % (P("ASSEMBLY_MAX_HITS"), P("ASSEMBLY_MAX_HITS")))
            if P("ASSEMBLY_MAX_HITS", "0") != "0" else "all equally good hits, -a", X.tool("samtools")),
        "Each alignment was assigned to the annotated feature it overlapped best (≥ %s %% of the read) on the same "
        "strand (sense) or, if there was none, on the opposite strand (antisense); otherwise it was counted as "
        "intergenic." % int(round(float(P("BINS_MIN_OVERLAP", "0.5")) * 100)),
        "Reads on protein-coding genes were summed per functional term (KEGG orthology, pathway and module, Pfam "
        "domain, COG), every term of a gene receiving the gene's reads.",
    ]
    rr = "the rRNA reads identified by RiboDetector"
    if X.on("RRNA_EXCLUDE_RFAM"):
        rr += (" and, because RiboDetector misses many short rRNA fragments, all reads whose sequence had an Rfam "
               "rRNA hit")
    items.append("Reads that are sense to unannotated or partial rRNA genes appear antisense to open reading "
                 "frames predicted on the opposite strand; sense and antisense counts were therefore also computed "
                 "without %s, and these counts were used for all antisense analyses." % rr)
    if X.on("RRNA_RESCUE") and X.on("RUN_CMSCAN"):
        items.append("In addition, rRNA genes missing from the annotation were recovered: contigs carrying protein-coding "
                     "genes with ≥ %s antisense reads and sense reads ≤ %s %% of the antisense reads were searched with "
                     "%s against the Rfam rRNA covariance models (5S, 5.8S, SSU, LSU; --cut_ga), hits not overlapping "
                     "an annotated rRNA were added as rRNA genes, and all counts were recomputed; protein-coding genes overlapping "
                     "a rescued rRNA gene on the opposite strand were flagged and excluded from functional and "
                     "multi-omic comparisons." % (
                         P("RRNA_RESCUE_MIN_ANTISENSE", "20"),
                         int(round(float(P("RRNA_RESCUE_MAX_SENSE_RATIO", "0.1")) * 100)), X.tool("infernal")))
    S.append(("Mapping to the sample's co-assembly and gene-level quantification", items))
    S.append(("Genomic bins", [
        "Contigs were assigned to metagenome-assembled genomes (%s) with taxonomy from %s; bins were classified as "
        "high, medium or low quality from completeness and contamination, and small RNA counts were summed per "
        "contig, bin and taxonomic rank." % (X.rc("binning"), X.rc("bin_taxonomy")),
    ]))
    items = [
        "Novel small RNA loci were detected from the assembly alignments%s: reads were merged per strand into "
        "clusters with %s (merge -s -d %s), and clusters of %s–%s nt with at least %s weighted reads were "
        "kept." % (
            " after removal of the rRNA reads defined above" if X.on("NOVEL_EXCLUDE_RRNA") else "",
            X.tool("bedtools"), P("NOVEL_MERGEDIST", "10"), P("NOVEL_MINLEN", "18"), P("NOVEL_MAXLEN", "250"),
            P("NOVEL_MINREADS", "10")),
        "The annotation was complemented with tRNA and tmRNA genes predicted by %s%s; clusters without a feature on "
        "their own strand were reported as intergenic (no feature on either strand) or antisense to a CDS." % (
            X.tool("aragorn"), (" and %s" % X.tool("trnascan")) if X.on("RUN_TRNASCAN", "0") else ""),
    ]
    if X.on("RUN_CMSCAN"):
        items.append("Candidate sequences were searched against the Rfam covariance models with %s (%s)." % (
            X.tool("infernal"), P("CMSCAN_OPTS", "--rfam --cut_ga")))
    if X.on("RUN_RNAFOLD"):
        items.append("Candidates without an Rfam hit were folded with %s (minimum free energy)." % X.tool("viennarna"))
    S.append(("Detection of novel small RNA loci", items))
    S.append(("tRNA and rRNA fragments and RNA integrity", [
        "Alignments lying ≥ %s %% within a tRNA, rRNA or tmRNA on the same strand were classified by their "
        "position on the parent RNA (tolerance %s nt; anticodon region %s–%s %% of the tRNA length) into 5′ and 3′ "
        "fragments, halves and internal fragments%s." % (
            int(round(float(P("FRAG_MIN_OVERLAP", "0.9")) * 100)), P("FRAG_TOL", "3"),
            int(round(float(P("FRAG_AC_FROM", "0.40")) * 100)), int(round(float(P("FRAG_AC_TO", "0.55")) * 100)),
            c("trf")),
        "Because size-selected small RNA libraries carry no 5′→3′ coverage information comparable to the RNA "
        "integrity number (RIN)%s, RNA integrity was assessed from how defined the fragment ends were: the share "
        "of tRNA reads with a mature 5′ or 3′ end, of internal tRNA fragments and of 3′ fragments carrying the "
        "terminal CCA, the normalised entropy and top-3 share of fragment end positions, the evenness of fragment "
        "midpoints along complete 16S and 23S rRNA genes, and the read length profile." % c("rin"),
        ("Integrity indicators were compared with the RIN of the long RNA of the same extractions." if P("RNA_QC")
         else ""),
    ]))
    S.append(("CRISPR arrays", [
        "Reads on annotated CRISPR arrays (≥ 50 % of the read inside) were counted per array; the transcribed "
        "strand was inferred from the read strands, reads on that strand were assigned to the spacer they "
        "overlapped by ≥ 10 nt (otherwise to the repeat), and CRISPR–Cas types were assigned from Cas genes within "
        "10 kb on the same contig.",
        "In type II and V loci, read clusters outside the array (≥ 10 weighted reads, 30–300 nt, within 1.5 kb of the "
        "locus) were tested as tracrRNA candidates by local alignment (Smith–Waterman, G·U pairs allowed; reads "
        "plus 30 nt upstream) with the "
        "reverse complement of the crRNA repeat; clusters with ≥ 15 aligned nt, ≥ 75 % base pairs and an empirical "
        "P ≤ 0.01 against 500 shuffled sequences were reported as candidates.",
    ]))
    S.append(("Comparison with metatranscriptomic and metagenomic data", [
        "For each sample, small RNA fractions, metatranscriptome and metagenome were compared pairwise per "
        "feature, contig and bin: entries with ≥ %s reads in at least one of the two sources were kept, counts "
        "were normalised to counts per million (CPM, also within feature type), log2 ratios were computed with a "
        "pseudocount of %s, and agreement was summarised with Spearman correlations%s; bin-level comparisons "
        "were repeated on CDS reads only and without rRNA-carrying contigs." % (
            P("COMPARE_MIN_COUNT", "5"), P("COMPARE_PSEUDOCOUNT", "1"), c("spearman")),
    ]))
    S.append(("Reports and software", [
        "Per-sample reports and all statistics were produced with %s, %s and %s; all tool versions are listed "
        "in %s." % (X.tool("python"), X.tool("numpy"), X.tool("matplotlib"),
                    "the Key resources table" if X.a.style == "star" else "Supplementary Table [CHECK]"),
    ]))
    return S


def cohort_facts(X):
    cd = X.a.cohort_dir or os.path.join(X.a.outdir or ".", "Cohort")
    st = read_tsv(os.path.join(cd, "tables", "cohort.sample_stats.tsv"))
    f = {"samples": len({r["sample"] for r in st}), "libs": collections.Counter(), "used": collections.Counter(),
         "groups": collections.defaultdict(collections.Counter), "dropped": collections.Counter()}
    for r in st:
        f["libs"][r["type"]] += 1
        if r.get("in_statistics", "yes") == "yes":
            f["used"][r["type"]] += 1
            f["groups"][r["type"]][r.get("group", "NA")] += 1
        else:
            f["dropped"][r["type"]] += 1
    return f


def cohort_sections(X):
    c, P = X.c, X.P
    f = cohort_facts(X)
    S = []
    q = P("COHORT_SIG_Q", "0.1")
    if f["samples"]:
        per = "; ".join("%s: %d libraries (%s)%s" % (
            t, f["used"][t], ", ".join("%s %d" % (g, n) for g, n in sorted(f["groups"][t].items())),
            (", %d excluded" % f["dropped"][t]) if f["dropped"][t] else "") for t in sorted(f["libs"]))
        comp = "The cohort analysis comprised %d samples (%s)." % (f["samples"], per)
    else:
        comp = "[CHECK: cohort composition — run with --cohort-dir]"
    drop = P("COHORT_QC_DROP", "FAIL")
    oz = P("COHORT_OUTLIER_Z", "4")
    S.append(("Cohort analysis: libraries and metadata", [
        comp,
        ("Libraries with a QC verdict of %s were excluded from all statistics." % (
            "FAIL" if drop == "FAIL" else "WARN or FAIL") if drop != "none" else
         "No library was excluded on the basis of the QC verdict."),
        ("Libraries deviating by more than %s robust z-scores (median and median absolute deviation within the "
         "fraction%s) in key metrics were flagged for inspection." % (oz, c("iglewicz")) if oz != "0" else ""),
        "Sample metadata were joined by sample identifier; the sRNA fractions were always analysed separately.",
    ]))
    S.append(("Feature tables", [
        "Per library, count tables were built for genera, families and species (Bracken estimates from rRNA "
        "reads), human miRNAs, Rfam families and types, RNA classes, tRNA amino acids, and sense reads per KEGG "
        "pathway, module and orthology, Pfam domain and COG.",
        "Features were analysed as relative abundances and retained if they had ≥ 10 reads in ≥ 10 % of the "
        "libraries of a fraction.",
        "miRNAs were analysed only in fractions with ≥ 6 libraries of ≥ %s human miRNA reads." % P(
            "COHORT_MIN_MIRNA_READS", "50"),
    ]))
    S.append(("Group comparisons", [
        "Library metrics and features were compared between groups with the Kruskal–Wallis test with tie "
        "correction%s and pairwise with the two-sided Mann–Whitney *U* test%s (exact for *n* ≤ 20; effect size "
        "P(X > Y)); *P* values were adjusted with the Benjamini–Hochberg procedure%s within each table, and "
        "*q* < %s was considered significant." % (c("kruskal"), c("mannwhitney"), c("bh"), q),
        "Paired comparisons of the intracellular and extracellular fraction of the same sample used the Wilcoxon "
        "signed-rank test%s; associations of integrity indicators with wet-lab RNA quality used Spearman "
        "correlations%s (BH within each QC measure); covariate balance between groups was tested with the "
        "Kruskal–Wallis test (numeric) or the χ² test (categorical)." % (c("wilcoxon"), c("spearman")),
    ]))
    items = [
        "Beta diversity was computed as Bray–Curtis dissimilarity%s of relative abundances, ordinated by principal "
        "coordinate analysis%s and tested with PERMANOVA%s (%s permutations, marginal per variable)." % (
            c("braycurtis"), c("pcoa"), c("permanova"), P("COHORT_PERMUTATIONS", "999")),
    ]
    if P("COHORT_NMDS", "genus"):
        items.append("Non-metric multidimensional scaling%s (two dimensions) of the Bray–Curtis dissimilarities was "
                     "fitted by SMACOF majorisation%s with monotone (isotonic) regression, %s random starts plus a "
                     "PCoA start, and Kruskal's stress-1 was reported." % (c("nmds"), c("smacof"),
                                                                        P("COHORT_NMDS_STARTS", "10")))
    S.append(("Community composition", items))
    S.append(("Antisense transcription", [
        "For each functional term and library with ≥ 10 reads, an antisense index was defined as antisense / "
        "(sense + antisense) reads, using the rRNA-filtered counts; indices were compared between groups as "
        "above.",
    ]))
    blk = P("COHORT_BLOCK") or P("COHORT_MATCH_COL")
    if blk:
        ref = P("COHORT_REFERENCE", "the first group")
        S.append(("Block design (%s)" % blk, [
            "Because samples sharing a %s are not independent, comparisons with the reference group (%s) were run "
            "both pooled (all libraries) and within %s levels containing both groups; the design was classified per "
            "group as nested (≤ 2 shared levels) or crossed (≥ 3)." % (blk, ref, blk),
            "Within-%s comparisons used the van Elteren stratified Wilcoxon test%s (weights 1/(*n* + 1), normal "
            "approximation) and PERMANOVA with permutations restricted to %s levels; nested designs were "
            "additionally contrasted between the reference subgroups to estimate the %s effect alone." % (
                blk, c("vanelteren"), blk, blk),
            "Poolability of the reference group across %s levels was tested by PERMANOVA within the reference "
            "group or, with too few reference libraries per level, on all libraries with %s after group (partial "
            "*R*²) and permutations within groups." % (blk, blk),
        ]))
    if P("COHORT_MODEL"):
        models = [m.strip() for m in P("COHORT_MODEL").split(",") if m.strip()]
        items = [
            "Covariate-adjusted associations were estimated with per-feature linear models, analogous to the "
            "default model of MaAsLin2%s, implemented in Python: %s." % (
                c("maaslin2"), "; ".join("M%d: feature ~ %s" % (i + 1, (m.split("|")[0].strip().replace("+", " + ") + (
                    " (%s only)" % m.split("|", 1)[1].split("=", 1)[-1].replace("/", ", ").strip() if "|" in m else "")))
                                         for i, m in enumerate(models))),
            "Outcomes were centred log-ratio-transformed abundances%s (pseudocount 0.5; features with ≥ 10 reads in "
            "≥ 10 %% of the libraries), library metrics and antisense indices; categorical terms were coded as "
            "treatment contrasts against %s and numeric terms scaled to standard-deviation units; complete cases "
            "only; features were modelled only with ≥ 20 libraries (≥ 40 for the antisense index; ≥ 3 per design level, "
            "≥ 5 distinct values, maximum leverage ≤ 0.95)." % (c("clr", "gloor"), P("COHORT_REFERENCE", "the reference group")),
            "Coefficients were tested with heteroskedasticity-consistent (HC3) standard errors%s and *t* tests, "
            "with Benjamini–Hochberg adjustment per term within each outcome set." % c("hc3"),
        ]
        if blk:
            items.append("For crossed block designs each model was additionally fitted with %s-clustered standard "
                         "errors (CR1, *t* with *G* − 1 degrees of freedom)%s and with %s fixed effects "
                         "(within-%s transformation)." % (blk, c("cluster"), blk, blk))
        S.append(("Adjusted models", items))
    S.append(("Statistical software", [
        "All cohort statistics were implemented in Python with %s (no R dependencies); the test implementations were "
        "checked against SciPy reference implementations. QC outputs of all libraries were aggregated with %s." % (
            X.tool("numpy"), X.tool("multiqc")),
    ]))
    return S


# ----------------------------------------------------------------------------------------------- rendering
def clean(items):
    return [re.sub(r"\s+", " ", x).strip() for x in items if x and x.strip()]


def references_md(X):
    out = ["## References", ""]
    for i, k in enumerate(X.c.order, 1):
        out.append("%d. %s" % (i, fmt_ref(k, X.a.style)))
    return out


def key_resources(X):
    rows = ["| REAGENT or RESOURCE | SOURCE | IDENTIFIER |", "|---|---|---|"]
    rows.append("| **Deposited data / reference databases** | | |")
    for k in X.used_res:
        r = X.res.get(k)
        if not r:
            continue
        name = r["name"] + (" (%s)" % r["version"] if r["version"] and not r["version"].startswith("set ") else "")
        src = r["source"] if r["source"] and not r["source"].startswith("set ") else "[CHECK]"
        src += X.c(*r["refs"]) if r["refs"] else ""
        ident = r["id"] if r["id"] and not r["id"].startswith("set ") else "[CHECK]"
        rows.append("| %s | %s | %s |" % (name, src, ident))
    rows.append("| **Software and algorithms** | | |")
    rows.append("| %s%s | this protocol | %s |" % (PIPELINE, " (%s)" % X.git if X.git else "", REPO_URL))
    for t in X.used_tools:
        name, _, refs, url = TOOLS[t]
        rows.append("| %s v%s | %s | %s |" % (name, X.v(t), ("ref." + X.c(*refs.split(","))) if refs else "—",
                                              url))
    return rows


def render_nature(X, sections, title):
    out = ["# %s" % title, ""]
    for head, items in sections:
        items = clean(items)
        if items:
            out += ["**%s.** %s" % (head, " ".join(items)), ""]
    out += ["**Data availability.** [CHECK: accession numbers of the raw sequencing data (e.g. ENA/SRA) and of the "
            "processed tables.]", "",
            "**Code availability.** The %s is available at %s%s; documentation: %s." % (
                PIPELINE, REPO_URL, " (version %s)" % X.git if X.git else "", DOCS_URL), ""]
    return out + references_md(X)


STAR_CMD = {
    "Read preprocessing and quality control": "qc,contam",
    "Separation of human and microbial reads": "prep,split",
    "Human small RNA classes": "human",
    "Microbial small RNA classes and rRNA-based taxonomy": "nonhuman",
    "Mapping to the sample's co-assembly and gene-level quantification": "assembly",
    "Genomic bins": "bins",
    "Detection of novel small RNA loci": "novel",
    "tRNA and rRNA fragments and RNA integrity": "fragments",
    "CRISPR arrays": "crispr",
    "Comparison with metatranscriptomic and metagenomic data": "--compare",
    "Reports and software": "--report",
}
PRESENT = [(r"\bwere\b", "are"), (r"\bwas\b", "is"), (r"\breceived\b", "receives"),
           (r"\bcontributed\b", "contributes"), (r"\bthat matched\b", "that match"),
           (r"\boverlapped\b", "overlaps"), (r"\bhad\b", "have"), (r"\bcomprised\b", "comprises"),
           (r"\bused\b", "use"), (r"\bwas considered\b", "is considered")]


def present(t):
    for a_, b_ in PRESENT:
        t = re.sub(a_, b_, t)
    return t


def render_star(X, sections, title, cohort):
    out = ["# %s" % title, "",
           "## Summary", "",
           "This protocol describes the processing of gut small RNA sequencing libraries (intracellular and "
           "extracellular fractions) with the %s together with the metagenome and metatranscriptome of the same "
           "samples%s. [CHECK: edit summary.]" % (PIPELINE, ", followed by a cross-sample cohort analysis" if cohort
                                                  else ""), "",
           "## Before you begin", "",
           "1. Process the metagenomic and metatranscriptomic reads of each sample with the integrated multi-omic "
           "workflow (co-assembly, annotation, binning).",
           "2. Install the software environments (envs/*.yml) and the reference databases listed in the key "
           "resources table, and set their paths in srna_pipeline.config.sh.",
           "3. Prepare a sample sheet with the small RNA FASTQ files and the multi-omic outputs per sample"
           "%s." % (" and a metadata table (one row per sample, group and covariate columns)" if cohort else ""), "",
           "## Key resources table", ""]
    # the table must be built after the text (citations and used tools), so render the body first
    body = ["## Step-by-step method details", ""]
    n = 0
    for head, items in sections:
        items = clean(items)
        if not items:
            continue
        if head.startswith(("Cohort analysis", "Feature tables", "Group comparisons", "Community", "Antisense",
                            "Block design", "Adjusted models", "Statistical software")):
            continue
        cmd = STAR_CMD.get(head)
        body += ["### %s" % head, ""]
        if cmd:
            n += 1
            run = ("srna_pipeline.sh -s samples.tsv -o <outdir> %s" % cmd if cmd.startswith("--") else
                   "srna_pipeline.sh -s samples.tsv -o <outdir> -m %s" % cmd)
            body.append("%d. Run `%s` (all modules in order: omit -m; on SLURM: submit_srna_pipeline.sh)." % (n, run)
                        if head.startswith("Read preprocessing") else "%d. Run `%s`." % (n, run))
            body += ["    %s. %s" % (chr(97 + i), present(it)) for i, it in enumerate(items)]
        else:
            for it in items:
                n += 1
                body.append("%d. %s" % (n, present(it)))
        body.append("")
    qa = [(h, clean(i)) for h, i in sections if h.startswith(("Cohort analysis", "Feature tables",
                                                             "Group comparisons", "Community", "Antisense",
                                                             "Block design", "Adjusted models",
                                                             "Statistical software"))]
    body += ["## Expected outcomes", "",
             "Per sample and fraction: read fate (trimming, human, rRNA, tRNA, Rfam, unknown), RNA class and "
             "taxonomic profiles, gene, functional-term and bin counts (sense and antisense), candidate novel loci, "
             "tRNA/rRNA fragment and integrity profiles, CRISPR array expression, comparisons with metaT and metaG, "
             "and a PDF/HTML report with a QC verdict%s." % (
                 "; for the cohort, group comparisons, ordinations and adjusted models with a cohort report"
                 if cohort else ""), "",
             "## Quantification and statistical analysis", ""]
    if cohort:
        body += ["Run `srna_pipeline.sh -s samples.tsv -o <outdir> --cohort` with COHORT_METADATA and the COHORT_* "
                 "settings of the configuration; the analysis below is applied to all finished samples.", ""]
    if qa:
        for h, items in qa:
            if items:
                body += ["**%s.** %s" % (h, present(" ".join(items))), ""]
    else:
        body += ["Counts are weighted 1/*n* for reads with *n* equally good alignments; relative measures are "
                 "compositional and compared within a fraction.", ""]
    body += ["## Limitations", "",
             "Short reads (< 20 nt) cannot be assigned reliably to the host or to single loci; multi-mapping reads are "
             "distributed fractionally; RNA-class shares are relative (compositional); antisense counts depend on "
             "the completeness of the rRNA annotation of the contigs; human miRNA counts in stool are low.", "",
             "## Troubleshooting", "",
             "Module status and failures are reported by srna_status.sh; failed rows can be resubmitted with "
             "--resubmit. See the pipeline documentation for module-specific problems.", "",
             "## Resource availability", "",
             "### Lead contact", "", "[CHECK: name and e-mail of the lead contact.]", "",
             "### Materials availability", "", "This study did not generate new materials.", "",
             "### Data and code availability", "",
             "[CHECK: accession numbers of the raw data.] The %s is available at %s%s; documentation: %s." % (
                 PIPELINE, REPO_URL, " (version %s)" % X.git if X.git else "", DOCS_URL), ""]
    body += references_md(X)
    return out + key_resources(X) + [""] + body


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["sample", "cohort", "all"], default="all")
    ap.add_argument("--style", choices=["nature", "star", "cell"], default="nature")
    ap.add_argument("--params", default="")
    ap.add_argument("--config", default="", help="srna_pipeline.config.sh to read the parameters from (instead of "
                    "--params; environment overrides apply)")
    ap.add_argument("--envs", default="", help="comma-separated environment prefixes (conda-meta)")
    ap.add_argument("--resources", default=os.path.join(HERE, "methods_resources.tsv"))
    ap.add_argument("--outdir", default="")
    ap.add_argument("--cohort-dir", default="")
    ap.add_argument("--samples", default="", help="(unused; kept for compatibility)")
    ap.add_argument("--types", default="", help="sRNA types (default: from the cohort table or the params)")
    ap.add_argument("--out", required=True, help="output Markdown")
    ap.add_argument("--docx", action="store_true", help="also write .docx with pandoc (if available)")
    a = ap.parse_args()
    if a.style == "cell":
        a.style = "star"
    if a.config and not a.params:
        dump = ('source "$1" >/dev/null 2>&1; for v in $(grep -oE "^[A-Z_][A-Z0-9_]*=" "$1" | tr -d = | sort -u); '
                'do printf "%s\\t%s\\n" "$v" "${!v-}"; done')
        pf = a.out + ".params.tsv"
        with open(pf, "w") as fh:
            subprocess.run(["bash", "-c", dump, "_", a.config], stdout=fh, check=False)
        a.params = pf
    X = Ctx(a)
    types = [t for t in a.types.split(",") if t]
    if not types and a.mode != "sample":
        types = sorted(cohort_facts(X)["libs"])
    sections = processing_sections(X, types)
    cohort = a.mode in ("cohort", "all")
    if cohort:
        sections += cohort_sections(X)
    title = "Methods" if a.style == "nature" else "STAR Protocols: %s" % PIPELINE
    lines = render_nature(X, sections, title) if a.style == "nature" else render_star(X, sections, title, cohort)
    stamp = "<!-- generated by srna_methods.py (%s, style %s, mode %s) on %s; check every [CHECK: …] -->" % (
        PIPELINE, a.style, a.mode, datetime.date.today().isoformat())
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(stamp + "\n\n" + "\n".join(lines).rstrip() + "\n")
    n_check = sum(l.count("[CHECK") for l in lines)
    sys.stderr.write("[methods] %s (%s, %d references, %d [CHECK] placeholders)\n" % (
        a.out, a.style, len(X.c.order), n_check))
    if a.docx:
        pd = shutil.which("pandoc")
        if pd:
            subprocess.run([pd, a.out, "-f", "markdown+superscript", "-o", os.path.splitext(a.out)[0] + ".docx"],
                           check=False)
        else:
            sys.stderr.write("[methods] pandoc not found — no .docx\n")


if __name__ == "__main__":
    main()
