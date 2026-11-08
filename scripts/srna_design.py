#!/usr/bin/env python3
"""srna_design.py — what a run contains, so that reports and Methods describe only what is there.

The study design is read from the inputs each finished library recorded
(<outdir>/<sample>/log_cluster_jobs/<type>/inputs.tsv):

  types      sRNA library types found (free labels, e.g. sRNA, sRNA-Int, sRNA-Ext)
  int_ext    True only if both sRNA-Int and sRNA-Ext are present (cell fraction vs cell-free fraction);
             fraction-specific text and the paired Int vs Ext analyses depend on it
  metaT      True if any library has metaT feature or contig counts
  metaG      True if any library has metaG feature or contig counts (or an IMP annotation directory)

Text helpers (used by srna_report.py, srna_cohort.py and srna_methods.py):

  D.unit       "fraction" when Int and Ext are present, else "library type"
  D.omics      "metagenome and metatranscriptome" / "metagenome" / "" (what the sRNA is compared with)
  D.assembly   "co-assembly of the metagenomic and metatranscriptomic reads" or "assembly of the metagenomic reads"
  D.libraries  "small RNA libraries (intracellular and extracellular fractions)" or "small RNA libraries"
  D.sample_type  free label of the sample material (SAMPLE_TYPE, default "stool")
  D.host       host word (HOST_NAME, default "human")
  D.host_only  library types run without an assembly (assembly and gff '-', e.g. blood sRNA)
"""
import glob
import os

FRACTIONS = ("sRNA-Ext", "sRNA-Int")


def _read_inputs(path):
    d = {}
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t", 1)
                if f and f[0]:
                    d[f[0]] = f[1].strip() if len(f) > 1 else ""
    except OSError:
        pass
    return d


def _given(v):
    return bool(v) and v not in ("-", "NA", "none")


class Design:
    def __init__(self, types=(), metaT=False, metaG=False, sample_type="stool", host="human", host_only=()):
        self.types = list(types)
        self.host_only = [t for t in host_only if t in self.types] if self.types else list(host_only)
        self.metaT = bool(metaT)
        self.metaG = bool(metaG)
        self.sample_type = sample_type or "stool"
        self.host = (host or "human").strip()

    @property
    def int_ext(self):
        return all(t in self.types for t in FRACTIONS)

    @property
    def unit(self):
        return "fraction" if self.int_ext else "library type"

    @property
    def units(self):
        return "fractions" if self.int_ext else "library types"

    @property
    def omics(self):
        if self.metaT and self.metaG:
            return "metagenome and metatranscriptome"
        if self.metaT:
            return "metatranscriptome"
        if self.metaG:
            return "metagenome"
        return ""

    @property
    def omics_short(self):
        return " / ".join(s for s, on in (("metaG", self.metaG), ("metaT", self.metaT)) if on)

    @property
    def assembly(self):
        if self.metaT:
            return "co-assembly of the metagenomic and metatranscriptomic reads"
        return "assembly of the metagenomic reads"

    @property
    def libraries(self):
        if self.int_ext:
            return "small RNA libraries (intracellular and extracellular fractions)"
        return "small RNA libraries"

    def __repr__(self):
        return "Design(types=%s, int_ext=%s, metaT=%s, metaG=%s, sample_type=%s, host=%s)" % (
            self.types, self.int_ext, self.metaT, self.metaG, self.sample_type, self.host)


def design(outdir="", types=None, samples=None, sample_type="stool", host="human"):
    """Design of a run: from the inputs.tsv of all (or the given) samples under outdir.
    types: library types to use instead of the ones found; samples: sample IDs to restrict to."""
    found, metaT, metaG, host_only, with_asm = [], False, False, set(), set()
    keep = set(samples) if samples else None
    if outdir:
        for p in sorted(glob.glob(os.path.join(outdir, "*", "log_cluster_jobs", "*", "inputs.tsv"))):
            s = p.split(os.sep)[-4]
            if keep is not None and s not in keep:
                continue
            d = _read_inputs(p)
            t = d.get("type") or p.split(os.sep)[-2]
            if t not in found:
                found.append(t)
            (host_only if d.get("library") == "host_only" or d.get("assembly") in ("", "-") else with_asm).add(t)
            metaT = metaT or _given(d.get("metaT_features")) or _given(d.get("metaT_contigs"))
            metaG = metaG or _given(d.get("metaG_features")) or _given(d.get("metaG_contigs"))
    t = [x for x in (types or found) if x]
    t.sort(key=lambda x: (x not in FRACTIONS, FRACTIONS.index(x) if x in FRACTIONS else 0, x))
    return Design(t, metaT, metaG, sample_type, host, sorted(host_only - with_asm))
