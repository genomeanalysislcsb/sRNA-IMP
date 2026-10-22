# Input & sample sheet

## Sample sheet

Tab- or space-separated text file; lines starting with `#` and a header line (`sample_id …`) are ignored.
**One row = one sample + one sRNA type.**

```text
sample_id  type      raw_fastq                    assembly                         gff                            contig2bin  gtdbtk  metaT_features  metaT_contigs  metaG_features  metaG_contigs
BXXX1      sRNA-Ext  /data/raw/Ext/BXXX1.fastq.gz /data/imp/bxxx1/assembly.fa.gz   /data/imp/bxxx1/annotation.gff.gz  /data/bxxx1.c2b  -  /data/imp/bxxx1/Analysis/annotation  /data/imp/bxxx1/Stats  -  -
BXXX1      sRNA-Int  /data/raw/Int/BXXX1.fastq.gz /data/imp/bxxx1/assembly.fa.gz   /data/imp/bxxx1/annotation.gff.gz
```

| column | content | required |
| --- | --- | --- |
| 1 `sample_id` | sample name (letters, digits, `.`, `_`, `-`) | yes |
| 2 `type` | sRNA type, e.g. `sRNA-Ext`, `sRNA-Int` | yes |
| 3 `raw_fastq` | raw single-end small-RNA reads, FASTQ or FASTQ.gz | yes, full path |
| 4 `assembly` | contigs of the sample (IMP metaG/metaT assembly), FASTA(.gz) | yes, full path |
| 5 `gff` | annotation of the contigs (GFF3, .gz allowed) | yes, full path |
| 6 `contig2bin` | `contig<TAB>bin` table (.gz allowed) | optional |
| 7 `gtdbtk` | bin table with the GTDB `classification` column: GTDB-Tk summary TSV(s), comma-separated for bac120 + ar53, or a combined bin-stats table (first column = bin name); `Completeness` / `Contamination` columns (CheckM/CheckM2) are used for the bin quality; .gz allowed | optional |
| 8 `metaT_features` | metaT read counts per GFF feature (featureCounts, `-g ID`): one table, a comma-separated list, the IMP `Analysis/annotation` directory, or a tar archive of IMP tables (`mt.counts.gz`, `per_contig_counts.tar.gz`) | optional |
| 9 `metaT_contigs` | metaT reads per contig: IMP `Stats/mt/mgmt.assembly.contig_read_count.txt.gz`, the IMP `Stats` directory, a tar of these, or a featureCounts table per contig; `-` = summed from the feature tables | optional |
| 10 `metaG_features` | metaG read counts per GFF feature (IMP `mg.ID.<type>_counts…`), same forms as column 8; `-` = the `mg.ID` tables found in column 8 | optional |
| 11 `metaG_contigs` | metaG reads per contig, same forms as column 9; `-` = the `mg/` file found in column 9, else summed from the metaG feature tables | optional |

* Columns 3–5 must be **absolute paths**; a row with a missing or relative path is reported and skipped,
  the other rows still run.
* Columns 6–11 belong to the *sample*, not the row: give them on any one row of the sample, `-` = not
  available. They are only used by `fragments`, `bins` and `--compare` (`srna_pipeline.sh` only).
* The same sample may appear with several types. Assembly and annotation are staged once
  (`<sample>/Assembly`, `<sample>/Analysis/annotation`) and shared by all its types.

## Raw reads

* Single-end small-RNA reads, one file per library, no UMIs.
* The 3′ adapter must be set in the config (`ADAPTER_3P`). The expobiome libraries use the
  Illumina TruSeq/NEBNext read-1 adapter `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` for both sRNA-Ext and
  sRNA-Int — **not** the TruSeq Small RNA adapter `TGGAATTCTCGG…` (see [QC](qc.md#wrong-adapter)).
* Reads are copied into `Preprocessing/<type>/` (`STAGE_MODE=copy`) or symlinked (`STAGE_MODE=symlink`).

## Assembly and annotation

* Contigs: the GFF and FASTA must use the same contig names.
* The GFF is normalised on staging (`?` strand → `.`) and converted to BED with `gff2bed`.
* Feature types counted by the `assembly` module: `ASSEMBLY_FEATURES`
  (default `CDS CRISPR crispr-repeat crispr-spacer ncRNA oriC regulatory_region rRNA tmRNA tRNA`).
* If the annotation contains no tRNAs, `fragments` uses the ARAGORN tRNA/tmRNA predictions of `novel`.

## metaT and metaG counts from IMP

IMP writes one featureCounts table per feature type and omic (`mt.ID.<type>_counts.tsv.gz` for metaT,
`mg.ID.<type>_counts.tsv.gz` for metaG, keyed by the GFF
`ID`, reverse-stranded `-s 2`, multi-mappers split with `-M --fraction`) into `Analysis/annotation/`.
Point column 8 at that directory — nothing has to be packed or converted:

```text
/path/to/imp/bxxx1/Analysis/annotation          # all mt.ID.* and mg.ID.* tables are read
```

A tar archive of the same files (e.g. `tar cfz mt.counts.gz …/Analysis/annotation/mt*counts.tsv.gz`) or a
comma-separated list of tables works as well. Intermediate files (`*.paired.tsv`, `*.single.tsv`) are
ignored, so no read is counted twice.

**metaG is picked up automatically:** the `mg.ID.*` tables in the column-8 directory or archive become
the source `metaG`; column 10 is only needed when the metaG tables are somewhere else. IMP's
`per_contig_counts.tar.gz` (`mt|mg.ID.<type>_per_contig.counts.tsv`) can be used instead — despite the
name these are the same per-feature tables (identical numbers), not counts per contig. If both forms are
present, each feature type is read once.

* All count columns are summed (IMP has one for paired and one for single reads).
* The feature type is taken from the file name (`mt.ID.CDS_counts…` → `CDS`).
* `*.Name.*` (keyed by product name) and `*.mantis.*` (functional categories) are not used.
### Contig counts

IMP counts all reads per contig in `Stats/<mt|mg>/mgmt.assembly.contig_read_count.txt.gz`
(`contig  start  end  reads`). These include intergenic reads, as the sRNA contig counts do, so use them
for the contig and bin levels:

```text
col 9   /path/to/imp/bxxx1/Stats                 # mt/ and mg/ file both found
col 9   /path/to/BXXX1/mt.contig.tar.gz          # or one tarball per omic ...
col 11  /path/to/BXXX1/mg.contig.tar.gz          # ... with the mg file in column 11
```

A file under an `mt/` or `mg/` directory (or tar member) is used only for that omic; column 11 may stay
`-` when column 9 already contains the `mg/` file. Without contig tables the contig counts are the sum of
the CDS, tRNA, rRNA, ncRNA, tmRNA and CRISPR feature counts per contig, which misses intergenic reads
(BXXX1 metaT: 45.2 M feature reads vs 56.2 M contig reads).
