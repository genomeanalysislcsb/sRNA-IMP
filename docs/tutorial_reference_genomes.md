# Tutorial 2: without a metagenome — MGnify reference genomes

sRNA-IMP is built around an assembly of the sample's own metagenome or metatranscriptome. Many small RNA studies
have neither. Then a set of annotated reference genomes can take the place of the assembly
([Inputs](inputs.md#reference-genomes-instead-of-an-assembly)): reads are still placed on genes, strands and
genomes, and novel small RNA loci are still called, on the reference genomes.

This tutorial picks the reference genomes from the libraries themselves: the rRNA reads of each library are
classified anyway (module `nonhuman`, Kraken 2 and Bracken), and the species found there are looked up in an
[MGnify genome catalogue](https://www.ebi.ac.uk/metagenomics/browse/genomes), which provides annotated
species-representative genomes for host-associated biomes (human gut, oral cavity and vagina, mouse gut, pig and
chicken gut, cow rumen and others).

It uses the stool small RNA libraries of [Tutorial 1](tutorial_tarallo.md) and ignores their metagenomes, so the
result can be compared with the run on each sample's own assembly.

| step | what | tool |
| --- | --- | --- |
| 1 | run the libraries without an assembly (host split, RNA classes, rRNA taxonomy) | `submit_srna_pipeline.sh`, assembly `-` |
| 2 | choose and download MGnify genomes from the libraries' rRNA species | `scripts/srna_fetch_mgnify.py` |
| 3 | turn the genomes into sRNA-IMP inputs | `scripts/srna_prepare_inputs.py --genomes` |
| 4 | sample sheet and configuration for the reference set | `ASSEMBLY_SOURCE=reference` |
| 5 | run the assembly-based modules on the reference genomes, then reports and cohort | `submit_srna_pipeline.sh` |
| 6 | compare with the run on the sample's own assembly | `paper/compare_reference_mode.py` |

You need what Tutorial 1 needs up to its step 3 (sRNA-IMP, databases, the small RNA reads; no IMP3, no
metagenomes), and internet access on the node where you run step 2.

## 0. Folders

```bash
P=/path/to/sRNA-IMP                  # sRNA-IMP checkout
B=/path/to/tarallo2019               # small RNA reads in fastq/srna/ (Tutorial 1, steps 1-3)
R=$B/reference_mode                  # this tutorial
mkdir -p $R && cd $R
```

## 1. Run the libraries without an assembly

A library whose `assembly` and `gff` columns are `-` runs the modules that need no assembly (`qc contam split
human nonhuman`): trimming and QC, the host split (without the contig rescue), human and microbial RNA classes, and
the rRNA taxonomy. Take the first three columns of the Tutorial 1 sample sheet and set the others to `-`:

```bash
awk -F'\t' -v OFS='\t' '/^#/ || NF < 3 {next} {print $1,$2,$3,"-","-","-","-","-","-","-","-"}' \
    $B/samples.tsv > samples.noasm.tsv                       # the Tutorial 1 sample sheet
cp $B/tarallo2019.config.sh reference.config.sh             # same databases and settings
$P/submit_srna_pipeline.sh -s samples.noasm.tsv -o $R/out -c reference.config.sh \
    -m qc,contam,split,human,nonhuman -n 16 -T 6:00:00 --mem 64G -J ref1 --no-compare
```

Each library now has its rRNA species in `out/<ID>/Analysis/sRNA/nonhuman/rrna_tax/<ID>.rRNA.bracken.S.tsv`.

## 2. MGnify genomes from the rRNA species

```bash
ls $R/out/*/Analysis/sRNA/nonhuman/rrna_tax/*.rRNA.bracken.S.tsv | paste -sd, > bracken.list
python3 $P/scripts/srna_fetch_mgnify.py --catalogue human-gut --version v2.0.2 --outdir genomes \
    --bracken "$(cat bracken.list)" --min-reads 10
```

Check the latest catalogue version first
(`https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/`); for other body sites or hosts use
another catalogue (`human-oral`, `human-vaginal`, `mouse-gut`, …). The script takes every species with at least
`--min-reads` Bracken reads in any library, one species-representative genome each, and downloads its GFF (Prokka
annotation with the sequence). It writes:

| file | content |
| --- | --- |
| `genomes/<MGYG…>.gff.gz` | one genome per species |
| `genome_taxonomy.tsv`, `genome_quality.tsv` | GTDB lineage, CheckM completeness and contamination from the catalogue |
| `selection.tsv` | why each genome was chosen (species and its rRNA reads) |
| `unmatched_species.txt` | species without a representative: eukaryotes from food, taxa outside the biome, or names that differ between GTDB releases (the Kraken 2 database and the catalogue) |

One reference set for the whole cohort keeps the libraries comparable. It stops above 500 genomes
(`--max-genomes`); raise `--min-reads` to keep the set small, since the bowtie index of the set is built for every
sample. `--species`, `--genomes` and `--taxa` add genomes by name, accession or taxon.

## 3. Inputs from the genomes

```bash
python3 $P/scripts/srna_prepare_inputs.py --sample REF --outdir prep --genomes genomes \
    --genome-taxonomy genomes/genome_taxonomy.tsv --checkm genomes/genome_quality.tsv
```

The genomes are concatenated into `prep/REF/assembly.fa.gz` and `annotation.gff.gz` (contigs and IDs renamed
`<genome>__<name>`), each genome becomes one bin, and the MGnify `ncRNA` features from the Rfam search are given
their own types (tRNA, rRNA, tmRNA, riboswitches and leaders as `regulatory_region`); the script prints how many.
The MGnify GFFs have no CRISPR arrays and no KEGG, Pfam or COG terms in the form sRNA-IMP reads, so module `crispr`
finds no arrays and the functional summaries stay empty.

## 4. Sample sheet and configuration

Every library gets the same reference inputs:

```bash
awk -F'\t' -v OFS='\t' -v D=$R/prep/REF '{$4=D"/assembly.fa.gz"; $5=D"/annotation.gff.gz"; $6=D"/contig2bin.tsv.gz";
    $7=D"/bins.tsv"; print}' samples.noasm.tsv > samples.ref.tsv
printf 'ASSEMBLY_SOURCE=reference\n' >> reference.config.sh
```

`ASSEMBLY_SOURCE=reference` changes only the wording of reports and generated Methods (a reference set, not an
assembly of the sample).

## 5. Run on the reference genomes

The host split is repeated because the contig rescue now has contigs to check; everything after it follows:

```bash
$P/submit_srna_pipeline.sh -s samples.ref.tsv -o $R/out -c reference.config.sh \
    -m prep,split,human,nonhuman,mirna_origin,assembly,novel,fragments,bins,crispr -f \
    -n 16 -T 12:00:00 --mem 64G -J ref2
```

As in Tutorial 1, a final job writes the reports and the cohort summary (`out/Cohort/`). Fig. 3 can be drawn from
this run as well (`python3 $P/paper/fig3_tarallo.py $R/out $R/fig3_reference`).

## 6. Compare with the sample's own assembly

```bash
python3 $P/paper/compare_reference_mode.py $B/sRNA-IMP $R/out VF02 sRNA $R/VF02_reference
```

`VF02_reference.metrics.tsv` lists, for the run on the sample's own assembly and the run on the reference genomes:
host reads and the host hits returned by the contig rescue, the share of non-host reads mapped, the shares of
mapped reads on tRNA, rRNA, CDS, antisense CDS, ncRNA, tmRNA and intergenic regions, novel loci and CRISPR arrays.
`VF02_reference.genera.tsv` has the reads per genus in both runs; the script prints their Spearman correlation over
the genera seen in both and the share of reads those genera carry.

What to expect: the reference set recovers the abundant species, so tRNA and rRNA reads and the genus profile
agree well, while fewer non-host reads are placed (strains and species without a catalogue genome stay unmapped),
fewer host hits are rescued, and novel loci are found only on the reference genomes. The assembly of the sample
remains the better input whenever a metagenome or metatranscriptome exists.
