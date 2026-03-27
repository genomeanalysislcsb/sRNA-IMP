# ncRNA Discovery

The novel ncRNA branch uses non-host reads together with the optional assembly FASTA and GFF.

Main stages:
- build Bowtie index on the assembly
- map reads to contigs
- extract strand-aware BED intervals
- cluster supported loci
- remove overlaps with existing annotations
- retain intergenic and antisense-to-CDS candidates
- optionally filter candidates against Rfam and fold them with RNAfold

Key summary output:
- `<sample>.novel.summary.tsv`
