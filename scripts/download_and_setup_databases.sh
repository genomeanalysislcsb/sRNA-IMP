#!/usr/bin/env bash
set -euo pipefail

# Database download and setup script for sRNA-IMP
# Edit config_examples/database_paths.example.sh and pass it here if needed:
#   bash scripts/download_and_setup_databases.sh config_examples/database_paths.example.sh

CONFIG="${1:-config_examples/database_paths.example.sh}"
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG" >&2
  exit 1
fi
source "$CONFIG"

for exe in wget curl gunzip tar bowtie-build kraken2-build bracken-build cmpress; do
  if ! command -v "$exe" >/dev/null 2>&1; then
    echo "ERROR: required executable missing: $exe" >&2
    exit 1
  fi
done

mkdir -p "$DB_ROOT" "$HUMAN_DIR" "$NONHUMAN_DIR" "$SHARED_DIR"

echo "[1/8] Rfam"
mkdir -p "$RFAM_DIR"
if [[ ! -f "$RFAM_DIR/Rfam.cm" ]]; then
  wget -O "$RFAM_DIR/Rfam.cm.gz" "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
  gunzip -f "$RFAM_DIR/Rfam.cm.gz"
fi
if [[ ! -f "$RFAM_DIR/Rfam.clanin" ]]; then
  wget -O "$RFAM_DIR/Rfam.clanin" "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
fi
if [[ ! -f "$RFAM_DIR/Rfam.cm.i1f" ]]; then
  (cd "$RFAM_DIR" && cmpress Rfam.cm)
fi

echo "[2/8] Kraken2 standard database"
mkdir -p "$KRAKEN_DIR"
if [[ ! -f "$KRAKEN_DIR/hash.k2d" ]]; then
  kraken2-build --standard --db "$KRAKEN_DIR" --threads "$THREADS"
fi

echo "[3/8] Bracken build for read length $READ_LENGTH"
mkdir -p "$BRACKEN_DIR"
# Use the Kraken DB directory directly; Bracken writes into the same DB tree by default.
if [[ ! -f "$KRAKEN_DIR/database${READ_LENGTH}mers.kmer_distrib" ]]; then
  bracken-build -d "$KRAKEN_DIR" -t "$THREADS" -k 35 -l "$READ_LENGTH"
fi
ln -snf "$KRAKEN_DIR" "$BRACKEN_DIR/db"

echo "[4/8] Human miRNA (miRBase mature.fa)"
mkdir -p "$HUMAN_MIRNA_DIR"
if [[ ! -f "$HUMAN_MIRNA_DIR/mature.fa" ]]; then
  wget -O "$HUMAN_MIRNA_DIR/mature.fa" "https://www.mirbase.org/download/CURRENT/mature.fa"
fi
grep -A1 '^>hsa' "$HUMAN_MIRNA_DIR/mature.fa" | sed '/^--$/d' > "$HUMAN_MIRNA_DIR/hsa_mature.fa"
if [[ ! -f "$HUMAN_MIRNA_DIR/hsa_mature.1.ebwt" ]]; then
  bowtie-build "$HUMAN_MIRNA_DIR/hsa_mature.fa" "$HUMAN_MIRNA_DIR/hsa_mature"
fi

echo "[5/8] Optional non-human miRNA reference from miRBase"
mkdir -p "$NONHUMAN_MIRNA_DIR"
grep -A1 -v '^>hsa' "$HUMAN_MIRNA_DIR/mature.fa" | sed '/^--$/d' > "$NONHUMAN_MIRNA_DIR/nonhuman_mature.fa"
if [[ ! -f "$NONHUMAN_MIRNA_DIR/nonhuman_mature.1.ebwt" ]]; then
  bowtie-build "$NONHUMAN_MIRNA_DIR/nonhuman_mature.fa" "$NONHUMAN_MIRNA_DIR/nonhuman_mature"
fi

echo "[6/8] Placeholder human rRNA / tRNA / other ncRNA directories"
mkdir -p "$HUMAN_TRNA_DIR" "$HUMAN_RRNA_DIR" "$HUMAN_OTHERRNA_DIR"
mkdir -p "$NONHUMAN_TRNA_DIR"

cat > "$HUMAN_TRNA_DIR/README.txt" <<'EOF'
Populate this directory with curated human tRNA FASTA files, e.g. nuclear + mitochondrial tRNAs, then build:
  cat human_tRNA.fa human_mt_tRNA.fa > hsa_tRNA_all.fa
  bowtie-build hsa_tRNA_all.fa hsa_tRNA_all
EOF

cat > "$HUMAN_RRNA_DIR/README.txt" <<'EOF'
Populate this directory with curated human cytosolic + mitochondrial rRNA FASTA files, then build:
  bowtie-build hsa_rRNA.fa hsa_rRNA
EOF

cat > "$HUMAN_OTHERRNA_DIR/README.txt" <<'EOF'
Populate this directory with curated human ncRNA FASTA (e.g. Ensembl ncRNA), then build:
  bowtie-build Homo_sapiens.GRCh38.ncrna.fa hsa_other_ncRNA
EOF

cat > "$NONHUMAN_TRNA_DIR/README.txt" <<'EOF'
Populate this directory with combined non-human tRNA FASTA (bacterial, fungal/eukaryotic, plant chloroplast, plant mitochondrial), then build:
  bowtie-build trna_all.fa trna_all
EOF

echo "[7/8] Environment and tool notes"
cat > "$DB_ROOT/tool_notes.txt" <<EOF
RiboDetector: install via Bioconda in your environment; no large external DB setup required.
ARAGORN: install via Bioconda.
tRNAscan-SE: install via Bioconda or source; depends on Infernal.
MirGeneDB: if preferred over miRBase for animals, download mature sequences manually and build Bowtie indices.
EOF

echo "[8/8] Database manifest"
cat > "$DB_ROOT/database_manifest.txt" <<EOF
Rfam_cm=$RFAM_DIR/Rfam.cm
Rfam_clanin=$RFAM_DIR/Rfam.clanin
Kraken2_db=$KRAKEN_DIR
Bracken_db=$KRAKEN_DIR
Bracken_read_length=$READ_LENGTH
Human_miRNA_index=$HUMAN_MIRNA_DIR/hsa_mature
Nonhuman_miRNA_index=$NONHUMAN_MIRNA_DIR/nonhuman_mature
EOF

echo "Done. Review README placeholders for the references that require curated local FASTA files."
