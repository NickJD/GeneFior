#!/bin/bash

# Genefíor Database Builder
# Usage: ./build_database.sh database_name dna.fasta[.gz] aa.fasta[.gz] [threads]
#
# Both DNA and AA FASTA files may be plain or gzip-compressed (.gz).
# The base name used for index files is derived by stripping ALL extensions
# from the input filename so that inputs like "seq.fas.fasta.gz" always
# produce clean index names such as "seq_blastdb", "seq_bowtie2db", etc.

set -e

DB_NAME=$1
DNA_FASTA=$2
AA_FASTA=$3
THREADS=${4:-4}

if [ -z "$DB_NAME" ] || [ -z "$DNA_FASTA" ] || [ -z "$AA_FASTA" ]; then
    echo "Usage: $0 <database_name> <dna.fasta[.gz]> <aa.fasta[.gz]> [threads]"
    exit 1
fi

# ---------------------------------------------------------------------------
# Helper: strip ALL extensions from a filename to get a clean base name.
# e.g.  VFDB_setA_nt.fas.fasta.gz  ->  VFDB_setA_nt
#       genes.faa                  ->  genes
# ---------------------------------------------------------------------------
strip_all_extensions() {
    local name
    name=$(basename "$1")
    # Remove .gz first if present
    name="${name%.gz}"
    # Then keep removing the last extension until no dot remains
    while [[ "$name" == *.* ]]; do
        name="${name%.*}"
    done
    echo "$name"
}

DNA_BASE=$(strip_all_extensions "$DNA_FASTA")
AA_BASE=$(strip_all_extensions "$AA_FASTA")

echo "Building Genefíor database: $DB_NAME"
echo "DNA sequences:  $DNA_FASTA  (base: $DNA_BASE)"
echo "AA sequences:   $AA_FASTA  (base: $AA_BASE)"
echo "Threads: $THREADS"
echo ""

# ---------------------------------------------------------------------------
# Decompress inputs if gzipped – most indexers don't accept .gz natively.
# We write to a temp directory and clean up at the end.
# ---------------------------------------------------------------------------
TMPDIR_LOCAL=$(mktemp -d)
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

decompress() {
    local src="$1"
    local dest="$2"
    if [[ "$src" == *.gz ]]; then
        echo "Decompressing $src ..."
        gzip -dc "$src" > "$dest"
    else
        # Use a symlink to avoid copying large files
        ln -sf "$(realpath "$src")" "$dest"
    fi
}

DNA_PLAIN="$TMPDIR_LOCAL/${DNA_BASE}.fasta"
AA_PLAIN="$TMPDIR_LOCAL/${AA_BASE}.fasta"
decompress "$DNA_FASTA" "$DNA_PLAIN"
decompress "$AA_FASTA"  "$AA_PLAIN"

# ---------------------------------------------------------------------------
# Create directory structure
# ---------------------------------------------------------------------------
mkdir -p "${DB_NAME}"/{blast_aa,blast_dna,bowtie2,bwa,diamond,minimap2}

# Store compressed copies of source files for reference
echo "Storing compressed source files..."
gzip -c "$DNA_PLAIN" > "${DB_NAME}/${DNA_BASE}.fasta.gz"
gzip -c "$AA_PLAIN"  > "${DB_NAME}/${AA_BASE}.fasta.gz"

# ---------------------------------------------------------------------------
# Build BLAST DNA database
# ---------------------------------------------------------------------------
echo "Building BLAST DNA database..."
makeblastdb -in "$DNA_PLAIN" -dbtype nucl \
    -out "${DB_NAME}/blast_dna/${DNA_BASE}_blastdb" \
    -parse_seqids

# ---------------------------------------------------------------------------
# Build BLAST AA database
# ---------------------------------------------------------------------------
echo "Building BLAST AA database..."
makeblastdb -in "$AA_PLAIN" -dbtype prot \
    -out "${DB_NAME}/blast_aa/${AA_BASE}_blastdb" \
    -parse_seqids

# ---------------------------------------------------------------------------
# Build DIAMOND database
# ---------------------------------------------------------------------------
echo "Building DIAMOND database..."
diamond makedb --in "$AA_PLAIN" \
    --db "${DB_NAME}/diamond/${AA_BASE}_diamonddb" \
    --threads "$THREADS"

# ---------------------------------------------------------------------------
# Build Bowtie2 index
# ---------------------------------------------------------------------------
echo "Building Bowtie2 index..."
bowtie2-build --threads "$THREADS" "$DNA_PLAIN" \
    "${DB_NAME}/bowtie2/${DNA_BASE}_bowtie2db"

# ---------------------------------------------------------------------------
# Build BWA index
# ---------------------------------------------------------------------------
echo "Building BWA index..."
bwa index -p "${DB_NAME}/bwa/${DNA_BASE}_bwadb" "$DNA_PLAIN"

# ---------------------------------------------------------------------------
# Build Minimap2 index
# ---------------------------------------------------------------------------
echo "Building Minimap2 index..."
minimap2 -d "${DB_NAME}/minimap2/${DNA_BASE}_minimap2db" "$DNA_PLAIN"

# ---------------------------------------------------------------------------
# Validation: check expected files exist
# ---------------------------------------------------------------------------
echo ""
echo "Validating output files..."
ERRORS=0
check_exists() {
    if ls "$1" 1>/dev/null 2>&1; then
        echo "  ✓  $1"
    else
        echo "  ✗  MISSING: $1"
        ERRORS=$((ERRORS+1))
    fi
}

check_exists "${DB_NAME}/blast_dna/${DNA_BASE}_blastdb.nhr"
check_exists "${DB_NAME}/blast_aa/${AA_BASE}_blastdb.phr"
check_exists "${DB_NAME}/diamond/${AA_BASE}_diamonddb.dmnd"
check_exists "${DB_NAME}/bowtie2/${DNA_BASE}_bowtie2db.1.bt2"
check_exists "${DB_NAME}/bwa/${DNA_BASE}_bwadb.amb"
check_exists "${DB_NAME}/minimap2/${DNA_BASE}_minimap2db"

echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo "Database build complete — all index files present."
else
    echo "WARNING: $ERRORS expected file(s) missing — check the build output above."
    exit 1
fi

echo ""
echo "Database location: $DB_NAME"
echo ""
echo "To use with GeneFíor:"
echo "  GeneFior -i reads.fasta -st Single-FASTA -o results/ \\"
echo "    --db-path $(pwd)/$DB_NAME"
