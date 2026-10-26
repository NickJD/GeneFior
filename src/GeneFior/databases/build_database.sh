#!/bin/bash

# Genefíor Database Builder
# Usage: ./build_database.sh -n <db_name> [-d <dna.fasta>] [-a <aa.fasta>] [-t <threads>]
#
# Flags:
#  -n, --name     Database name (REQUIRED)
#  -d, --dna      DNA FASTA file (Optional)
#  -a, --aa       AA FASTA file (if identical to the DNA file, the DNA is
#                 translated with codon table 11 to build the AA database)
#  -t, --threads  Number of threads (Default: 8)

set -e

# Default values
DB_NAME=""
DNA_FASTA=""
AA_FASTA=""
THREADS=8

# ---------------------------------------------------------------------------
# Parse Arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--name)
            DB_NAME="$2"
            shift 2
            ;;
        -d|--dna)
            DNA_FASTA="$2"
            shift 2
            ;;
        -a|--aa)
            AA_FASTA="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 -n <db_name> [-d <dna.fasta>] [-a <aa.fasta>] [-t <threads>]"
            exit 1
            ;;
    esac
done

# Validate mandatory name
if [ -z "$DB_NAME" ]; then
    echo "Error: Database name (-n) is required."
    echo "Usage: $0 -n <db_name> [-d <dna.fasta>] [-a <aa.fasta>] [-t <threads>]"
    exit 1
fi

# ---------------------------------------------------------------------------
# Helper: strip ALL extensions from a filename
# ---------------------------------------------------------------------------
strip_all_extensions() {
    local name
    name=$(basename "$1")
    name="${name%.gz}"
    while [[ "$name" == *.* ]]; do
        name="${name%.*}"
    done
    echo "$name"
}

DNA_BASE=$(strip_all_extensions "$DNA_FASTA")
AA_BASE=$(strip_all_extensions "$AA_FASTA")

# If the same file is given for DNA and AA, translate DNA (codon table 11)
TRANSLATE_AA=0
if [ -n "$DNA_FASTA" ] && [ -n "$AA_FASTA" ] && \
   [ "$(realpath "$DNA_FASTA")" = "$(realpath "$AA_FASTA")" ]; then
    TRANSLATE_AA=1
    AA_BASE="${DNA_BASE}_translated"
    echo "Same FASTA given for DNA and AA: AA database will be translated from DNA (codon table 11)."
fi

echo "Building Genefíor database: $DB_NAME"
[ -n "$DNA_FASTA" ] && echo "DNA sequences:  $DNA_FASTA  (base: $DNA_BASE)"
[ -n "$AA_FASTA" ]  && echo "AA sequences:   $AA_FASTA  (base: $AA_BASE)"
echo "Threads: $THREADS"
echo ""

# ---------------------------------------------------------------------------
# Decompress inputs and setup temp directory
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
        ln -sf "$(realpath "$src")" "$dest"
    fi
}

DNA_PLAIN="$TMPDIR_LOCAL/${DNA_BASE}.fasta"
AA_PLAIN="$TMPDIR_LOCAL/${AA_BASE}.fasta"

if [ -n "$DNA_FASTA" ]; then
    decompress "$DNA_FASTA" "$DNA_PLAIN"
fi

translate_fasta() {
    python3 - "$1" "$2" <<'PYEOF'
import sys
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
table = {a + b + c: aas[i] for i, (a, b, c) in
         enumerate((a, b, c) for a in bases for b in bases for c in bases)}

def translate(seq):
    seq = seq.upper().replace("U", "T")
    prot = "".join(table.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return prot[:-1] if prot.endswith("*") else prot

def flush(out, header, parts):
    if header is not None:
        prot = translate("".join(parts))
        out.write(header + "\n")
        for i in range(0, len(prot), 60):
            out.write(prot[i:i + 60] + "\n")

with open(sys.argv[1]) as fin, open(sys.argv[2], "w") as out:
    header, parts = None, []
    for line in fin:
        line = line.strip()
        if line.startswith(">"):
            flush(out, header, parts)
            header, parts = line, []
        elif line:
            parts.append(line)
    flush(out, header, parts)
PYEOF
}

if [ "$TRANSLATE_AA" -eq 1 ]; then
    echo "Translating DNA to protein (codon table 11)..."
    translate_fasta "$DNA_PLAIN" "$AA_PLAIN"
elif [ -n "$AA_FASTA" ]; then
    decompress "$AA_FASTA" "$AA_PLAIN"
    if ! grep -v '^>' "$AA_PLAIN" | head -n 1000 | grep -qv '^[ACGTUNacgtun[:space:]]*$'; then
        echo "Error: the -a file contains only nucleotide letters (looks like DNA)."
        echo "To translate it (codon table 11), pass the same file to both flags: -d <file> -a <file>"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Create directory structure
# ---------------------------------------------------------------------------
mkdir -p "${DB_NAME}"/{blast_aa,blast_dna,bowtie2,bwa,diamond,minimap2}

# Store compressed copies of source files for reference
if [ -n "$DNA_FASTA" ]; then
    echo "Storing DNA source file..."
    gzip -c "$DNA_PLAIN" > "${DB_NAME}/${DNA_BASE}.fasta.gz"
fi

if [ -n "$AA_FASTA" ]; then
    echo "Storing AA source file..."
    gzip -c "$AA_PLAIN"  > "${DB_NAME}/${AA_BASE}.fasta.gz"
fi

# ---------------------------------------------------------------------------
# Build BLAST DNA database
# ---------------------------------------------------------------------------
if [ -n "$DNA_FASTA" ]; then
    echo "Building BLAST DNA database..."
    makeblastdb -in "$DNA_PLAIN" -dbtype nucl \
        -out "${DB_NAME}/blast_dna/${DNA_BASE}_blastdb" \
        -parse_seqids
fi

# ---------------------------------------------------------------------------
# Build BLAST AA database & DIAMOND
# ---------------------------------------------------------------------------
if [ -n "$AA_FASTA" ]; then
    echo "Building BLAST AA database..."
    makeblastdb -in "$AA_PLAIN" -dbtype prot \
        -out "${DB_NAME}/blast_aa/${AA_BASE}_blastdb" \
        -parse_seqids

    echo "Building DIAMOND database..."
    diamond makedb --in "$AA_PLAIN" \
        --db "${DB_NAME}/diamond/${AA_BASE}_diamonddb" \
        --threads "$THREADS"
fi

# ---------------------------------------------------------------------------
# Build DNA-specific indices
# ---------------------------------------------------------------------------
if [ -n "$DNA_FASTA" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build --threads "$THREADS" "$DNA_PLAIN" \
        "${DB_NAME}/bowtie2/${DNA_BASE}_bowtie2db"

    echo "Building BWA index..."
    bwa index -p "${DB_NAME}/bwa/${DNA_BASE}_bwadb" "$DNA_PLAIN"

    echo "Building Minimap2 index..."
    minimap2 -d "${DB_NAME}/minimap2/${DNA_BASE}_minimap2db" "$DNA_PLAIN"
fi

# ---------------------------------------------------------------------------
# Validation
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

if [ -n "$DNA_FASTA" ]; then
    check_exists "${DB_NAME}/blast_dna/${DNA_BASE}_blastdb.nhr"
    check_exists "${DB_NAME}/bowtie2/${DNA_BASE}_bowtie2db.1.bt2"
    check_exists "${DB_NAME}/bwa/${DNA_BASE}_bwadb.amb"
    check_exists "${DB_NAME}/minimap2/${DNA_BASE}_minimap2db"
fi

if [ -n "$AA_FASTA" ]; then
    check_exists "${DB_NAME}/blast_aa/${AA_BASE}_blastdb.phr"
    check_exists "${DB_NAME}/diamond/${AA_BASE}_diamonddb.dmnd"
fi

echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo "Database build complete — all expected index files present."
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
