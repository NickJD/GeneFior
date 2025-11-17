#!/bin/bash

# AMRfíor Database Builder
# Usage: ./build_amrfior_database.sh database_name dna.fasta aa.fasta threads

set -e

DB_NAME=$1
DNA_FASTA=$2
AA_FASTA=$3
THREADS=${4:-4}

if [ -z "$DB_NAME" ] || [ -z "$DNA_FASTA" ] || [ -z "$AA_FASTA" ]; then
    echo "Usage: $0 <database_name> <dna.fasta> <aa.fasta> [threads]"
    exit 1
fi

echo "Building AMRfíor database: $DB_NAME"
echo "DNA sequences: $DNA_FASTA"
echo "AA sequences: $AA_FASTA"
echo "Threads: $THREADS"
echo ""

# Create directory structure
mkdir -p ${DB_NAME}/{blast_aa,blast_dna,bowtie2,bwa,diamond,minimap2}

# Copy and compress source files
echo "Copying source files..."
gzip -c $DNA_FASTA > ${DB_NAME}/$(basename "${DNA_FASTA%.*}").gz
gzip -c $AA_FASTA > ${DB_NAME}/$(basename "${AA_FASTA%.*}").gz

# Build BLAST DNA database
echo "Building BLAST DNA database..."
makeblastdb -in $DNA_FASTA -dbtype nucl \
    -out ${DB_NAME}/blast_dna/$(basename "${DNA_FASTA%.*}")_blastdb \
    -parse_seqids

# Build BLAST AA database
echo "Building BLAST AA database..."
makeblastdb -in $AA_FASTA -dbtype prot \
    -out ${DB_NAME}/blast_aa/$(basename "${AA_FASTA%.*}")_blastdb \
    -parse_seqids

# Build DIAMOND database
echo "Building DIAMOND database..."
diamond makedb --in $AA_FASTA \
    --db ${DB_NAME}/diamond/$(basename "${AA_FASTA%.*}")_diamonddb \
    --threads $THREADS

# Build Bowtie2 index
echo "Building Bowtie2 index..."
bowtie2-build --threads $THREADS $DNA_FASTA \
    ${DB_NAME}/bowtie2/$(basename "${DNA_FASTA%.*}")_bowtie2db

# Build BWA index
echo "Building BWA index..."
bwa index -p ${DB_NAME}/bwa/$(basename "${DNA_FASTA%.*}")_bwadb $DNA_FASTA

# Build Minimap2 index
echo "Building Minimap2 index..."
minimap2 -d ${DB_NAME}/minimap2/$(basename "${DNA_FASTA%.*}")_minimap2db $DNA_FASTA

echo ""
echo "Database build complete!"
echo "Database location: $DB_NAME"
echo ""
echo "To use with AMRfíor:"
echo "AMRfior -i reads.fasta -st Single-FASTA -o results/ \\"
echo "  --databases user-provided \\"
echo "  --user-db-path $(pwd)/$DB_NAME"