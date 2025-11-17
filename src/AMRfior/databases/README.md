
# AMRfíor Custom Database Setup Guide

This guide explains how to prepare your own antimicrobial resistance gene (or any gene-collection) database for use with AMRfíor.

## Directory Structure Overview

AMRfíor expects databases to follow a specific directory structure. Each database should be organised as follows:

```
your_database_name/
├── README                          # Optional: Database description
├── sequence_aa.fasta.gz           # Source amino acid sequences (optional)
├── sequence_dna.fasta.gz          # Source DNA sequences (optional)
├── blast_aa/                      # BLAST protein database
│   └── sequence_aa_blastdb.*
├── blast_dna/                     # BLAST nucleotide database
│   └── sequence_dna_blastdb.*
├── bowtie2/                       # Bowtie2 index files
│   └── sequence_dna_bowtie2db.*
├── bwa/                           # BWA index files
│   └── sequence_dna_bwadb.*
├── diamond/                       # DIAMOND database
│   └── sequence_aa_diamonddb.dmnd
└── minimap2/                      # Minimap2 index
    └── sequence_dna_minimap2db
```

## Prerequisites

Ensure the following tools are installed and accessible in your PATH:
- `makeblastdb` (BLAST+)
- `diamond`
- `bowtie2-build`
- `bwa`
- `minimap2`

## Step-by-Step Database Creation

### 1. Prepare Your Source Sequences

Create two FASTA files containing your AMR gene sequences:

- **sequence_dna.fasta** - Nucleotide sequences
- **sequence_aa.fasta** - Amino acid sequences (translated ORFs)

**Optional:** Compress them with gzip to save space:
```bash
gzip sequence_dna.fasta
gzip sequence_aa.fasta
```

### 2. Create Directory Structure

```bash
mkdir -p your_database_name/{blast_aa,blast_dna,bowtie2,bwa,diamond,minimap2}
```

### 3. Build Tool-Specific Databases

#### BLAST Nucleotide Database (DNA-based searches)

```bash

makeblastdb \
    -in sequence_dna.fasta \
    -dbtype nucl \
    -out sequence_dna_blastdb 

```

**Expected output files:**
- `sequence_dna_blastdb.ndb`
- `sequence_dna_blastdb.nhr`
- `sequence_dna_blastdb.nin`
- `sequence_dna_blastdb.not`
- `sequence_dna_blastdb.nsq`
- `sequence_dna_blastdb.ntf`
- `sequence_dna_blastdb.nto`

#### BLAST Protein Database (Protein-based searches)

```bash
makeblastdb \
    -in sequence_aa.fasta\
    -dbtype prot \
    -out sequence_aa_blastdb 

```

**Expected output files:**
- `sequence_aa_blastdb.pdb`
- `sequence_aa_blastdb.phr`
- `sequence_aa_blastdb.pin`
- `sequence_aa_blastdb.psq`
- `sequence_aa_blastdb.pot`
- `sequence_aa_blastdb.ptf`
- `sequence_aa_blastdb.pto`

#### DIAMOND Database

```bash
diamond makedb \
    --in sequence_aa.fasta \
    --db sequence_aa_diamonddb \

```

**Expected output:** `sequence_aa_diamonddb.dmnd`

#### Bowtie2 Index

```bash

bowtie2-build \
    sequence_dna.fasta \
    sequence_dna_bowtie2db


```

**Expected output files:**
- `sequence_dna_bowtie2db.1.bt2`
- `sequence_dna_bowtie2db.2.bt2`
- `sequence_dna_bowtie2db.3.bt2`
- `sequence_dna_bowtie2db.4.bt2`
- `sequence_dna_bowtie2db.rev.1.bt2`
- `sequence_dna_bowtie2db.rev.2.bt2`

#### BWA Index

```bash
bwa index \
    -p sequence_dna_bwadb \
    sequence_dna.fasta

```

**Expected output files:**
- `sequence_dna_bwadb.amb`
- `sequence_dna_bwadb.ann`
- `sequence_dna_bwadb.bwt`
- `sequence_dna_bwadb.pac`
- `sequence_dna_bwadb.sa`

#### Minimap2 Index

```bash
minimap2 \
    -d sequence_dna_minimap2db \
    sequence_dna.fasta

```

**Expected output:** `sequence_dna_minimap2db`

## 4. Verify Your Database Structure

After building all databases, verify the structure:

```bash
tree your_database_name/
```

Or check manually:
```bash
ls -R your_database_name/
```

Ensure all expected files are present in their respective directories.

## 5. Using Your Custom Database with AMRfíor

### Option A (Receommended!!): Use with `--user-db-path`

Keep your database anywhere and specify the path at runtime:

```bash
AMRfior \
    -i reads.fasta \
    -st Single-FASTA \
    -o results/ \
    --databases user-provided \
    --user-db-path /path/to/your_database_name
```

### Option B: Install in AMRfíor Package Directory

Copy your database to AMRfíor's database directory:

```bash
cp -r your_database_name /path/to/AMRfíor/databases/
```

Then update `databases.py` to register your database (see "Registering Custom Databases" below).



## Registering Custom Databases (Option B)

If you want your database permanently available, edit `src/AMRfíor/databases.py`:

```python
YOUR_DATABASE_NAME = {
    "blastn": "databases/your_database_name/blast_dna/sequence_dna_blastdb",
    "blastx": "databases/your_database_name/blast_aa/sequence_aa_blastdb",
    "diamond": "databases/your_database_name/diamond/sequence_aa_diamonddb.dmnd",
    "bowtie2": "databases/your_database_name/bowtie2/sequence_dna_bowtie2db",
    "bwa": "databases/your_database_name/bwa/sequence_dna_bwadb",
    "minimap2": "databases/your_database_name/minimap2/sequence_dna_minimap2db",
}
```

Then use it:
```bash
AMRfior -i reads.fasta -st Single-FASTA -o results/ --databases your_database_name
```

## Naming Conventions (IMPORTANT!!!)

AMRfíor expects specific naming patterns. **Prefix is YOUR DNA/AA filename - SUFFIX IS NON-NEGOTIABLE - Do not change!:**

| Tool | Prefix | Suffix |
|------|--------|--------|
| BLAST (DNA) | `sequence_dna_` | `blastdb` |
| BLAST (AA) | `sequence_aa_` | `blastdb` |
| DIAMOND | `sequence_aa_` | `diamonddb.dmnd` |
| Bowtie2 | `sequence_dna_` | `bowtie2db` |
| BWA | `sequence_dna_` | `bwadb` |
| Minimap2 | `sequence_dna_` | `minimap2db` |

### Complete Automation Script - 'Should' work

Save this as `build_amrfior_database.sh`:

```bash
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
gzip -c $DNA_FASTA > ${DB_NAME}/$(basename "$DNA_FASTA").fasta.gz
gzip -c $AA_FASTA > ${DB_NAME}/$(basename "$DNA_FASTA").fasta.gz

# Build BLAST DNA database
echo "Building BLAST DNA database..."
makeblastdb -in $DNA_FASTA -dbtype nucl \
    -out ${DB_NAME}/blast_dna/$(basename "$DNA_FASTA")_blastdb \
    -parse_seqids

# Build BLAST AA database
echo "Building BLAST AA database..."
makeblastdb -in $AA_FASTA -dbtype prot \
    -out ${DB_NAME}/blast_aa/$(basename "$DNA_FASTA")_blastdb \
    -parse_seqids

# Build DIAMOND database
echo "Building DIAMOND database..."
diamond makedb --in $AA_FASTA \
    --db ${DB_NAME}/diamond/$(basename "$DNA_FASTA")_diamonddb \
    --threads $THREADS

# Build Bowtie2 index
echo "Building Bowtie2 index..."
bowtie2-build --threads $THREADS $DNA_FASTA \
    ${DB_NAME}/bowtie2/$(basename "$DNA_FASTA")_bowtie2db

# Build BWA index
echo "Building BWA index..."
bwa index -p ${DB_NAME}/bwa/$(basename "$DNA_FASTA")_bwadb $DNA_FASTA

# Build Minimap2 index
echo "Building Minimap2 index..."
minimap2 -d ${DB_NAME}/minimap2/$(basename "$DNA_FASTA")_minimap2db $DNA_FASTA

echo ""
echo "Database build complete!"
echo "Database location: $DB_NAME"
echo ""
echo "To use with AMRfíor:"
echo "AMRfior -i reads.fasta -st Single-FASTA -o results/ \\"
echo "  --databases user-provided \\"
echo "  --user-db-path $(pwd)/$DB_NAME"
```

Make it executable:
```bash
chmod +x build_amrfior_database.sh
```

Run it:
```bash
./build_amrfior_database.sh my_custom_db genes_dna.fasta genes_aa.fasta 8
```

## Troubleshooting

**Error: "Database not found"**
- Check directory structure matches exactly
- Verify file naming conventions
- Ensure database path is correct

**Error: "makeblastdb not found"**
- Install BLAST+: `conda install -c bioconda blast`

**Error: "diamond not found"**
- Install DIAMOND: `conda install -c bioconda diamond`

**Large database memory issues**
- For very large databases (>10GB), consider splitting into smaller subsets
- Increase system memory or use a compute cluster

## Best Practices

1. **Version control**: Include version/date in README file
2. **Validation**: Test your database with a known positive control
3. **Documentation**: Document your sequence sources and curation methods
4. **Backup**: Keep source FASTA files separate from indices
5. **Updates**: Rebuild all indices (databases) when updating source sequences

---

For questions or issues, please open an issue on the AMRfíor GitHub repository.
