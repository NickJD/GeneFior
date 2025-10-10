"""Database path configuration for AMR detection tools."""

from pathlib import Path

# Get the directory where this file is located
PACKAGE_DIR = Path(__file__).parent
DB_ROOT = PACKAGE_DIR / "databases"

CARD_DATABASES = {
    "diamond": str(DB_ROOT / "card/diamond/protein_fasta_protein_homolog_model_SID_diamonddb.dmnd"),
    "blastn": str(DB_ROOT / "card/blast_dna/nucleotide_fasta_protein_homolog_model_SID_blastdb"),
    "blastx": str(DB_ROOT / "card/blast_aa/protein_fasta_protein_homolog_model_SID_blastdb"),
    "bowtie2": str(DB_ROOT / "card/bowtie2/nucleotide_fasta_protein_homolog_model_SID_bowtie2db"),
    "bwa": str(DB_ROOT / "card/bwa/nucleotide_fasta_protein_homolog_model_SID_bwadb"),
    "minimap2": str(DB_ROOT / "card/minimap2/nucleotide_fasta_protein_homolog_model_SID_minimap2db"),
}

RESFINDER_DATABASES = {
    "diamond": str(DB_ROOT / "resfinder/diamond/all_aa_diamonddb.dmnd"),
    "blastn": str(DB_ROOT / "resfinder/blast_dna/all_blastdb"),
    "blastx": str(DB_ROOT / "resfinder/blast_aa/all_aa_blastdb"),
    "bowtie2": str(DB_ROOT / "resfinder/bowtie2/all_bowtie2db"),
    "bwa": str(DB_ROOT / "resfinder/bwa/all_bwadb"),
    "minimap2": str(DB_ROOT / "resfinder/minimap2/all_minimap2db"),
}