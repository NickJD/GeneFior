# Database path configuration for AMR detection tools.

from pathlib import Path
from typing import Dict

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

NCBI_DATABASES = {
    "diamond": str(DB_ROOT / "ncbi/diamond/sequence_aa_diamonddb.dmnd"),
    "blastn": str(DB_ROOT / "ncbi/blast_dna/sequence_dna_blastdb"),
    "blastx": str(DB_ROOT / "ncbi/blast_aa/sequence_aa_blastdb"),
    "bowtie2": str(DB_ROOT / "ncbi/bowtie2/sequence_dna_bowtie2db"),
    "bwa": str(DB_ROOT / "ncbi/bwa/sequence_dna_bwadb"),
    "minimap2": str(DB_ROOT / "ncbi/minimap2/sequence_dna_minimap2db"),
}


def gather_databases(base_dir: Path, tools: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    # Automatically gathers database files from a given directory based on unique tags for each tool.
    # Define the unique tags for each tool
    tool_tags = {
        "diamond": "diamonddb.dmnd",
        "blastn": "blastdb",
        "blastx": "blastdb",
        "bowtie2": "bowtie2db",
        "bwa": "bwadb",
        "minimap2": "minimap2db",
    }
    tool_tags = {tool: tool_tags[tool] for tool in tools if tool in tool_tags}  # Filter tags based on provided tools
    # Initialise the result dictionary
    databases = {}
    # Iterate over subdirectories in the base directory
    base_dir = Path(base_dir)  # Ensure base_dir is a Path object
    for category_dir in base_dir.iterdir():
        if category_dir.is_dir():
            category_name = category_dir.name
            databases[category_name] = {}
            if category_name == "blast_dna":
                # Special handling for blast_dna category
                for tool, tag in tool_tags.items():
                    if tool == "blastn":
                        matching_files = list(category_dir.rglob(f"*{tag}*"))
                        nucleotide_files = [f for f in matching_files if "dna" in f.parts or "blast_dna" in f.parts]
                        if nucleotide_files:
                            databases[category_name][tool] = str(nucleotide_files[0].with_suffix(''))
                continue
            if category_name == "blast_aa":
                # Special handling for blast_aa category
                # Expose discovered protein BLAST databases as both blastx (for translated searches)
                # and blastp (for direct protein searches) so downstream code can prefer blastp when
                # protein inputs are provided.
                matching_files = list(category_dir.rglob(f"*{tool_tags.get('blastx','blastdb')}*"))
                amino_acid_files = [f for f in matching_files if "aa" in f.parts or "blast_aa" in f.parts]
                if amino_acid_files:
                    db_base = str(amino_acid_files[0].with_suffix(''))
                    databases[category_name]['blastx'] = db_base
                    # Also expose the same base for blastp (protein BLAST) so the workflow can call BLASTp
                    databases[category_name]['blastp'] = db_base
                continue
            if category_name == "diamond":
                # Special handling for diamond category
                for tool, tag in tool_tags.items():
                    if tool == "diamond":
                        matching_files = list(category_dir.rglob(f"*{tag}"))
                        if matching_files:
                            databases[category_name][tool] = str(matching_files[0])
                continue
            if category_name == "bwa":
                # Special handling for bwa category
                for tool, tag in tool_tags.items():
                    if tool == "bwa":
                        matching_files = list(category_dir.rglob(f"*{tag}*"))
                        if matching_files:
                            databases[category_name][tool] = str(matching_files[0].with_suffix(''))
                continue
            if category_name == "bowtie2":
                # Special handling for bowtie2 category
                for tool, tag in tool_tags.items():
                    if tool == "bowtie2":
                        matching_files = list(category_dir.rglob(f"*{tag}*"))
                        if matching_files: # Not clean
                            matching_file = str(matching_files[0].with_suffix('')).split('.')[0]
                            databases[category_name][tool] = matching_file
                continue
            if category_name == "minimap2":
                # Special handling for minimap2 category
                for tool, tag in tool_tags.items():
                    if tool == "minimap2":
                        matching_files = list(category_dir.rglob(f"*{tag}"))
                        if matching_files:
                            databases[category_name][tool] = str(matching_files[0].with_suffix(''))
                continue

    # Normalise category keys to expected names if present
    if "blast_dna" in databases:
        databases["blastn"] = databases.pop("blast_dna")
    if "blast_aa" in databases:
        databases["blastx"] = databases.pop("blast_aa")

    # Flatten nested mapping into tool -> path mapping so callers receive a consistent
    # dict[str, str] (matching the shape used for packaged DB constants).
    flat = {}
    for category, mapping in databases.items():
        if isinstance(mapping, dict):
            for tool_key, path_val in mapping.items():
                # Only set if we have a truthy path
                if path_val:
                    flat[tool_key] = path_val

    return flat
