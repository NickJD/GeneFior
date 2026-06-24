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

    # ---------- HMMER early discovery (handled separately, not via tool_tags loop) ----------
    # Accept any of the common directory names for protein/dna HMM databases.
    # Also detect any annotations CSV in the same directory.
    hmmer_flat: Dict[str, str] = {}
    base_dir_path = Path(base_dir)
    _hmmer_protein_dirs = ('hmmer', 'hmmer_protein', 'hmm', 'hmm_protein')
    _hmmer_dna_dirs = ('hmmer_dna', 'nhmmer', 'hmm_dna')
    for _dirname, _key in (
        *[(_d, 'hmmer_protein') for _d in _hmmer_protein_dirs],
        *[(_d, 'hmmer_dna') for _d in _hmmer_dna_dirs],
    ):
        _candidate = base_dir_path / _dirname
        if _candidate.is_dir():
            _hmm_files = list(_candidate.glob('*.hmm'))
            if _hmm_files:
                hmmer_flat[_key] = str(_hmm_files[0])
                # Look for an annotations CSV in the same directory
                _csv_files = list(_candidate.glob('*.csv'))
                if _csv_files and 'hmmer_annotations' not in hmmer_flat:
                    hmmer_flat['hmmer_annotations'] = str(_csv_files[0])
            break  # only use first matching directory per mode
    # Also check the root base_dir itself for a .hmm file (flat layout)
    if 'hmmer_protein' not in hmmer_flat and ('hmmer_protein' in (tools or []) or 'all' in (tools or [])):
        _root_hmm = list(base_dir_path.glob('*.hmm'))
        if _root_hmm:
            hmmer_flat['hmmer_protein'] = str(_root_hmm[0])
            _root_csv = list(base_dir_path.glob('*.csv'))
            if _root_csv and 'hmmer_annotations' not in hmmer_flat:
                hmmer_flat['hmmer_annotations'] = str(_root_csv[0])
    # ---------- end HMMER discovery ----------
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
                # Special handling for bwa category.
                # BWA index files share a common base name with single-part extensions
                # (.amb, .ann, .bwt, .pac, .sa).  Use .amb (always present) so
                # with_suffix('') reliably yields the correct basename, e.g.:
                #   VFDB_setA_nt.fas_bwadb.amb  ->  VFDB_setA_nt.fas_bwadb
                for tool, tag in tool_tags.items():
                    if tool == "bwa":
                        # Prefer .amb as the anchor file; fall back to any bwadb file
                        anchor_files = list(category_dir.rglob(f"*{tag}.amb"))
                        if not anchor_files:
                            anchor_files = sorted(category_dir.rglob(f"*{tag}*"))
                        if anchor_files:
                            databases[category_name][tool] = str(
                                Path(anchor_files[0]).with_suffix('')
                            )
                continue
            if category_name == "bowtie2":
                # Special handling for bowtie2 category.
                # Bowtie2 index files are named <base>.1.bt2, <base>.2.bt2,
                # <base>.rev.1.bt2, etc.  The previous code used split('.')[0]
                # which incorrectly discarded everything after the first dot in the
                # filename (e.g. VFDB_setA_nt.fas_bowtie2db -> VFDB_setA_nt).
                # Instead, locate the .1.bt2 primary index file and strip that
                # known suffix to recover the exact basename bowtie2 expects.
                for tool, tag in tool_tags.items():
                    if tool == "bowtie2":
                        # Primary index always ends in <tag>.1.bt2
                        primary_files = list(category_dir.rglob(f"*{tag}.1.bt2"))
                        if not primary_files:
                            # Fallback: any .bt2 file – derive base via tag position
                            primary_files = sorted(category_dir.rglob(f"*{tag}*.bt2"))
                        if primary_files:
                            base = str(primary_files[0])
                            if base.endswith('.1.bt2'):
                                # Clean strip: remove the '.1.bt2' suffix
                                databases[category_name][tool] = base[:-len('.1.bt2')]
                            else:
                                # Fallback: find the tag in the path and cut there
                                tag_end = base.rfind(tag)
                                if tag_end != -1:
                                    databases[category_name][tool] = base[:tag_end + len(tag)]
                continue
            if category_name == "minimap2":
                # Special handling for minimap2 category.
                # A minimap2 index is a single file (e.g. VFDB_setA_nt.fas_minimap2db or .mmi).
                # Do NOT call with_suffix('') – if the file is named with a compound tag
                # like ".fas_minimap2db", pathlib treats the whole thing as the suffix and
                # with_suffix('') would strip it back to just "VFDB_setA_nt".
                # Pass the full file path as minimap2 accepts the index file directly.
                for tool, tag in tool_tags.items():
                    if tool == "minimap2":
                        matching_files = list(category_dir.rglob(f"*{tag}"))
                        if not matching_files:
                            # Also accept .mmi extension (standard minimap2 index extension)
                            matching_files = list(category_dir.rglob(f"*{tag}.mmi"))
                        if matching_files:
                            databases[category_name][tool] = str(matching_files[0])
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

    # Merge HMMER results (only keys not already set by the loop above)
    for k, v in hmmer_flat.items():
        if k not in flat and v:
            flat[k] = v

    return flat
