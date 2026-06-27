from pathlib import Path

import pytest

from GeneFior.databases import (
    gather_multiple_databases,
    parse_database_paths,
)
from GeneFior.utils import diamond_mode_label


def _make_diamond_database(root):
    diamond_dir = root / 'diamond'
    diamond_dir.mkdir(parents=True)
    (diamond_dir / 'test_diamonddb.dmnd').write_text('database')


def test_parse_database_paths_accepts_comma_separated_directories(tmp_path):
    first = tmp_path / 'first'
    second = tmp_path / 'second'

    paths = parse_database_paths(f" {first}, {second} ")

    assert paths == [str(first.resolve()), str(second.resolve())]


def test_single_database_keeps_legacy_label(tmp_path):
    database = tmp_path / 'custom_amr'
    _make_diamond_database(database)

    gathered = gather_multiple_databases(str(database), ['diamond'])

    assert list(gathered) == ['user-provided-db']
    assert gathered['user-provided-db']['diamond'].endswith(
        'test_diamonddb.dmnd'
    )


def test_multiple_databases_use_directory_names(tmp_path):
    first = tmp_path / 'environmental_amr'
    second = tmp_path / 'virulence'
    _make_diamond_database(first)
    _make_diamond_database(second)

    gathered = gather_multiple_databases(
        f"{first},{second}", ['diamond']
    )

    assert list(gathered) == ['environmental_amr', 'virulence']
    assert all(mapping.get('diamond') for mapping in gathered.values())


def test_duplicate_database_directory_names_are_rejected(tmp_path):
    first = tmp_path / 'one' / 'database'
    second = tmp_path / 'two' / 'database'
    _make_diamond_database(first)
    _make_diamond_database(second)

    with pytest.raises(ValueError, match='same directory name'):
        gather_multiple_databases(f"{first},{second}", ['diamond'])


@pytest.mark.parametrize(
    'sequence_type,genes_type,expected',
    [
        ('Paired-FASTQ', None, 'DIAMOND-BLASTX'),
        ('Single-FASTA', None, 'DIAMOND-BLASTX'),
        ('Genes-FASTA', 'dna', 'DIAMOND-BLASTX'),
        ('Genes-FASTA', 'aa', 'DIAMOND-BLASTP'),
        ('Genes-FASTA', 'protein', 'DIAMOND-BLASTP'),
    ],
)
def test_diamond_mode_label(sequence_type, genes_type, expected):
    assert diamond_mode_label(sequence_type, genes_type) == expected
