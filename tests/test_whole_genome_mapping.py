import csv
import logging

from GeneFior.evidence import (
    EvidenceConfig,
    WHOLE_GENOME_COMPLETE,
    WHOLE_GENOME_NEAR_COMPLETE,
    WHOLE_GENOME_PARTIAL,
)
from GeneFior.gene_stats import GeneStats
from GeneFior.workflow import Workflow, split_whole_genome_reference


def _stats(name, length=100, depth=3, identity=100.0, end=None):
    stats = GeneStats(name)
    for _ in range(depth):
        stats.add_hit(1, end or length, identity, length)
    stats.add_read_support(mapped=True, passing=True, best=True)
    stats.finalise()
    return stats


def test_whole_genome_reference_names_split_only_on_explicit_contig_forms():
    assert split_whole_genome_reference('genome_A|contig_1') == (
        'genome_A', 'contig_1'
    )
    assert split_whole_genome_reference('genome_A::scaffold-2') == (
        'genome_A', 'scaffold-2'
    )
    assert split_whole_genome_reference('genome_A_contig_1') == (
        'genome_A', 'contig_1'
    )
    assert split_whole_genome_reference('unstructured_reference_1') == (
        'unstructured_reference_1', 'unstructured_reference_1'
    )


def test_genome_summary_aggregates_contigs_and_keeps_genomes_separate(tmp_path):
    workflow = Workflow.__new__(Workflow)
    workflow.output_dir = tmp_path
    workflow.logger = logging.getLogger('test-whole-genome-mapping')
    workflow.gene_stats = {
        'db': {
            'minimap2': {
                'genome_A_contig_1': _stats('genome_A_contig_1'),
                'genome_A_contig_2': _stats('genome_A_contig_2'),
                'genome_B_contig_1': _stats('genome_B_contig_1', end=90),
                'genome_C_contig_1': _stats('genome_C_contig_1', end=30),
            }
        }
    }
    workflow.evidence_calls = {'db': {}}
    workflow.evidence_config = EvidenceConfig()
    workflow.is_genes_fasta = False
    workflow.detection_min_coverage = 80.0
    workflow.detection_min_identity = 80.0
    workflow.detection_min_num_reads = 1

    workflow.generate_genome_mapping_summary('db')

    with open(tmp_path / 'db_genome_mapping_summary.tsv', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert {row['Genome_ID'] for row in rows} == {'genome_A', 'genome_B', 'genome_C'}
    genome_a = next(row for row in rows if row['Genome_ID'] == 'genome_A')
    assert genome_a['Contig_Count'] == '2'
    assert genome_a['Reference_Length'] == '200'
    assert genome_a['Coverage_3x'] == '100.00'
    assert genome_a['Mean_Depth'] == '3.00'
    assert genome_a['Evidence_Status'] == WHOLE_GENOME_COMPLETE
    genome_b = next(row for row in rows if row['Genome_ID'] == 'genome_B')
    assert genome_b['Evidence_Status'] == WHOLE_GENOME_NEAR_COMPLETE
    genome_c = next(row for row in rows if row['Genome_ID'] == 'genome_C')
    assert genome_c['Evidence_Status'] == WHOLE_GENOME_PARTIAL

    with open(tmp_path / 'db_contig_mapping_summary.tsv', newline='') as handle:
        contigs = list(csv.DictReader(handle, delimiter='\t'))
    assert len(contigs) == 4
    assert {row['Genome_ID'] for row in contigs} == {'genome_A', 'genome_B', 'genome_C'}
