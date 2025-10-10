from typing import List, Set
from dataclasses import dataclass, field
import sys


@dataclass
class GeneStats:
    """Statistics for a gene detection.

    Only sequences meeting min_identity threshold are counted.
    - gene_coverage: Percentage of the gene/subject sequence covered by all alignments combined
    - avg_identity: Average identity across all qualifying sequences for this gene
    """
    gene_name: str
    gene_length: int = 0  # Length of the gene in the database
    num_sequences: int = 0  # Number of sequences passing min_identity
    gene_coverage: float = 0.0  # Percentage of gene covered by all alignments
    avg_identity: float = 0.0
    identities: List[float] = field(default_factory=list)
    covered_positions: Set[int] = field(default_factory=set)  # All positions covered on the gene

    def add_hit(self, sstart: int, send: int, identity: float, gene_len: int = 0):
        """Add a hit to the statistics (only called for sequences passing min_identity).

        Args:
            sstart: Start position on subject/gene (1-based, inclusive)
            send: End position on subject/gene (1-based, inclusive)
            identity: Percent identity of this alignment
            gene_len: Length of the gene (if known)
        """
        self.num_sequences += 1
        self.identities.append(identity)

        # Add all positions covered by this alignment (convert to 0-based for consistency)
        start = min(sstart, send) - 1  # Convert to 0-based
        end = max(sstart, send)  # Inclusive end in 1-based becomes exclusive in 0-based
        for pos in range(start, end):
            self.covered_positions.add(pos)

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

    def add_positions(self, positions: Set[int], identity: float, gene_len: int = 0):
        """Add a set of positions directly (for BAM parsing where we track exact positions).

        Args:
            positions: Set of 0-based positions covered on the gene
            identity: Percent identity of this alignment
            gene_len: Length of the gene (if known)
        """
        self.num_sequences += 1
        self.identities.append(identity)
        self.covered_positions.update(positions)

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

    def finalise(self):
        """Calculate final statistics."""
        if self.num_sequences > 0:
            self.avg_identity = sum(self.identities) / self.num_sequences

        # Calculate gene coverage as percentage of gene length covered
        if self.gene_length > 0:
            self.gene_coverage = (len(self.covered_positions) / self.gene_length) * 100


# @dataclass
# class GeneStats:
#     """Statistics for a gene detection.
#     Only sequences meeting min_identity threshold are counted.
#     - gene_coverage: Percentage of the gene/subject sequence covered by all alignments combined
#     - avg_identity: Average identity across all qualifying sequences for this gene
#     """
#     gene_name: str
#     gene_length: int = 0  # Length of the gene in the database
#     num_sequences: int = 0  # Number of sequences passing min_identity
#     gene_coverage: float = 0.0  # Percentage of gene covered by all alignments
#     avg_identity: float = 0.0
#     identities: List[float] = field(default_factory=list)
#     covered_positions: Set[int] = field(default_factory=set)  # All positions covered on the gene
#
#     def add_hit(self, sstart: int, send: int, identity: float, gene_len: int = 0):
#         """Add a hit to the statistics (only called for sequences passing min_identity).
#         Args:
#             sstart: Start position on subject/gene
#             send: End position on subject/gene
#             identity: Percent identity of this alignment
#             gene_len: Length of the gene (if known)
#         """
#         self.num_sequences += 1
#         self.identities.append(identity)
#
#         # Add all positions covered by this alignment
#         start = min(sstart, send)
#         end = max(sstart, send)
#         for pos in range(start, end + 1):
#             self.covered_positions.add(pos)
#
#         if gene_len > 0:
#             self.gene_length = max(self.gene_length, gene_len)
#         else:
#             sys.exit("gene length not known")
#
#     def finalise(self):
#         """Calculate final statistics."""
#         if self.num_sequences > 0:
#             self.avg_identity = sum(self.identities) / self.num_sequences
#
#         # Calculate gene coverage as percentage of gene length covered
#         if self.gene_length > 0:
#             self.gene_coverage = (len(self.covered_positions) / self.gene_length) * 100
#
