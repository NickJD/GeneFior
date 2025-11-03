from typing import List, Set
from dataclasses import dataclass, field
from typing import  Dict


@dataclass
class GeneStats:
    """Statistics for a gene detection.

    Only sequences meeting thresholds are counted.
    - gene_coverage: Percentage of the gene/subject sequence covered by all alignments combined
    - avg_identity: Average identity across all qualifying sequences for this gene
    - base_depth: Average depth across the gene (including uncovered positions as 0)
    - base_depth_hit: Average depth across covered positions only
    0-based positions are used throughout.
    - gene_length: Length of the gene in the database
    - num_sequences: Number of sequences passing the identity threshold
    - identities: List of individual sequence identities
    - covered_positions: Set of all positions covered on the gene (0-based)
    - position_depths: Dictionary tracking depth at each position on the gene (0-based)
    """
    gene_name: str
    gene_length: int = 0  # Length of the gene in the database
    num_sequences: int = 0  # Number of sequences passing min_identity
    gene_coverage: float = 0.0  # Percentage of gene covered by all alignments
    base_depth: float = 0.0 # Average 'depth' across the gene
    base_depth_hit: float = 0.0 # Average 'depth' across covered positions only
    avg_identity: float = 0.0 # Average identity across all qualifying sequences
    identities: List[float] = field(default_factory=list) # List of individual sequence identities
    covered_positions: Set[int] = field(default_factory=set)  # All positions covered on the gene
    position_depths: Dict[int, int] = field(default_factory=dict) # Track depth per position


    def add_hit(self, sstart: int, send: int, identity: float, gene_len: int = 0):
        # Add a hit to the statistics (only called for sequences passing min_identity).

        self.num_sequences += 1
        self.identities.append(identity)

        # Add all positions covered by this alignment (convert to 0-based for consistency)
        start = min(sstart, send) - 1  # Convert to 0-based
        end = max(sstart, send)  # Inclusive end in 1-based becomes exclusive in 0-based
        for pos in range(start, end):
            self.covered_positions.add(pos)
            # Track depth at each position
            self.position_depths[pos] = self.position_depths.get(pos, 0) + 1

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

    def add_positions(self, positions: Set[int], identity: float, gene_len: int = 0):
        # Add a set of positions directly (for BAM parsing).
        self.num_sequences += 1
        self.identities.append(identity)
        self.covered_positions.update(positions)

        # Track depth at each position
        for pos in positions:
            self.position_depths[pos] = self.position_depths.get(pos, 0) + 1

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

    def finalise(self):
        # Calculate final statistics.
        if self.num_sequences > 0:
            self.avg_identity = sum(self.identities) / self.num_sequences

        # Calculate gene coverage as percentage of gene length covered
        if self.gene_length > 0:
            self.gene_coverage = (len(self.covered_positions) / self.gene_length) * 100
            # Calculate base depth (average depth across covered positions)
            if self.position_depths:
                total_depth = sum(self.position_depths.values())
                # Average depth across ALL gene positions (including uncovered = 0)
                self.base_depth = total_depth / self.gene_length
                # Average depth across positions with at least one read mapped
                covered_depth = sum(depth for pos, depth in self.position_depths.items() if depth > 0)
                self.base_depth_hit = covered_depth / len(self.position_depths)

