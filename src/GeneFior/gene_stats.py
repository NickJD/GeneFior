import math
from typing import Dict, Iterator, List, Set, Tuple
from dataclasses import dataclass, field


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
    - identity_sum: Running sum of individual sequence identities
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
    identity_sum: float = 0.0 # Running identity total; avoids storing one value per hit
    identities: List[float] = field(default_factory=list) # Deprecated compatibility field; no longer populated
    covered_positions: Set[int] = field(default_factory=set)  # All positions covered on the gene
    position_depths: Dict[int, int] = field(default_factory=dict) # Track depth per position
    compact: bool = False # Store interval boundary deltas instead of one entry per covered base
    _depth_deltas: Dict[int, int] = field(default_factory=dict, repr=False)
    coverage_2x: float = 0.0
    coverage_3x: float = 0.0
    coverage_5x: float = 0.0
    coverage_10x: float = 0.0
    median_depth: float = 0.0
    depth_cv: float = 0.0
    num_gaps: int = 0
    num_internal_gaps: int = 0
    longest_gap: int = 0
    longest_internal_gap: int = 0
    longest_covered_run: int = 0
    mapped_read_support: int = 0
    passing_read_support: int = 0
    best_read_support: int = 0
    unique_best_read_support: int = 0
    ambiguous_best_read_support: int = 0
    high_confidence_read_support: int = 0
    support_score_sum: float = 0.0
    _coverage_by_depth: Dict[int, float] = field(default_factory=dict, repr=False)

    def _add_interval(self, start: int, end: int):
        """Add a 0-based half-open covered interval."""
        start = max(0, int(start))
        end = int(end)
        if self.gene_length > 0:
            end = min(end, self.gene_length)
        if end <= start:
            return
        if self.compact:
            self._depth_deltas[start] = self._depth_deltas.get(start, 0) + 1
            self._depth_deltas[end] = self._depth_deltas.get(end, 0) - 1
            return
        for pos in range(start, end):
            self.covered_positions.add(pos)
            self.position_depths[pos] = self.position_depths.get(pos, 0) + 1


    def add_hit(self, sstart: int, send: int, identity: float, gene_len: int = 0):
        # Add a hit to the statistics (only called for sequences passing min_identity).

        self.num_sequences += 1
        self.identity_sum += identity

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

        # Add all positions covered by this alignment (convert to 0-based for consistency)
        start = min(sstart, send) - 1  # Convert to 0-based
        end = max(sstart, send)  # Inclusive end in 1-based becomes exclusive in 0-based
        self._add_interval(start, end)

    def add_positions(self, positions: Set[int], identity: float, gene_len: int = 0):
        # Add a set of positions directly (for BAM parsing).
        self.num_sequences += 1
        self.identity_sum += identity

        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)
        if self.compact:
            sorted_positions = sorted(positions)
            if not sorted_positions:
                return
            start = previous = sorted_positions[0]
            for pos in sorted_positions[1:]:
                if pos != previous + 1:
                    self._add_interval(start, previous + 1)
                    start = pos
                previous = pos
            self._add_interval(start, previous + 1)
            return

        self.covered_positions.update(positions)
        for pos in positions:
            self.position_depths[pos] = self.position_depths.get(pos, 0) + 1

    def add_intervals(self, intervals: List[Tuple[int, int]], identity: float,
                      gene_len: int = 0):
        """Add one read represented by 0-based half-open reference intervals."""
        self.num_sequences += 1
        self.identity_sum += identity
        if gene_len > 0:
            self.gene_length = max(self.gene_length, gene_len)

        normalised_intervals = []
        for start, end in intervals:
            start = max(0, int(start))
            end = int(end)
            if self.gene_length > 0:
                end = min(end, self.gene_length)
            if end > start:
                normalised_intervals.append((start, end))
        if not normalised_intervals:
            return

        # add_positions() deduplicates positions within one read via a set.
        # Merge here so the compact interval path reports identical depth
        # metrics even if a caller supplies overlapping intervals.
        normalised_intervals.sort()
        merged_start, merged_end = normalised_intervals[0]
        for start, end in normalised_intervals[1:]:
            if start <= merged_end:
                merged_end = max(merged_end, end)
                continue
            self._add_interval(merged_start, merged_end)
            merged_start, merged_end = start, end
        self._add_interval(merged_start, merged_end)

    def add_read_support(self, mapped: bool = False, passing: bool = False,
                         best: bool = False, unique_best: bool = False,
                         ambiguous_best: bool = False,
                         high_confidence: bool = False,
                         score: float = 0.0):
        """Record one read-level support event without retaining its name."""
        if mapped:
            self.mapped_read_support += 1
        if passing:
            self.passing_read_support += 1
        if best:
            self.best_read_support += 1
        if unique_best:
            self.unique_best_read_support += 1
        if ambiguous_best:
            self.ambiguous_best_read_support += 1
        if high_confidence:
            self.high_confidence_read_support += 1
        try:
            self.support_score_sum += float(score)
        except (TypeError, ValueError):
            pass

    def _iter_depth_segments(self) -> Iterator[Tuple[int, int, int]]:
        """Yield 0-based half-open spans with constant depth."""
        if self.gene_length <= 0:
            return
        if self.compact:
            current_depth = 0
            previous = 0
            for raw_position, delta in sorted(self._depth_deltas.items()):
                position = max(0, min(raw_position, self.gene_length))
                if position > previous:
                    yield previous, position, current_depth
                current_depth += delta
                previous = position
            if previous < self.gene_length:
                yield previous, self.gene_length, current_depth
            return

        previous_depth = None
        segment_start = 0
        for position in range(self.gene_length):
            depth = self.position_depths.get(position, 0)
            if previous_depth is None:
                previous_depth = depth
                segment_start = position
            elif depth != previous_depth:
                yield segment_start, position, previous_depth
                segment_start = position
                previous_depth = depth
        if previous_depth is not None:
            yield segment_start, self.gene_length, previous_depth

    def coverage_at_depth(self, depth: int) -> float:
        try:
            threshold = max(1, int(depth))
        except (TypeError, ValueError):
            threshold = 1
        if threshold == 1:
            return self.gene_coverage
        return self._coverage_by_depth.get(threshold, 0.0)

    def finalise(self):
        # Calculate final statistics.
        if self.num_sequences > 0:
            # Older GeneStats instances may have been populated through the
            # identities list before identity_sum existed. Prefer the bounded
            # running total, but retain a compatibility fallback.
            total_identity = self.identity_sum if self.identity_sum else sum(self.identities)
            self.avg_identity = total_identity / self.num_sequences

        if self.gene_length <= 0:
            return

        thresholds = (1, 2, 3, 5, 10)
        covered_by_depth = {threshold: 0 for threshold in thresholds}
        total_depth = 0
        total_depth_squared = 0
        depth_lengths: Dict[int, int] = {}
        gap_spans: List[Tuple[int, int]] = []
        current_covered_run = 0

        for start, end, depth in self._iter_depth_segments():
            span = end - start
            if span <= 0:
                continue
            total_depth += depth * span
            total_depth_squared += depth * depth * span
            depth_lengths[depth] = depth_lengths.get(depth, 0) + span
            for threshold in thresholds:
                if depth >= threshold:
                    covered_by_depth[threshold] += span
            if depth > 0:
                current_covered_run += span
                self.longest_covered_run = max(
                    self.longest_covered_run,
                    current_covered_run,
                )
            else:
                gap_spans.append((start, end))
                current_covered_run = 0

        covered = covered_by_depth[1]
        self._coverage_by_depth = {
            threshold: length / self.gene_length * 100
            for threshold, length in covered_by_depth.items()
        }
        self.gene_coverage = self._coverage_by_depth[1]
        self.coverage_2x = self._coverage_by_depth[2]
        self.coverage_3x = self._coverage_by_depth[3]
        self.coverage_5x = self._coverage_by_depth[5]
        self.coverage_10x = self._coverage_by_depth[10]
        self.base_depth = total_depth / self.gene_length
        self.base_depth_hit = total_depth / covered if covered else 0.0

        mean_square = total_depth_squared / self.gene_length
        variance = max(0.0, mean_square - (self.base_depth * self.base_depth))
        self.depth_cv = (
            math.sqrt(variance) / self.base_depth
            if self.base_depth > 0
            else 0.0
        )

        midpoint = (self.gene_length + 1) / 2
        cumulative = 0
        for depth in sorted(depth_lengths):
            cumulative += depth_lengths[depth]
            if cumulative >= midpoint:
                self.median_depth = float(depth)
                break

        self.num_gaps = len(gap_spans)
        self.longest_gap = max(
            (end - start for start, end in gap_spans),
            default=0,
        )
        internal_gaps = [
            (start, end)
            for start, end in gap_spans
            if start > 0 and end < self.gene_length
        ]
        self.num_internal_gaps = len(internal_gaps)
        self.longest_internal_gap = max(
            (end - start for start, end in internal_gaps),
            default=0,
        )
