import gzip
import os
import re
import sqlite3
from pathlib import Path
from typing import Iterator, Optional, Set, Tuple


def _safe_filename(value: str) -> str:
    return re.sub(r"[^\w.-]", "_", value)


def _open_text(path: Path):
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, "rt")
    return open(path, "r")


def _iter_sequences(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (record_id, sequence) from FASTA or conventional four-line FASTQ."""
    with _open_text(path) as handle:
        first = ""
        while not first:
            first = handle.readline()
            if not first:
                return
            first = first.strip()

        if first.startswith(">"):
            record_id = first[1:].split()[0]
            parts = []
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    yield record_id, "".join(parts)
                    record_id = line[1:].split()[0]
                    parts = []
                elif line:
                    parts.append(line)
            yield record_id, "".join(parts)
            return

        if first.startswith("@"):
            header = first
            while header:
                sequence = handle.readline().strip()
                plus = handle.readline()
                quality = handle.readline()
                if not plus or not quality:
                    return
                yield header[1:].split()[0], sequence
                header = handle.readline().strip()
            return

        raise ValueError(f"Unsupported query sequence format: {path}")


class ReadMatchStore:
    """Disk-backed, deduplicated read-to-gene match storage.

    This is used only when read-name or FASTA output is explicitly requested.
    It keeps large recompute and gene-stats runs bounded in memory.
    """

    def __init__(self, path: Path, batch_size: int = 10000):
        self.path = Path(path)
        self.batch_size = max(1, int(batch_size))
        self._buffer = []
        self._connection = sqlite3.connect(str(self.path))
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute("PRAGMA cache_size=-32768")
        self._connection.execute(
            """
            CREATE TABLE IF NOT EXISTS matches (
                database_name TEXT NOT NULL,
                tool_name TEXT NOT NULL,
                gene_name TEXT NOT NULL,
                read_name TEXT NOT NULL,
                passing INTEGER NOT NULL DEFAULT 0,
                sequence TEXT,
                PRIMARY KEY (database_name, tool_name, gene_name, read_name)
            ) WITHOUT ROWID
            """
        )

    def add(self, database: str, tool: str, gene: str, read_name: str,
            passing: bool = False, sequence: Optional[str] = None):
        self._buffer.append((
            str(database), str(tool), str(gene), str(read_name),
            1 if passing else 0, sequence,
        ))
        if len(self._buffer) >= self.batch_size:
            self.flush()

    def flush(self):
        if not self._buffer:
            return
        self._connection.executemany(
            """
            INSERT INTO matches (
                database_name, tool_name, gene_name, read_name, passing, sequence
            ) VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT (database_name, tool_name, gene_name, read_name)
            DO UPDATE SET
                passing = CASE
                    WHEN excluded.passing > matches.passing
                    THEN excluded.passing ELSE matches.passing END,
                sequence = CASE
                    WHEN matches.sequence IS NULL OR matches.sequence = ''
                    THEN excluded.sequence ELSE matches.sequence END
            """,
            self._buffer,
        )
        self._connection.commit()
        self._buffer.clear()

    def populate_sequences(self, query_path: Path,
                           suffix: Optional[str] = None) -> Tuple[int, int]:
        """Fill missing sequences by streaming a FASTA/FASTQ file once."""
        self.flush()
        self._connection.execute(
            "CREATE INDEX IF NOT EXISTS idx_matches_read_name ON matches(read_name)"
        )
        matched_records = 0
        updated_rows = 0
        for read_name, sequence in _iter_sequences(Path(query_path)):
            candidate_names = [read_name]
            if suffix and not read_name.endswith(('/1', '/2')):
                candidate_names.insert(0, f"{read_name}{suffix}")
            row_updates = 0
            for candidate_name in candidate_names:
                cursor = self._connection.execute(
                    """
                    UPDATE matches
                    SET sequence = ?
                    WHERE read_name = ? AND (sequence IS NULL OR sequence = '')
                    """,
                    (sequence, candidate_name),
                )
                row_updates += cursor.rowcount
                if cursor.rowcount:
                    break
            if row_updates:
                matched_records += 1
                updated_rows += row_updates
            if matched_records and matched_records % self.batch_size == 0:
                self._connection.commit()
        self._connection.commit()
        return matched_records, updated_rows

    def _iter_selected(self, mode: str, detected_genes: Optional[Set[str]] = None):
        self.flush()
        detected = detected_genes or set()
        cursor = self._connection.execute(
            """
            SELECT database_name, tool_name, gene_name, read_name, passing, sequence
            FROM matches
            ORDER BY database_name, tool_name, gene_name, read_name
            """
        )
        for row in cursor:
            gene = row[2]
            passing = bool(row[4])
            if mode == "all":
                yield row
            elif mode == "passing":
                if passing:
                    yield row
            elif mode == "detected":
                if gene in detected and passing:
                    yield row
            elif mode == "detected-all":
                if gene in detected:
                    yield row
            else:
                raise ValueError(f"Unsupported read selection mode: {mode}")

    def export_names(self, output_dir: Path, mode: str,
                     detected_genes: Optional[Set[str]] = None) -> int:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        current_key = None
        handle = None
        exported = 0
        try:
            for database, tool, gene, read_name, _, _ in self._iter_selected(mode, detected_genes):
                key = (database, tool, gene)
                if key != current_key:
                    if handle:
                        handle.close()
                    filename = (
                        f"{_safe_filename(database)}_{_safe_filename(tool)}_"
                        f"{_safe_filename(gene)}_read_names.txt"
                    )
                    handle = open(output_dir / filename, "w")
                    current_key = key
                handle.write(f"{read_name}\n")
                exported += 1
        finally:
            if handle:
                handle.close()
        return exported

    def export_fasta(self, output_dir: Path, mode: str,
                     detected_genes: Optional[Set[str]] = None) -> Tuple[int, int]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        current_key = None
        handle = None
        exported = 0
        missing = 0
        try:
            for database, tool, gene, read_name, _, sequence in self._iter_selected(
                    mode, detected_genes):
                if not sequence:
                    missing += 1
                    continue
                key = (database, tool, gene)
                if key != current_key:
                    if handle:
                        handle.close()
                    filename = (
                        f"{_safe_filename(database)}_{_safe_filename(tool)}_"
                        f"{_safe_filename(gene)}_reads.fasta"
                    )
                    handle = open(output_dir / filename, "w")
                    current_key = key
                handle.write(f">{read_name}\n")
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\n")
                exported += 1
        finally:
            if handle:
                handle.close()
        return exported, missing

    def close(self, delete: bool = False):
        try:
            self.flush()
            self._connection.close()
        finally:
            if delete:
                try:
                    os.unlink(self.path)
                except FileNotFoundError:
                    pass
