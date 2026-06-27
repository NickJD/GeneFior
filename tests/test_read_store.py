from GeneFior.read_store import ReadMatchStore


def test_read_store_deduplicates_and_upgrades_passing(tmp_path):
    store = ReadMatchStore(tmp_path / "matches.sqlite", batch_size=2)
    try:
        store.add("db", "BLASTx", "geneA", "read1", passing=False)
        store.add("db", "BLASTx", "geneA", "read1", passing=True)
        store.add("db", "BLASTx", "geneA", "read2", passing=False)
        store.flush()

        output = tmp_path / "names"
        assert store.export_names(output, "all") == 2
        assert store.export_names(output, "passing") == 1
        names = (output / "db_BLASTx_geneA_read_names.txt").read_text().splitlines()
        assert names == ["read1"]
    finally:
        store.close()


def test_read_store_populates_sequences_and_exports_fasta(tmp_path):
    store = ReadMatchStore(tmp_path / "matches.sqlite")
    try:
        store.add("db", "BLASTn", "geneA", "read1", passing=True)
        store.add("db", "BLASTn", "geneA", "read2", passing=True)
        query = tmp_path / "reads.fasta"
        query.write_text(">read1 description\nACGT\n>read2\nTTAA\n")

        matched, updated = store.populate_sequences(query)
        assert matched == 2
        assert updated == 2

        output = tmp_path / "fasta"
        exported, missing = store.export_fasta(output, "passing")
        assert exported == 2
        assert missing == 0
        content = (output / "db_BLASTn_geneA_reads.fasta").read_text()
        assert ">read1\nACGT\n" in content
        assert ">read2\nTTAA\n" in content
    finally:
        store.close()


def test_detected_modes_filter_by_gene(tmp_path):
    store = ReadMatchStore(tmp_path / "matches.sqlite")
    try:
        store.add("db", "BWA", "geneA", "read1", passing=True, sequence="AAAA")
        store.add("db", "BWA", "geneA", "read2", passing=False, sequence="CCCC")
        store.add("db", "BWA", "geneB", "read3", passing=True, sequence="GGGG")

        output = tmp_path / "detected"
        assert store.export_names(output, "detected", {"geneA"}) == 1
        assert store.export_names(output, "detected-all", {"geneA"}) == 2
    finally:
        store.close()
