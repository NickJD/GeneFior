#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
# The test data directories (Genes-AA, Genes-DNA, Paired-FASTQ, Single-FASTA)
# live in the same directory as this script. Use ROOT_DIR directly.
TEST_DATA_DIR="$ROOT_DIR"
TMP_DB_ROOT="$ROOT_DIR/tmp_user_dbs/resfinder"
OUTPUT_ROOT="$ROOT_DIR/outputs"

mkdir -p "$OUTPUT_ROOT"
# Use a relative path to the GeneFior CLI from this script's directory
# script is in test_data/, GeneFior.py lives in ../src/GeneFior/
GENEFIOR_PY="$ROOT_DIR/../src/GeneFior/GeneFior.py"
PYTHON=${PYTHON:-python3}

# Mode flags
COMPARE_ONLY=0
if [ "${1:-}" = "--compare-only" ] || [ "${1:-}" = "-c" ]; then
    COMPARE_ONLY=1
fi

fail_count=0
pass_count=0

log_and_fail() {
    echo "$1"
    fail_count=$((fail_count+1))
}

log_and_pass() {
    echo "$1"
    pass_count=$((pass_count+1))
}

prepare_db() {
    # Build minimal DBs from provided fasta/gz files into $TMP_DB_ROOT
    mkdir -p "$TMP_DB_ROOT/blast_aa" "$TMP_DB_ROOT/blast_dna" "$TMP_DB_ROOT/diamond" "$TMP_DB_ROOT/bowtie2" "$TMP_DB_ROOT/bwa"

    # Genes-AA -> protein DB
    if [ -f "$TEST_DATA_DIR/Genes-AA/all_aa.gz" ]; then
        echo "Preparing protein DB from Genes-AA/all_aa.gz"
        gunzip -c "$TEST_DATA_DIR/Genes-AA/all_aa.gz" > "$TMP_DB_ROOT/blast_aa/all_aa_blastdb.fa"
        if command -v makeblastdb >/dev/null 2>&1; then
            makeblastdb -in "$TMP_DB_ROOT/blast_aa/all_aa_blastdb.fa" -dbtype prot -out "$TMP_DB_ROOT/blast_aa/all_aa_blastdb" >/dev/null 2>&1 || true
        fi
        if command -v diamond >/dev/null 2>&1; then
            diamond makedb --in "$TMP_DB_ROOT/blast_aa/all_aa_blastdb.fa" -d "$TMP_DB_ROOT/diamond/all_aa_diamonddb" >/dev/null 2>&1 || true
            cp "$TMP_DB_ROOT/diamond/all_aa_diamonddb.dmnd" "$TMP_DB_ROOT/diamond/all_aa_diamonddb.dmnd" 2>/dev/null || true
        fi
    else
        echo "Protein genes file not found: $TEST_DATA_DIR/Genes-AA/all_aa.gz"
    fi

    # Genes-DNA -> nucleotide DB
    if [ -f "$TEST_DATA_DIR/Genes-DNA/all.gz" ]; then
        echo "Preparing nucleotide DB from Genes-DNA/all.gz"
        gunzip -c "$TEST_DATA_DIR/Genes-DNA/all.gz" > "$TMP_DB_ROOT/blast_dna/all_blastdb.fa"
        if command -v makeblastdb >/dev/null 2>&1; then
            makeblastdb -in "$TMP_DB_ROOT/blast_dna/all_blastdb.fa" -dbtype nucl -out "$TMP_DB_ROOT/blast_dna/all_blastdb" >/dev/null 2>&1 || true
        fi
        # Build indices for bowtie2/bwa if available
        if command -v bowtie2-build >/dev/null 2>&1; then
            bowtie2-build "$TMP_DB_ROOT/blast_dna/all_blastdb.fa" "$TMP_DB_ROOT/bowtie2/all_bowtie2db" >/dev/null 2>&1 || true
        fi
        if command -v bwa >/dev/null 2>&1; then
            bwa index -p "$TMP_DB_ROOT/bwa/all_bwadb" "$TMP_DB_ROOT/blast_dna/all_blastdb.fa" >/dev/null 2>&1 || true
        fi
    else
        echo "Nucleotide genes file not found: $TEST_DATA_DIR/Genes-DNA/all.gz"
    fi
}

run_case() {
    local name="$1"; shift
    local cmd=("$@")
    local outdir="$OUTPUT_ROOT/out_${name}"
    # If a previous output directory exists, move it aside so we can compare
    if [ -d "$outdir" ]; then
        if [ -d "${outdir}.prev" ]; then
            rm -rf "${outdir}.prev"
        fi
        mv "$outdir" "${outdir}.prev"
    fi
    mkdir -p "$outdir"

    echo "\n=== Running test: $name ==="
    echo "Command: ${cmd[*]}"
    set +e
    "${cmd[@]}" > "$outdir/run.log" 2>&1
    RC=$?
    set -e
    if [ $RC -ne 0 ]; then
        echo "Command failed with exit code $RC - see $outdir/run.log"
        echo "--- last 200 lines of $outdir/run.log ---"
        tail -n 200 "$outdir/run.log" || true
        echo "--- end log ---"
        return 1
    fi

    echo "Command succeeded"

    # If there is a previous run directory, compute checksums for both and compare
    if [ -d "${outdir}.prev" ]; then
        compare_with_prev "$outdir"
        return $?
    fi

    return 0
}

# Compare an existing output directory to its .prev directory using sha256 checksums
compare_with_prev() {
    local outdir="$1"
    local name
    name=$(basename "$outdir")
    # Choose checksum command
    if command -v sha256sum >/dev/null 2>&1; then
        SHACMD="sha256sum"
    elif command -v shasum >/dev/null 2>&1; then
        SHACMD="shasum -a 256"
    else
        echo "No sha256sum/shasum available; skipping checksum comparison for $name"
        return 0
    fi

    if [ ! -d "${outdir}.prev" ]; then
        echo "No previous run directory to compare for: $outdir"
        return 1
    fi

    TMP_NEW="$(mktemp -t checksums_new.XXXXXX)"
    TMP_OLD="$(mktemp -t checksums_old.XXXXXX)"
    # Exclude transient logs and generated checksum files from comparison
    (cd "$outdir" && find . -type f ! -name 'pipeline_*' ! -name 'run.log' ! -name 'checksum_diff.txt' -print0 | xargs -0 -I{} bash -c "$SHACMD '{}'" 2>/dev/null | sort) > "$TMP_NEW" 2>/dev/null || true
    (cd "${outdir}.prev" && find . -type f ! -name 'pipeline_*' ! -name 'run.log' ! -name 'checksum_diff.txt' -print0 | xargs -0 -I{} bash -c "$SHACMD '{}'" 2>/dev/null | sort) > "$TMP_OLD" 2>/dev/null || true

    if [ -s "$TMP_NEW" ] && [ -s "$TMP_OLD" ]; then
        if diff -u "$TMP_OLD" "$TMP_NEW" > "${outdir}/checksum_diff.txt"; then
            echo "Outputs identical to previous run for test: $(basename "$outdir")"
            rm -f "$TMP_NEW" "$TMP_OLD"
            return 0
        else
            echo "Outputs differ from previous run for test: $(basename "$outdir")"
            echo "Saved unified diff to: ${outdir}/checksum_diff.txt"
            echo "--- checksum diff (first 200 lines) ---"
            head -n 200 "${outdir}/checksum_diff.txt" || true
            echo "--- end diff ---"
            rm -f "$TMP_NEW" "$TMP_OLD"
            return 2
        fi
    else
        echo "Could not compute checksums for both runs (one side empty) for test: $(basename "$outdir")"
        rm -f "$TMP_NEW" "$TMP_OLD"
        return 2
    fi
}

#################### Prepare DBs from provided test inputs ########################
echo "Preparing temporary user DBs from files in $TEST_DATA_DIR"
prepare_db

if [ "$COMPARE_ONLY" -eq 1 ]; then
    echo "Running in compare-only mode: comparing existing outputs against their .prev counterparts"
    for t in genes_aa genes_dna paired single_fasta; do
        outdir="$OUTPUT_ROOT/out_${t}"
        if [ -d "$outdir" ]; then
            compare_with_prev "$outdir"
            RC=$?
            if [ $RC -eq 0 ]; then
                log_and_pass "Compare-only: $t identical to previous"
            elif [ $RC -eq 1 ]; then
                echo "Compare-only: $t skipped (no .prev)"
            else
                log_and_fail "Compare-only: $t differs from previous"
            fi
        else
            echo "Compare-only: output dir missing for $t: $outdir"
        fi
    done

    echo "\n=== Compare-only summary ==="
    echo "Passed: $pass_count"
    echo "Failed: $fail_count"
    if [ $fail_count -ne 0 ]; then
        exit 2
    else
        exit 0
    fi
fi

#################### Test A: Genes-FASTA (aa) ########################
GENES_AA_F="${TMP_DB_ROOT}/blast_aa/all_aa_blastdb.fa"
if [ ! -f "$GENES_AA_F" ]; then
    GENES_AA_F="$TEST_DATA_DIR/Genes-AA/all_aa.gz"
fi
run_case "genes_aa" $PYTHON "$GENEFIOR_PY" -i "$GENES_AA_F" -st Genes-FASTA --genes-type aa --db-path "$ROOT_DIR/tmp_user_dbs/resfinder" -o "$OUTPUT_ROOT/out_genes_aa" --tools all || log_and_fail "Genes-AA test failed"
if [ -f "$OUTPUT_ROOT/out_genes_aa/run_parameters.json" ]; then
    if grep -q 'blastp' "$OUTPUT_ROOT/out_genes_aa/run_parameters.json" 2>/dev/null || grep -q 'diamond' "$OUTPUT_ROOT/out_genes_aa/run_parameters.json" 2>/dev/null; then
        log_and_pass "Genes-AA tools expansion OK"
    else
        log_and_fail "Genes-AA tools expansion did not include blastp/diamond"
    fi
else
    log_and_fail "Genes-AA run_parameters.json missing"
fi

# Ensure BLASTp perc_identity error is absent
if grep -q "Unknown argument: \"perc_identity\"" "$OUTPUT_ROOT/out_genes_aa/run.log" 2>/dev/null; then
    log_and_fail "Detected blastp 'perc_identity' unknown-argument error in genes_aa run log"
fi

#################### Test B: Genes-FASTA (dna) ########################
GENES_DNA_F="${TMP_DB_ROOT}/blast_dna/all_blastdb.fa"
if [ ! -f "$GENES_DNA_F" ]; then
    GENES_DNA_F="$TEST_DATA_DIR/Genes-DNA/all.gz"
fi
run_case "genes_dna" $PYTHON "$GENEFIOR_PY" -i "$GENES_DNA_F" -st Genes-FASTA --genes-type dna --db-path "$ROOT_DIR/tmp_user_dbs/resfinder" -o "$OUTPUT_ROOT/out_genes_dna" --tools all || log_and_fail "Genes-DNA test failed"
if [ -f "$OUTPUT_ROOT/out_genes_dna/run_parameters.json" ]; then
    if grep -q 'blastn' "$OUTPUT_ROOT/out_genes_dna/run_parameters.json" 2>/dev/null || grep -q 'blastx' "$OUTPUT_ROOT/out_genes_dna/run_parameters.json" 2>/dev/null || grep -q 'diamond' "$OUTPUT_ROOT/out_genes_dna/run_parameters.json" 2>/dev/null; then
        log_and_pass "Genes-DNA tools expansion OK"
    else
        log_and_fail "Genes-DNA tools expansion did not include blastn/blastx/diamond"
    fi
else
    log_and_fail "Genes-DNA run_parameters.json missing"
fi

#################### Test C: Paired-FASTQ ########################
PAIRED_R1="$TEST_DATA_DIR/Paired-FASTQ/Test_R1.fastq.gz"
PAIRED_R2="$TEST_DATA_DIR/Paired-FASTQ/Test_R2.fastq.gz"
if [ -f "$PAIRED_R1" ] && [ -f "$PAIRED_R2" ]; then
    run_case "paired" $PYTHON "$GENEFIOR_PY" -i "$PAIRED_R1,$PAIRED_R2" -st Paired-FASTQ --db-path "$ROOT_DIR/tmp_user_dbs/resfinder" -o "$OUTPUT_ROOT/out_paired" --tools bowtie2 || log_and_fail "Paired-FASTQ test failed"
    # If bowtie2 and samtools present check for BAM file
    if command -v bowtie2 >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then
        BAMFILE=$(find "$OUTPUT_ROOT/out_paired" -type f -name '*_bowtie2_results_sorted.bam' | head -n1 || true)
        if [ -n "$BAMFILE" ]; then
            log_and_pass "Paired-FASTQ produced BAM output"
        else
            log_and_fail "Paired-FASTQ expected BAM output not found"
        fi
    else
        log_and_pass "Paired-FASTQ CLI ran; bowtie2/samtools not available to verify BAM"
    fi
else
    echo "Paired-FASTQ test files not found; skipping Paired-FASTQ test"
fi

#################### Test D: Single-FASTA ########################
SINGLE_FASTA_GZ=$(ls "$TEST_DATA_DIR/Single-FASTA"/*.gz 2>/dev/null | head -n1 || true)
if [ -n "$SINGLE_FASTA_GZ" ]; then
    run_case "single_fasta" $PYTHON "$GENEFIOR_PY" -i "$SINGLE_FASTA_GZ" -st Single-FASTA --db-path "$ROOT_DIR/tmp_user_dbs/resfinder" -o "$OUTPUT_ROOT/out_single_fasta" --tools all || log_and_fail "Single-FASTA test failed"
    if [ -f "$OUTPUT_ROOT/out_single_fasta/run_parameters.json" ]; then
        log_and_pass "Single-FASTA run completed and produced run_parameters.json"
    else
        log_and_fail "Single-FASTA run_parameters.json missing"
    fi
else
    echo "Single-FASTA test file not found; skipping Single-FASTA test"
fi

#################### Summary ########################
echo "\n=== Test summary ==="
echo "Passed: $pass_count"
echo "Failed: $fail_count"

# Cleanup temporary DBs
rm -rf "$ROOT_DIR/tmp_user_dbs"

if [ $fail_count -ne 0 ]; then
    echo "Some tests failed. See logs in $OUTPUT_ROOT/out_*"
    exit 2
else
    echo "All tests passed"
    exit 0
fi



