import gzip
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from GeneFior import utils


def _write_fastq(path, read_id):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        handle.write(f"@{read_id}\nACGT\n+\nIIII\n")


def _options(tmp_path, r1, r2, **overrides):
    values = {
        'sequence_type': 'Paired-FASTQ',
        'input': f"{r1},{r2}",
        'output': str(tmp_path / 'output'),
        'temp_directory': None,
        'tools': ['bowtie2'],
        'threads': 8,
        'fastp_trim': False,
        'force_modify_fastq': False,
        'no_cleanup': False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def _fake_fastp_run(commands):
    def run(command, **kwargs):
        commands.append((command, kwargs))
        out1 = Path(command[command.index('--out1') + 1])
        out2 = Path(command[command.index('--out2') + 1])
        json_report = Path(command[command.index('--json') + 1])
        html_report = Path(command[command.index('--html') + 1])
        _write_fastq(out1, 'trimmed/1')
        _write_fastq(out2, 'trimmed/2')
        json_report.write_text('{"summary": {}}\n')
        html_report.write_text('<html></html>\n')
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    return run


def test_fastp_helper_builds_paired_command_and_reports(tmp_path, monkeypatch):
    r1 = tmp_path / 'reads_R1.fastq.gz'
    r2 = tmp_path / 'reads_R2.fastq.gz'
    _write_fastq(r1, 'raw/1')
    _write_fastq(r2, 'raw/2')
    options = _options(tmp_path, r1, r2, fastp_trim=True, threads=32)
    commands = []

    monkeypatch.setattr(utils.shutil, 'which', lambda executable: '/usr/bin/fastp')
    monkeypatch.setattr(utils.subprocess, 'run', _fake_fastp_run(commands))

    trimmed_r1, trimmed_r2 = utils.run_fastp_for_paired_reads(
        options, str(r1), str(r2), logging.getLogger('test-fastp-helper')
    )

    command, kwargs = commands[0]
    assert command[0] == '/usr/bin/fastp'
    assert command[command.index('--in1') + 1] == str(r1)
    assert command[command.index('--in2') + 1] == str(r2)
    assert command[command.index('--thread') + 1] == '16'
    assert '--detect_adapter_for_pe' in command
    assert kwargs['check'] is False
    assert Path(trimmed_r1).is_file()
    assert Path(trimmed_r2).is_file()
    assert Path(options.fastp_report_json).is_file()
    assert Path(options.fastp_report_html).is_file()


def test_fastp_requested_without_executable_fails_clearly(tmp_path, monkeypatch):
    r1 = tmp_path / 'reads_R1.fastq'
    r2 = tmp_path / 'reads_R2.fastq'
    _write_fastq(r1, 'raw/1')
    _write_fastq(r2, 'raw/2')
    options = _options(tmp_path, r1, r2, fastp_trim=True)
    monkeypatch.setattr(utils.shutil, 'which', lambda executable: None)

    with pytest.raises(RuntimeError, match='fastp.*not found'):
        utils.run_fastp_for_paired_reads(
            options, str(r1), str(r2), logging.getLogger('test-fastp-missing')
        )


def test_fastp_failure_does_not_fall_back_to_raw_reads(tmp_path, monkeypatch):
    r1 = tmp_path / 'reads_R1.fastq'
    r2 = tmp_path / 'reads_R2.fastq'
    _write_fastq(r1, 'raw/1')
    _write_fastq(r2, 'raw/2')
    options = _options(tmp_path, r1, r2, fastp_trim=True)
    monkeypatch.setattr(utils.shutil, 'which', lambda executable: '/usr/bin/fastp')
    monkeypatch.setattr(
        utils.subprocess,
        'run',
        lambda command, **kwargs: SimpleNamespace(
            returncode=1, stdout='', stderr='invalid input'
        ),
    )

    with pytest.raises(RuntimeError, match='fastp failed.*invalid input'):
        utils.handle_all_input_files(
            options, logging.getLogger('test-fastp-failure')
        )

    assert not hasattr(options, 'input_fastq')


def test_untrimmed_paired_reads_are_logged_and_used(tmp_path, caplog):
    r1 = tmp_path / 'reads_R1.fastq'
    r2 = tmp_path / 'reads_R2.fastq'
    _write_fastq(r1, 'raw/1')
    _write_fastq(r2, 'raw/2')
    options = _options(tmp_path, r1, r2)

    with caplog.at_level(logging.INFO):
        utils.handle_all_input_files(options, logging.getLogger('test-fastp-disabled'))

    assert options.input_fastq == (str(r1), str(r2))
    assert options.input_fastq_original == (str(r1), str(r2))
    assert 'non-trimmed reads were used' in caplog.text


def test_fastp_outputs_feed_alignment_and_blast_conversion(
    tmp_path, monkeypatch
):
    r1 = tmp_path / 'reads_R1.fastq.gz'
    r2 = tmp_path / 'reads_R2.fastq.gz'
    _write_fastq(r1, 'raw/1')
    _write_fastq(r2, 'raw/2')
    options = _options(
        tmp_path, r1, r2, fastp_trim=True, tools=['bowtie2', 'blastn']
    )
    commands = []
    conversion_inputs = []

    monkeypatch.setattr(utils.shutil, 'which', lambda executable: '/usr/bin/fastp')
    monkeypatch.setattr(utils.subprocess, 'run', _fake_fastp_run(commands))

    def fake_fastq_to_fasta(current_options, logger):
        conversion_inputs.append(current_options.input)
        current_options.input_fasta = str(tmp_path / 'combined.fasta.gz')

    monkeypatch.setattr(utils, 'FASTQ_to_FASTA', fake_fastq_to_fasta)

    utils.handle_all_input_files(options, logging.getLogger('test-fastp-flow'))

    trimmed_r1, trimmed_r2 = options.input_fastq_trimmed
    assert options.input_fastq == (trimmed_r1, trimmed_r2)
    assert options.input_fastq_original == (str(r1), str(r2))
    assert conversion_inputs == [f"{trimmed_r1},{trimmed_r2}"]
    assert options.input == f"{r1},{r2}"
    assert trimmed_r1 in options.temp_files_to_cleanup
    assert trimmed_r2 in options.temp_files_to_cleanup


def test_batch_samples_receive_isolated_temporary_directories(tmp_path):
    first = utils.sample_temp_directory(
        tmp_path / 'scratch', tmp_path / 'output',
        tmp_path / 'output' / 'sample_a', 'sample_a',
        user_specified=True,
    )
    second = utils.sample_temp_directory(
        tmp_path / 'scratch', tmp_path / 'output',
        tmp_path / 'output' / 'sample_b', 'sample_b',
        user_specified=True,
    )

    assert first != second
    assert Path(first).name == 'sample_a'
    assert Path(second).name == 'sample_b'


def test_default_batch_temp_is_the_sample_output(tmp_path):
    sample_output = tmp_path / 'output' / 'sample_a'

    selected = utils.sample_temp_directory(
        tmp_path / 'output', tmp_path / 'output',
        sample_output, 'sample_a',
        user_specified=False,
    )

    assert selected == str(sample_output.resolve())
