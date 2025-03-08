from unittest.mock import MagicMock, patch

import pytest
from quicat.cli.extract import ReadsProcessorFASTQ


@pytest.fixture
def mock_fastq():
    mock_fastq = MagicMock()
    mock_fastq.__iter__.return_value = iter(
        [("read1", "ACGT", "IIII"), ("read2", "ACGT", "####"), ("read3", "TGCA", "IIII")]
    )
    return mock_fastq


@pytest.fixture
def mock_pattern_matching():
    mock_pm = MagicMock()
    mock_pm.find_sequence.side_effect = lambda seq: seq if seq == "ACGT" else None
    return mock_pm


@patch("pyfastx.Fastx")
def test_reads_processor_single_end_with_filter_and_match(mock_pyfastx, mock_fastq, mock_pattern_matching):
    mock_pyfastx.return_value = mock_fastq

    processor = ReadsProcessorFASTQ(
        paired_end=False,
        read_qc_threshold=20,
        read_qc_percentage=75,
        phred33=True,
        fastq_path_r1="mock_path",
        fastq_path_r2=None,
        reference="ACGT",
        flanked_pattern=False,
        min_read_length=4,
        max_read_length=0,
        read_length=0,
        sequencing_technology="DNA",
    )
    processor.fastq1 = mock_fastq
    processor.pattern_matching = mock_pattern_matching

    matched_reads = processor._process_dna()
    assert matched_reads == ["ACGT"]
