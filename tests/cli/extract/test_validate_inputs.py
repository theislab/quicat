import os
from tempfile import NamedTemporaryFile

import pandas as pd
import pytest
from quicat.cli.utils import ValidateInputs


# Fixtures to create temporary CSV files for testing
@pytest.fixture
def tmp_csv_file():
    """Creates a temporary CSV file for testing."""
    with NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as tmp:
        yield tmp.name
    os.unlink(tmp.name)


@pytest.fixture
def mock_exists(monkeypatch):
    """Mocks os.path.exists to always return True."""
    monkeypatch.setattr(os.path, "exists", lambda path: True)


# Test cases for ValidateInputs
def test_validate_dna_single_end_valid(tmp_csv_file, mock_exists):
    # Configuration for DNA single-end
    config_dict = {"sequencing_technology": "DNA", "paired_end": False, "input": "fastq"}

    # Create a valid CSV
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "condition": [None, None],
            "replicate": [None, None],
            "fastq_path_r1": ["path/to/fastq1.fastq", "path/to/fastq2.fastq"],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    df_validated, use_condition, use_replicate = validator._validate()

    assert not df_validated.empty
    assert "sample" in df_validated.columns
    assert "fastq_path_r1" in df_validated.columns


def test_validate_dna_paired_end_valid(tmp_csv_file, mock_exists):
    # Configuration for DNA paired-end
    config_dict = {"sequencing_technology": "DNA", "paired_end": True, "input": "fastq"}

    # Create a valid CSV
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "condition": [None, None],
            "replicate": [None, None],
            "fastq_path_r1": ["path/to/fastq1_r1.fastq", "path/to/fastq2_r1.fastq"],
            "fastq_path_r2": ["path/to/fastq1_r2.fastq", "path/to/fastq2_r2.fastq"],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    df_validated, use_condition, use_replicate = validator._validate()

    assert not df_validated.empty
    assert "sample" in df_validated.columns
    assert "fastq_path_r1" in df_validated.columns
    assert "fastq_path_r2" in df_validated.columns


def test_validate_10x_fastq_valid(tmp_csv_file, mock_exists):
    # Configuration for 10x with fastq input
    config_dict = {"sequencing_technology": "10x", "input": "fastq"}

    # Create a valid CSV
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "condition": [None, None],
            "replicate": [None, None],
            "fastq_path_r1": ["path/to/fastq1_r1.fastq", "path/to/fastq2_r1.fastq"],
            "fastq_path_r2": ["path/to/fastq1_r2.fastq", "path/to/fastq2_r2.fastq"],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    df_validated, use_condition, use_replicate = validator._validate()

    assert not df_validated.empty
    assert "sample" in df_validated.columns
    assert "fastq_path_r1" in df_validated.columns
    assert "fastq_path_r2" in df_validated.columns


def test_validate_10x_bam_valid(tmp_csv_file, mock_exists):
    # Configuration for 10x with bam input
    config_dict = {"sequencing_technology": "10x", "input": "bam"}

    # Create a valid CSV
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "condition": [None, None],
            "replicate": [None, None],
            "bam_path": ["path/to/bam1.bam", "path/to/bam2.bam"],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    df_validated, use_condition, use_replicate = validator._validate()

    assert not df_validated.empty
    assert "sample" in df_validated.columns
    assert "bam_path" in df_validated.columns


def test_build_config_with_condition_replicate(tmp_csv_file, mock_exists):
    # Configuration for 10x with fastq input
    config_dict = {"sequencing_technology": "10x", "input": "fastq", "cells_barcodes_whitelist": False}

    # Create a valid CSV with 'condition' and 'replicate'
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "fastq_path_r1": ["path/to/fastq1_r1.fastq", "path/to/fastq2_r1.fastq"],
            "fastq_path_r2": ["path/to/fastq1_r2.fastq", "path/to/fastq2_r2.fastq"],
            "condition": ["treated", "control"],
            "replicate": [1, 2],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    final_config = validator.build_config()

    assert "sample1" in final_config
    assert final_config["sample1"]["condition"] == "treated"
    assert final_config["sample1"]["replicate"] == "1"  # Converted to string in _add_metadata


def test_build_config_missing_optional_columns(tmp_csv_file, mock_exists):
    # Configuration for 10x with bam input
    config_dict = {"sequencing_technology": "10x", "input": "bam", "cells_barcodes_whitelist": False}

    # Create a valid CSV without 'condition' and 'replicate'
    df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "condition": [None, None],
            "replicate": [None, None],
            "bam_path": ["path/to/bam1.bam", "path/to/bam2.bam"],
        }
    )
    df.to_csv(tmp_csv_file, index=False)

    validator = ValidateInputs(tmp_csv_file, config_dict)
    final_config = validator.build_config()

    assert "sample1" in final_config
