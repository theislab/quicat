import numpy as np
import pandas as pd
import pytest
from quicat.cli.simulation import GenerateBarcodeSequencesDNA


# Mocking the simulator_utils functions and constants
def mock_generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist):
    # Generate dummy barcodes of the specified length
    barcodes = ["A" * barcode_length for _ in range(num_barcodes)]
    return barcodes


def mock_generate_pcr_replicates(barcode, num_pcr_chimeras, max_hamming_dist, pcr_error_rate):
    # Generate dummy PCR replicates with slight modifications
    replicates = [barcode[:-1] + "G" for _ in range(num_pcr_chimeras)]
    return replicates


ALPHABET = "ACGT"


# Patching the dependencies in the GenerateBarcodeSequencesDNA class
@pytest.fixture(autouse=True)
def patch_simulator_utils(monkeypatch):
    monkeypatch.setattr(
        "quicat.cli.simulation.simulator_utils.generate_random_barcodes_with_threshold",
        mock_generate_random_barcodes_with_threshold,
    )
    monkeypatch.setattr("quicat.cli.simulation.simulator_utils.generate_pcr_replicates", mock_generate_pcr_replicates)
    monkeypatch.setattr("quicat.cli.simulation.simulator_utils.ALPHABET", ALPHABET)


# Test cases for GenerateBarcodeSequencesDNA
def test_generate_barcodes_uniform():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 5
    barcode_length = 10
    min_hamming_dist = 2
    counts = 100

    df = generator.generate_barcodes_uniform(num_barcodes, barcode_length, min_hamming_dist, counts)

    assert len(df) == num_barcodes
    assert all(df["barcode"].apply(lambda x: len(x) == barcode_length))
    assert all(df["count"] == counts)
    assert np.isclose(df["frequency"].sum(), 100.0)


def test_generate_barcodes_random():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 5
    barcode_length = 10
    min_hamming_dist = 2
    min_frequency = 0.01
    max_count = 1000

    df = generator.generate_barcodes_random(num_barcodes, barcode_length, min_hamming_dist, min_frequency, max_count)

    assert len(df) == num_barcodes
    assert all(df["barcode"].apply(lambda x: len(x) == barcode_length))
    assert all(df["count"] >= int(min_frequency * df["count"].sum() / 100))
    assert np.isclose(df["frequency"].sum(), 100.0)


def test_generate_barcodes_normal():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 5
    barcode_length = 8
    min_hamming_dist = 2
    mean = 50
    std_dev = 10
    min_frequency = 0.01

    df = generator.generate_barcodes_normal(
        num_barcodes, barcode_length, min_hamming_dist, mean, std_dev, min_frequency
    )

    assert len(df) == num_barcodes
    assert all(df["barcode"].apply(lambda x: len(x) == barcode_length))
    assert all(df["count"] >= 1)
    assert np.isclose(df["frequency"].sum(), 100.0)


def test_generate_barcodes_powerlaw():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 5
    barcode_length = 6
    min_hamming_dist = 2
    exponent = 1.5
    min_frequency = 0.01

    df = generator.generate_barcodes_powerlaw(num_barcodes, barcode_length, min_hamming_dist, exponent, min_frequency)

    assert len(df) == num_barcodes
    assert all(df["barcode"].apply(lambda x: len(x) == barcode_length))
    assert all(df["count"] >= int(min_frequency * df["count"].sum() / 100))
    assert np.isclose(df["frequency"].sum(), 100.0)


def test_simulate_pcr_with_errors():
    generator = GenerateBarcodeSequencesDNA()
    barcode_df = pd.DataFrame({"barcode": ["AAAAAA", "CCCCCC"], "count": [100, 200], "frequency": [33.33, 66.67]})
    num_pcr_chimeras = 2
    max_hamming_dist = 1
    pcr_error_rate = 0.01
    pcr_error_fraction = 0.5

    pcr_df = generator.simulate_pcr_with_errors(
        barcode_df, num_pcr_chimeras, max_hamming_dist, pcr_error_rate, pcr_error_fraction
    )

    expected_num_rows = len(barcode_df) + len(barcode_df) * num_pcr_chimeras
    assert len(pcr_df) == expected_num_rows
    assert "cluster" in pcr_df.columns
    assert "total_count" in pcr_df.columns


def test_generate_fastq_generator_input():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 3
    barcode_length = 6
    min_hamming_dist = 1
    num_pcr_chimeras = 2
    max_pcr_hamming_dist = 1
    pcr_error_rate = 0.01
    pcr_error_fraction = 0.5
    sequence_length = 20
    flank_left = "GGG"
    flank_right = "TTT"
    distribution_params = {"distribution_type": "uniform", "counts": 10}

    ground_truth, extended_reads = generator.generate_fastq_generator_input(
        num_barcodes=num_barcodes,
        barcode_length=barcode_length,
        min_hamming_dist=min_hamming_dist,
        num_pcr_chimeras=num_pcr_chimeras,
        max_pcr_hamming_dist=max_pcr_hamming_dist,
        pcr_error_rate=pcr_error_rate,
        pcr_error_fraction=pcr_error_fraction,
        sequence_length=sequence_length,
        flank_left=flank_left,
        flank_right=flank_right,
        distribution_params=distribution_params,
        art=False,
    )

    assert len(ground_truth) == num_barcodes
    assert len(extended_reads) > 0
    assert all(extended_reads["read"].apply(lambda x: len(x) == sequence_length))


def test_generate_fastq_generator_input_invalid_distribution():
    generator = GenerateBarcodeSequencesDNA()
    num_barcodes = 3
    barcode_length = 6
    min_hamming_dist = 1
    num_pcr_chimeras = 2
    max_pcr_hamming_dist = 1
    pcr_error_rate = 0.01
    pcr_error_fraction = 0.5
    sequence_length = 20
    flank_left = "GGG"
    flank_right = "TTT"
    distribution_params = {"distribution_type": "invalid_type"}

    with pytest.raises(ValueError, match="Unsupported distribution type: invalid_type"):
        generator.generate_fastq_generator_input(
            num_barcodes=num_barcodes,
            barcode_length=barcode_length,
            min_hamming_dist=min_hamming_dist,
            num_pcr_chimeras=num_pcr_chimeras,
            max_pcr_hamming_dist=max_pcr_hamming_dist,
            pcr_error_rate=pcr_error_rate,
            pcr_error_fraction=pcr_error_fraction,
            sequence_length=sequence_length,
            flank_left=flank_left,
            flank_right=flank_right,
            distribution_params=distribution_params,
            art=False,
        )
