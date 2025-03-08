import pandas as pd
import pytest
from quicat.cli.simulation import GenerateBarcodeSequencesSC


# Mocking the simulator_utils functions and constants
def mock_generate_unique_sequences(num_sequences, sequence_length):
    # Generate dummy unique sequences of specified length
    return [f"SEQ{str(i).zfill(sequence_length - 3)}" for i in range(num_sequences)]


def mock_generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist):
    # Generate dummy barcodes of the specified length
    return [f"BRC{str(i).zfill(barcode_length - 3)}" for i in range(num_barcodes)]


def mock_generate_pcr_replicates(barcode, num_pcr_chimeras, max_hamming_dist, pcr_error_rate):
    # Generate dummy PCR replicates by modifying the last character
    replicates = [
        barcode[:-1] + chr(((ord(barcode[-1]) - ord("A") + i) % 26) + ord("A")) for i in range(num_pcr_chimeras)
    ]
    return replicates


ALPHABET = "ACGT"


# Patching the dependencies in the GenerateBarcodeSequencesSC class
@pytest.fixture(autouse=True)
def patch_simulator_utils(monkeypatch):
    monkeypatch.setattr(
        "quicat.cli.simulation.simulator_utils.generate_unique_sequences", mock_generate_unique_sequences
    )
    monkeypatch.setattr(
        "quicat.cli.simulation.simulator_utils.generate_random_barcodes_with_threshold",
        mock_generate_random_barcodes_with_threshold,
    )
    monkeypatch.setattr("quicat.cli.simulation.simulator_utils.generate_pcr_replicates", mock_generate_pcr_replicates)
    monkeypatch.setattr("quicat.cli.simulation.simulator_utils.ALPHABET", ALPHABET)


# Test cases for GenerateBarcodeSequencesSC


def test_generate_clone_distribution_uniform():
    generator = GenerateBarcodeSequencesSC()
    num_clones = 5
    num_cells = 100
    distribution = generator._generate_clone_distribution_uniform(num_clones, num_cells)
    assert len(distribution) == num_clones
    assert all(count == num_cells // num_clones for count in distribution)


def test_generate_clone_distribution_random():
    generator = GenerateBarcodeSequencesSC()
    num_clones = 5
    num_cells = 100
    distribution = generator._generate_clone_distribution_random(num_clones, num_cells)
    assert len(distribution) == num_clones
    assert all(count >= 1 for count in distribution)


def test_generate_clone_distribution_normal():
    generator = GenerateBarcodeSequencesSC()
    num_clones = 5
    num_cells = 100
    mean = 20
    std_dev = 5
    distribution = generator._generate_clone_distribution_normal(num_clones, num_cells, mean, std_dev)
    assert len(distribution) == num_clones
    assert all(count >= 1 for count in distribution)


def test_generate_clone_distribution_powerlaw():
    generator = GenerateBarcodeSequencesSC()
    num_clones = 5
    num_cells = 100
    exponent = 1.5
    distribution = generator._generate_clone_distribution_powerlaw(num_clones, num_cells, exponent)
    assert len(distribution) == num_clones


def test_simulate_10x_data():
    generator = GenerateBarcodeSequencesSC()
    num_cells = 10
    num_clones = 2
    cell_bc_length = 8
    umi_length = 10
    barcode_length = 6
    umis_per_cell = 5
    distribution_params = {"distribution_type": "uniform"}

    df = generator.simulate_10x_data(
        num_cells=num_cells,
        num_clones=num_clones,
        cell_bc_length=cell_bc_length,
        umi_length=umi_length,
        barcode_length=barcode_length,
        umis_per_cell=umis_per_cell,
        distribution_params=distribution_params,
    )

    expected_rows = num_cells * umis_per_cell
    assert len(df) == expected_rows
    assert all(len(cb) == cell_bc_length for cb in df["cell_bc"])
    assert all(len(umi) == umi_length for umi in df["umi"])
    assert all(len(bc) == barcode_length for bc in df["barcode"])
    assert df["clone"].nunique() == num_clones


def test_introduce_pcr_errors():
    generator = GenerateBarcodeSequencesSC()
    # Create a simple DataFrame
    data = {
        "cell_bc": ["CELL001", "CELL001", "CELL002"],
        "umi": ["UMI0001", "UMI0002", "UMI0003"],
        "barcode": ["BRC001", "BRC001", "BRC001"],
        "clone": ["clone_1", "clone_1", "clone_1"],
    }
    df = pd.DataFrame(data)
    pcr_error_rate = 0.01
    max_pcr_hamming_dist = 1
    num_pcr_chimeras = 2
    pcr_error_fraction = 0.5

    perturbed_df = generator.introduce_pcr_errors(
        df,
        pcr_error_rate=pcr_error_rate,
        max_pcr_hamming_dist=max_pcr_hamming_dist,
        num_pcr_chimeras=num_pcr_chimeras,
        pcr_error_fraction=pcr_error_fraction,
    )

    assert len(perturbed_df) == len(df)
    assert "barcode" in perturbed_df.columns
    assert perturbed_df["barcode"].nunique() > 1


def test_extend_barcodes():
    generator = GenerateBarcodeSequencesSC()
    # Create a simple DataFrame
    data = {
        "cell_bc": ["CELL001", "CELL002"],
        "umi": ["UMI0001", "UMI0002"],
        "barcode": ["BRC001", "BRC002"],
        "clone": ["clone_1", "clone_2"],
    }
    df = pd.DataFrame(data)
    sequence_length = 100
    flank_left = "GGG"
    flank_right = "TTT"

    extended_df = generator.extend_barcodes(df, sequence_length, flank_left, flank_right)

    assert "read_1" in extended_df.columns
    assert "read_2" in extended_df.columns
    assert len(extended_df["read_2"][0]) == sequence_length


def test_generate_fastq_generator_input():
    generator = GenerateBarcodeSequencesSC()
    num_cells = 10
    num_clones = 2
    cell_bc_length = 8
    umi_length = 10
    barcode_length = 6
    umis_per_cell = 5
    pcr_error_rate = 0.01
    max_pcr_hamming_dist = 1
    num_pcr_chimeras = 2
    pcr_error_fraction = 0.5
    sequence_length = 100
    flank_left = "GGG"
    flank_right = "TTT"
    distribution_params = {"distribution_type": "uniform"}

    ground_truth, fastq_input = generator.generate_fastq_generator_input(
        num_cells=num_cells,
        num_clones=num_clones,
        cell_bc_length=cell_bc_length,
        umi_length=umi_length,
        barcode_length=barcode_length,
        umis_per_cell=umis_per_cell,
        pcr_error_rate=pcr_error_rate,
        max_pcr_hamming_dist=max_pcr_hamming_dist,
        num_pcr_chimeras=num_pcr_chimeras,
        pcr_error_fraction=pcr_error_fraction,
        sequence_length=sequence_length,
        flank_left=flank_left,
        flank_right=flank_right,
        distribution_params=distribution_params,
    )

    assert not ground_truth.empty
    assert not fastq_input.empty
    assert "read_1" in fastq_input.columns
    assert "read_2" in fastq_input.columns
    assert len(fastq_input) == num_cells * umis_per_cell
    assert len(fastq_input["read_2"][0]) == sequence_length


def test_generate_fastq_generator_input_invalid_distribution():
    generator = GenerateBarcodeSequencesSC()
    num_cells = 10
    num_clones = 2
    cell_bc_length = 8
    umi_length = 10
    barcode_length = 6
    umis_per_cell = 5
    pcr_error_rate = 0.01
    max_pcr_hamming_dist = 1
    num_pcr_chimeras = 2
    pcr_error_fraction = 0.5
    sequence_length = 100
    flank_left = "GGG"
    flank_right = "TTT"
    distribution_params = {"distribution_type": "invalid_type"}

    with pytest.raises(ValueError, match="Unsupported distribution type: invalid_type"):
        generator.generate_fastq_generator_input(
            num_cells=num_cells,
            num_clones=num_clones,
            cell_bc_length=cell_bc_length,
            umi_length=umi_length,
            barcode_length=barcode_length,
            umis_per_cell=umis_per_cell,
            pcr_error_rate=pcr_error_rate,
            max_pcr_hamming_dist=max_pcr_hamming_dist,
            num_pcr_chimeras=num_pcr_chimeras,
            pcr_error_fraction=pcr_error_fraction,
            sequence_length=sequence_length,
            flank_left=flank_left,
            flank_right=flank_right,
            distribution_params=distribution_params,
        )
