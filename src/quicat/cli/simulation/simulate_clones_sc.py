import random
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from .simulator_utils import (
    ALPHABET,
    generate_pcr_replicates,
    generate_random_barcodes_with_threshold,
    generate_unique_sequences,
)


class GenerateBarcodeSequencesSC:
    """
    A class to simulate single-cell sequencing data with various barcode distributions and PCR errors.

    This class provides methods to generate cell barcodes, UMI sequences, and clone-specific barcodes for
    single-cell RNA sequencing (scRNA-seq) data, simulate PCR errors, and generate paired-end reads
    formatted for FASTQ files.

    Attributes:
        alphabet (str): The alphabet used for generating random nucleotides.
    """

    def __init__(self) -> None:
        """Initializes the GenerateBarcodeSequencesSC class."""
        self.alphabet: str = ALPHABET

    def _generate_clone_distribution_uniform(self, num_clones: int, num_cells: int) -> List[int]:
        """
        Generate a uniform distribution of clones across cells.

        Parameters:
            num_clones (int): Number of clones to distribute.
            num_cells (int): Total number of cells.

        Returns:
            List[int]: A list of cell counts for each clone, evenly distributed.
        """
        return [num_cells // num_clones] * num_clones

    def _generate_clone_distribution_random(self, num_clones: int, num_cells: int) -> List[int]:
        """
        Generate a random distribution of clones across cells.

        Parameters:
            num_clones (int): Number of clones to distribute.
            num_cells (int): Total number of cells.

        Returns:
            List[int]: A list of random cell counts for each clone.
        """
        counts = np.array([random.randint(1, num_cells // num_clones * 2) for _ in range(num_clones)])
        counts = np.clip(counts, 1, None).astype(int)
        counts_casted: List[int] = list((counts / counts.sum() * num_cells).astype(int))
        return counts_casted

    def _generate_clone_distribution_normal(
        self, num_clones: int, num_cells: int, mean: float, std_dev: float
    ) -> List[int]:
        """
        Generate a normal distribution of clones across cells.

        Parameters:
            num_clones (int): Number of clones to distribute.
            num_cells (int): Total number of cells.
            mean (float): Mean value for the normal distribution.
            std_dev (float): Standard deviation for the normal distribution.

        Returns:
            List[int]: A list of cell counts for each clone following a normal distribution.
        """
        counts = np.random.normal(loc=mean, scale=std_dev, size=num_clones).astype(int)
        counts = np.clip(counts, 1, None)
        counts_casted: List[int] = list((counts / counts.sum() * num_cells).astype(int))
        return counts_casted

    def _generate_clone_distribution_powerlaw(self, num_clones: int, num_cells: int, exponent: float) -> List[int]:
        """
        Generate a power-law distribution of clones across cells.

        Parameters:
            num_clones (int): Number of clones to distribute.
            num_cells (int): Total number of cells.
            exponent (float): Exponent for the power-law distribution.

        Returns:
            List[int]: A list of cell counts for each clone following a power-law distribution.
        """
        counts = np.random.zipf(a=exponent, size=num_clones).astype(int)
        counts = np.clip(counts, 1, None)
        counts_casted: List[int] = list((counts / counts.sum() * num_cells).astype(int))
        return counts_casted

    def simulate_10x_data(
        self,
        num_cells: int,
        num_clones: int,
        cell_bc_length: int,
        umi_length: int,
        barcode_length: int,
        umis_per_cell: int,
        distribution_params: Dict[str, Any],
    ) -> DataFrame:
        """
        Simulate 10x Genomics single-cell RNA sequencing data with the specified distribution.

        Parameters:
            num_cells (int): Total number of cells to simulate.
            num_clones (int): Number of clones to distribute across cells.
            cell_bc_length (int): Length of the cell barcodes.
            umi_length (int): Length of the UMI sequences.
            barcode_length (int): Length of the clone-specific barcodes.
            umis_per_cell (int): Number of UMIs per cell.
            distribution_params (Dict[str, Any]): Additional parameters for the distribution.

        Returns:
            DataFrame: A DataFrame containing simulated cell barcodes, UMIs, and clone-specific barcodes.
        """
        distribution_type = distribution_params["distribution_type"]
        if distribution_type == "uniform":
            clone_distribution = self._generate_clone_distribution_uniform(num_clones, num_cells)
        elif distribution_type == "random":
            clone_distribution = self._generate_clone_distribution_random(num_clones, num_cells)
        elif distribution_type == "normal":
            mean = distribution_params.get("mean", 50.0)
            std_dev = distribution_params.get("std_dev", 10.0)
            clone_distribution = self._generate_clone_distribution_normal(num_clones, num_cells, mean, std_dev)
        elif distribution_type == "powerlaw":
            exponent = distribution_params.get("exponent", 1.5)
            clone_distribution = self._generate_clone_distribution_powerlaw(num_clones, num_cells, exponent)
        else:
            e = f"Unsupported distribution type: {distribution_type}"
            raise ValueError(e)

        cell_barcodes = generate_unique_sequences(num_cells, cell_bc_length)
        unique_umis = generate_unique_sequences(num_cells * umis_per_cell, umi_length)
        clone_barcodes = generate_random_barcodes_with_threshold(num_clones, barcode_length, 3)

        data = []
        umi_index = 0
        used_cell_barcodes = set()

        for clone_idx, num_clone_cells in enumerate(clone_distribution):
            clone_bc = clone_barcodes[clone_idx]
            available_cell_barcodes = [cb for cb in cell_barcodes if cb not in used_cell_barcodes]

            if len(available_cell_barcodes) < num_clone_cells:
                e = "Not enough unique cell barcodes to assign to all clones."
                raise ValueError(e)

            selected_cell_barcodes = random.sample(available_cell_barcodes, num_clone_cells)
            for cell_bc in selected_cell_barcodes:
                used_cell_barcodes.add(cell_bc)
                for _ in range(umis_per_cell):
                    umi = unique_umis[umi_index]
                    data.append((cell_bc, umi, clone_bc, f"clone_{clone_idx + 1}"))
                    umi_index += 1

        return pd.DataFrame(data, columns=["cell_bc", "umi", "barcode", "clone"])

    def introduce_pcr_errors(
        self,
        df: DataFrame,
        pcr_error_rate: float,
        max_pcr_hamming_dist: int,
        num_pcr_chimeras: int,
        pcr_error_fraction: float,
    ) -> DataFrame:
        """
        Introduce PCR errors into the simulated single-cell RNA sequencing data.

        Parameters:
            df (DataFrame): DataFrame containing the original simulated data.
            pcr_error_rate (float): Error rate for the PCR process.
            max_pcr_hamming_dist (int): Maximum Hamming distance for error-prone PCR replicates.
            num_pcr_chimeras (int): Number of PCR chimeras to generate.
            pcr_error_fraction (float): Fraction of UMIs to introduce PCR errors.

        Returns:
            DataFrame: A DataFrame containing original and error-introduced barcodes with their respective clones.
        """
        perturbed_data = []
        if pcr_error_fraction == 0:
            return df
        for clone in df["clone"].unique():
            clone_df = df[df["clone"] == clone].copy()
            major_barcode = clone_df["barcode"].iloc[0]

            error_barcodes = generate_pcr_replicates(
                major_barcode, num_pcr_chimeras, max_pcr_hamming_dist, pcr_error_rate
            )
            sampled_cells = clone_df["cell_bc"].unique()
            for cell in sampled_cells:
                cell_umis = clone_df[clone_df["cell_bc"] == cell].sample(frac=pcr_error_fraction, replace=False)
                for _, row in cell_umis.iterrows():
                    perturbed_data.append((cell, row["umi"], random.sample(error_barcodes, 1)[0], clone))
                remaining_umis = clone_df[clone_df["cell_bc"] == cell].drop(cell_umis.index)
                perturbed_data.extend(remaining_umis.values.tolist())
        return pd.DataFrame(perturbed_data, columns=["cell_bc", "umi", "barcode", "clone"])

    def extend_barcodes(self, df: DataFrame, sequence_length: int, flank_left: str, flank_right: str) -> DataFrame:
        """
        Extend barcodes with left and right flanks and random nucleotides to reach the specified sequence length.

        Parameters:
            df (DataFrame): DataFrame containing the barcodes to extend.
            sequence_length (int): The desired total length of the sequence.
            flank_left (str): Left flank to add to each barcode.
            flank_right (str): Right flank to add to each barcode.

        Returns:
            DataFrame: A DataFrame with paired-end reads ("read_1" and "read_2").
        """
        reads1 = []
        reads2 = []

        for _, row in df.iterrows():
            CB = row["cell_bc"]
            UMI = row["umi"]
            reads1.append(f"{CB}{UMI}{''.join(random.choices(self.alphabet, k=74))}")

            read2 = row["barcode"]
            full_sequence = f"{flank_left}{read2}{flank_right}"
            remaining_length = sequence_length - len(full_sequence)

            if remaining_length > 0:
                extension_each_side = remaining_length // 2
                left_extension = "".join(random.choices(self.alphabet, k=extension_each_side))
                right_extension = "".join(random.choices(self.alphabet, k=remaining_length - extension_each_side))
                read2 = f"{left_extension}{full_sequence}{right_extension}"

            reads2.append(read2)

        df["read_1"], df["read_2"] = reads1, reads2
        return df

    def generate_fastq_generator_input(
        self,
        num_cells: int,
        num_clones: int,
        cell_bc_length: int,
        umi_length: int,
        barcode_length: int,
        umis_per_cell: int,
        pcr_error_rate: float,
        max_pcr_hamming_dist: int,
        num_pcr_chimeras: int,
        pcr_error_fraction: float,
        sequence_length: int,
        flank_left: str,
        flank_right: str,
        distribution_params: Dict[str, Any],
    ) -> Tuple[DataFrame, DataFrame]:
        """
        Generate input data for a FASTQ generator based on simulated 10x Genomics scRNA-seq data.

        Parameters:
            num_cells (int): Total number of cells to simulate.
            num_clones (int): Number of clones to distribute across cells.
            cell_bc_length (int): Length of the cell barcodes.
            umi_length (int): Length of the UMI sequences.
            barcode_length (int): Length of the clone-specific barcodes.
            umis_per_cell (int): Number of UMIs per cell.
            pcr_error_rate (float): Error rate for the PCR process.
            max_pcr_hamming_dist (int): Maximum Hamming distance for error-prone PCR replicates.
            num_pcr_chimeras (int): Number of PCR chimeras to generate.
            pcr_error_fraction (float): Fraction of UMIs to introduce PCR errors.
            sequence_length (int): The total length of the final sequence after extending with flanks.
            flank_left (str): Left flank to add to each barcode.
            flank_right (str): Right flank to add to each barcode.
            distribution_params (Dict[str, Any]): Additional parameters for the distribution type.

        Returns:
            Tuple[DataFrame, DataFrame]: A tuple containing the ground truth DataFrame and the extended sequences DataFrame for FASTQ generation.
        """
        simulated_data = self.simulate_10x_data(
            num_cells=num_cells,
            num_clones=num_clones,
            cell_bc_length=cell_bc_length,
            umi_length=umi_length,
            barcode_length=barcode_length,
            umis_per_cell=umis_per_cell,
            distribution_params=distribution_params,
        )
        perturbed_data = self.introduce_pcr_errors(
            simulated_data,
            pcr_error_rate=pcr_error_rate,
            max_pcr_hamming_dist=max_pcr_hamming_dist,
            num_pcr_chimeras=num_pcr_chimeras,
            pcr_error_fraction=pcr_error_fraction,
        )

        extended_data = self.extend_barcodes(
            perturbed_data, sequence_length=sequence_length, flank_left=flank_left, flank_right=flank_right
        )
        ground_truth = simulated_data[["cell_bc", "barcode", "clone"]].drop_duplicates(subset="cell_bc")
        fastq_generator_input = extended_data[["read_1", "read_2"]]

        return ground_truth, fastq_generator_input
