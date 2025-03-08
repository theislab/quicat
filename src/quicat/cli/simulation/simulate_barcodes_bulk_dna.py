import random
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from .simulator_utils import (
    ALPHABET,
    generate_pcr_replicates,
    generate_random_barcodes_with_threshold,
)


class GenerateBarcodeSequencesDNA:
    """
    A class to generate and simulate DNA barcode sequences with various distributions.

    This class provides methods to generate barcodes using different statistical distributions,
    simulate the effects of PCR replication with errors, and extend barcodes with additional flanking sequences.

    Attributes:
        alphabet (str): The alphabet used for generating random nucleotides.
    """

    def __init__(self) -> None:
        """Initializes the GenerateBarcodeSequencesDNA class."""
        self.alphabet: str = ALPHABET

    def generate_barcodes_uniform(
        self, num_barcodes: int, barcode_length: int, min_hamming_dist: int, counts: int
    ) -> DataFrame:
        """
        Generate barcodes with a uniform distribution.

        Parameters:
            num_barcodes (int): Number of barcodes to generate.
            barcode_length (int): Length of each barcode.
            min_hamming_dist (int): Minimum Hamming distance between generated barcodes.
            counts (int): Count assigned to each barcode.

        Returns:
            DataFrame: A DataFrame containing barcodes, their counts, and frequencies.
        """
        barcodes = generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist)
        counts_array = np.full(num_barcodes, counts)
        frequencies = counts_array / counts_array.sum() * 100
        return pd.DataFrame({"barcode": barcodes, "count": counts_array, "frequency": frequencies})

    def generate_barcodes_random(
        self, num_barcodes: int, barcode_length: int, min_hamming_dist: int, min_frequency: float, max_count: int
    ) -> DataFrame:
        """
        Generate barcodes with a random distribution.

        Parameters:
            num_barcodes (int): Number of barcodes to generate.
            barcode_length (int): Length of each barcode.
            min_hamming_dist (int): Minimum Hamming distance between generated barcodes.
            min_frequency (float): Minimum frequency threshold for barcodes.
            max_count (int): Maximum count for each barcode.

        Returns:
            DataFrame: A DataFrame containing barcodes, their counts, and frequencies.
        """
        barcodes = generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist)
        counts = np.array([random.randint(1, max_count) for _ in range(num_barcodes)])
        total_count = np.sum(counts)
        min_count = int(min_frequency * total_count / 100)
        counts = np.clip(counts, min_count, None)
        frequencies = counts / np.sum(counts) * 100
        return pd.DataFrame({"barcode": barcodes, "count": counts, "frequency": frequencies})

    def generate_barcodes_normal(
        self,
        num_barcodes: int,
        barcode_length: int,
        min_hamming_dist: int,
        mean: float,
        std_dev: float,
        min_frequency: float,
    ) -> DataFrame:
        """
        Generate barcodes with a normal distribution.

        Parameters:
            num_barcodes (int): Number of barcodes to generate.
            barcode_length (int): Length of each barcode.
            min_hamming_dist (int): Minimum Hamming distance between generated barcodes.
            mean (float): Mean value for the normal distribution of counts.
            std_dev (float): Standard deviation for the normal distribution of counts.
            min_frequency (float): Minimum frequency threshold for barcodes.

        Returns:
            DataFrame: A DataFrame containing barcodes, their counts, and frequencies.
        """
        barcodes = generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist)
        counts = np.random.normal(loc=mean, scale=std_dev, size=num_barcodes).astype(int)
        counts = np.clip(counts, 1, None)
        total_count = np.sum(counts)
        min_count = int(min_frequency * total_count / 100)
        counts = np.clip(counts, min_count, None)
        frequencies = counts / np.sum(counts) * 100
        return pd.DataFrame({"barcode": barcodes, "count": counts, "frequency": frequencies})

    def generate_barcodes_powerlaw(
        self, num_barcodes: int, barcode_length: int, min_hamming_dist: int, exponent: float, min_frequency: float
    ) -> DataFrame:
        """
        Generate barcodes with a power-law distribution.

        Parameters:
            num_barcodes (int): Number of barcodes to generate.
            barcode_length (int): Length of each barcode.
            min_hamming_dist (int): Minimum Hamming distance between generated barcodes.
            exponent (float): Exponent for the power-law distribution.
            min_frequency (float): Minimum frequency threshold for barcodes.

        Returns:
            DataFrame: A DataFrame containing barcodes, their counts, and frequencies.
        """
        barcodes = generate_random_barcodes_with_threshold(num_barcodes, barcode_length, min_hamming_dist)
        counts = np.random.zipf(a=exponent, size=num_barcodes).astype(int)
        total_count = np.sum(counts)
        min_count = int(min_frequency * total_count / 100)
        counts = np.clip(counts, min_count, None)
        frequencies = counts / np.sum(counts) * 100
        return pd.DataFrame({"barcode": barcodes, "count": counts, "frequency": frequencies})

    def simulate_pcr_with_errors(
        self,
        barcode_df: DataFrame,
        num_pcr_chimeras: int,
        max_hamming_dist: int,
        pcr_error_rate: float,
        pcr_error_fraction: float,
    ) -> DataFrame:
        """
        Simulate the PCR process with error-prone replicates.

        Parameters:
            barcode_df (DataFrame): DataFrame containing the original barcodes.
            num_pcr_chimeras (int): Number of PCR chimeras to generate.
            max_hamming_dist (int): Maximum Hamming distance for error-prone PCR replicates.
            pcr_error_rate (float): Error rate for the PCR process.
            pcr_error_fraction (float): Fraction of the original count to use for generating PCR replicates.

        Returns:
            DataFrame: A DataFrame containing original and PCR-replicated barcodes with their counts and clusters.
        """
        if pcr_error_fraction == 0:
            return barcode_df
        all_barcodes = []
        clusters = []
        total_counts = []

        for idx, row in barcode_df.iterrows():
            barcode = row["barcode"]
            count = row["count"]

            # Add original barcode
            all_barcodes.append(barcode)
            clusters.append(f"cluster_{idx + 1}")
            total_counts.append(count)

            # Generate PCR replicates
            max_freq = max(1, int(count * pcr_error_fraction)) // num_pcr_chimeras
            if max_freq > 0:
                pcr_replicates = generate_pcr_replicates(barcode, num_pcr_chimeras, max_hamming_dist, pcr_error_rate)
                for replicate in pcr_replicates:
                    all_barcodes.append(replicate)
                    clusters.append(f"cluster_{idx + 1}")
                    total_counts.append(max_freq)

        pcr_df = pd.DataFrame({"barcode": all_barcodes, "cluster": clusters, "count": total_counts})
        cluster_totals = pcr_df.groupby("cluster")["count"].sum().reset_index().rename(columns={"count": "total_count"})
        return pcr_df.merge(cluster_totals, on="cluster")

    def extend_barcodes(self, df: DataFrame, sequence_length: int, flank_left: str, flank_right: str) -> DataFrame:
        """
        Extend barcodes with left and right flanks and random nucleotides to reach the specified sequence length.

        Parameters:
            df (DataFrame): DataFrame containing the barcodes to extend.
            sequence_length (int): The desired total length of the sequence.
            flank_left (str): Left flank to add to each barcode.
            flank_right (str): Right flank to add to each barcode.

        Returns:
            DataFrame: A DataFrame with the extended sequences as reads.
        """
        reads = []
        for _, row in df.iterrows():
            barcode = row["barcode"]
            count = row["count"]
            full_sequence = f"{flank_left}{barcode}{flank_right}"
            remaining_length = sequence_length - len(full_sequence)

            if remaining_length > 0:
                extension_each_side = remaining_length // 2
                left_extension = "".join(random.choices(self.alphabet, k=extension_each_side))
                right_extension = "".join(random.choices(self.alphabet, k=remaining_length - extension_each_side))
                full_sequence = f"{left_extension}{full_sequence}{right_extension}"

            for _ in range(count):
                reads.append(full_sequence)
        df = pd.DataFrame()
        df["read"] = reads
        return df

    def generate_fastq_generator_input(
        self,
        art: bool,
        num_barcodes: int,
        barcode_length: int,
        min_hamming_dist: int,
        num_pcr_chimeras: int,
        max_pcr_hamming_dist: int,
        pcr_error_rate: float,
        pcr_error_fraction: float,
        sequence_length: int,
        flank_left: str,
        flank_right: str,
        distribution_params: Dict[str, Any],
    ) -> Tuple[DataFrame, DataFrame]:
        """
        Generate input data for a FASTQ generator based on the specified barcode distribution and PCR parameters.

        Parameters:
            distribution_type (str): The type of distribution for generating barcodes ('uniform', 'random', 'normal', 'powerlaw').
            num_barcodes (int): Number of barcodes to generate.
            barcode_length (int): Length of each barcode.
            min_hamming_dist (int): Minimum Hamming distance between generated barcodes.
            num_pcr_chimeras (int): Number of PCR chimeras to generate.
            max_pcr_hamming_dist (int): Maximum Hamming distance for PCR replicates.
            pcr_error_rate (float): Error rate for the PCR process.
            pcr_error_fraction (float): Fraction of the original count to use for generating PCR replicates.
            sequence_length (int): The total length of the final sequence after extending with flanks.
            flank_left (str): Left flank to add to each barcode.
            flank_right (str): Right flank to add to each barcode.
            distribution_params (Dict[str, Any]): Additional parameters specific to the distribution type.

        Returns:
            Tuple[DataFrame, DataFrame]: A tuple containing the ground truth DataFrame and the extended sequences DataFrame.
        """
        if art:
            barcode_df = self.generate_barcodes_uniform(num_barcodes, barcode_length, min_hamming_dist, 1)
            ground_truth = barcode_df.loc[:, ["barcode", "count"]].drop_duplicates()
            extended_df = self.extend_barcodes(barcode_df, sequence_length, flank_left, flank_right)
            return barcode_df, extended_df[["read"]]
        distribution_type = distribution_params["distribution_type"]
        if distribution_type == "uniform":
            counts = distribution_params.get("counts", 1)
            barcode_df = self.generate_barcodes_uniform(num_barcodes, barcode_length, min_hamming_dist, counts)
        elif distribution_type == "random":
            min_frequency = distribution_params.get("min_frequency", 0.001)
            max_count = distribution_params.get("max_count", 100)
            barcode_df = self.generate_barcodes_random(
                num_barcodes, barcode_length, min_hamming_dist, min_frequency, max_count
            )
        elif distribution_type == "normal":
            mean = distribution_params.get("mean", 50.0)
            std_dev = distribution_params.get("std_dev", 10.0)
            min_frequency = distribution_params.get("min_frequency", 0.001)
            barcode_df = self.generate_barcodes_normal(
                num_barcodes, barcode_length, min_hamming_dist, mean, std_dev, min_frequency
            )
        elif distribution_type == "powerlaw":
            exponent = distribution_params.get("exponent", 1.5)
            min_frequency = distribution_params.get("min_frequency", 0.001)
            barcode_df = self.generate_barcodes_powerlaw(
                num_barcodes, barcode_length, min_hamming_dist, exponent, min_frequency
            )
        else:
            e = f"Unsupported distribution type: {distribution_type}"
            raise ValueError(e)
        pcr_df = self.simulate_pcr_with_errors(
            barcode_df, num_pcr_chimeras, max_pcr_hamming_dist, pcr_error_rate, pcr_error_fraction
        )
        ground_truth = barcode_df.loc[:, ["barcode", "count"]].drop_duplicates()
        extended_df = self.extend_barcodes(pcr_df, sequence_length, flank_left, flank_right)
        return ground_truth, extended_df[["read"]]
