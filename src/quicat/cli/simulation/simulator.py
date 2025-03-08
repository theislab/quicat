import multiprocessing
import os
from typing import Any, Dict, Optional

import pandas as pd
import yaml
from joblib import Parallel, delayed

from .fastq_writer import FastqWriter
from .simulate_barcodes_bulk_dna import GenerateBarcodeSequencesDNA
from .simulate_clones_sc import GenerateBarcodeSequencesSC


class Simulator:
    """
    The Simulator class is responsible for generating synthetic sequencing data based on user-defined configurations.

    This class supports two types of simulations: DNA barcodes and single-cell RNA sequencing (scRNA-seq).
    The simulation is configured using a YAML file, which specifies parameters such as barcode distribution,
    PCR error rates, and sequence lengths.

    Attributes:
        config (Dict[str, Any]): Configuration settings loaded from the YAML file.
    """

    def __init__(self, config_file: str) -> None:
        """
        Initializes the Simulator with the provided YAML configuration file.

        Parameters:
            config_file (str): Path to the YAML configuration file.
        """
        self.config = self._load_config(config_file)
        self.n_threads = self.config["general"].get("n_threads", 1)
        available_threads = multiprocessing.cpu_count()
        if self.n_threads > available_threads:  # safety check
            self.n_threads = available_threads
        self.samples = self.config["general"]["samples"]

    def _load_config(self, config_file: str) -> Dict[str, Any]:
        """
        Load the YAML configuration file.

        Parameters:
            config_file (str): Path to the YAML configuration file.

        Returns:
            Dict[str, Any]: Loaded configuration settings as a dictionary.
        """
        with open(config_file) as file:
            return yaml.safe_load(file) or {}

    def generate_synthetic_data(self) -> None:
        """
        Run the simulation based on the configuration file.

        This method determines the type of simulation (DNA or single-cell) and delegates
        the processing to the appropriate method.
        """
        simulation_type = self.config["general"]["simulation_type"]
        simulator = self._run_sc_simulation if simulation_type == "sc" else self._run_dna_simulation

        Parallel(n_jobs=self.n_threads)(delayed(simulator)(sample) for sample in self.samples)
        return None

    def _run_dna_simulation(self, sample: str) -> None:
        """
        Run the DNA barcode simulation based on the configuration settings.

        This method generates synthetic DNA barcode sequences, simulates PCR errors,
        and prepares the data for FASTQ generation. The ground truth and FASTQ data are saved to files.
        """
        dna_simulator = GenerateBarcodeSequencesDNA()
        output_path = self.config["general"]["output_path"]
        csv_output_path = os.path.join(output_path, f"{sample}_ground_truth.csv")
        fastq_prefix = os.path.join(output_path, f"{sample}")
        distribution_params = self._get_distribution_params()

        ground_truth, fastq_input = dna_simulator.generate_fastq_generator_input(
            art=self.config["common"].get("art"),
            barcode_length=self.config["common"].get("barcode_length"),
            min_hamming_dist=self.config["common"].get("min_hamming_dist"),
            pcr_error_rate=self.config["common"].get("pcr_error_rate"),
            pcr_error_fraction=self.config["common"].get("pcr_error_fraction"),
            max_pcr_hamming_dist=self.config["common"].get("max_pcr_hamming_dist"),
            num_pcr_chimeras=self.config["common"].get("num_pcr_chimeras"),
            sequence_length=self.config["common"].get("sequence_length"),
            flank_left=self.config["common"].get("flank_left"),
            flank_right=self.config["common"].get("flank_right"),
            num_barcodes=self.config["dna"].get("num_barcodes"),
            distribution_params=distribution_params,
        )

        if self.config["common"]["art"]:
            fasta_file = os.path.join(output_path, f"{sample}.fa")
            with open(fasta_file, "w") as fasta:
                i = 1
                for seq in fastq_input["read"]:
                    fasta.write(f">barcode_{i}\n{seq}\n")
                    i += 1
        else:
            writer = FastqWriter(fastq_prefix)
            writer.write_fastq(fastq_input, is_dna=True)
        self._save_ground_truth(ground_truth, csv_output_path, is_dna=True, barcodes_csv_path=None)

    def _run_sc_simulation(self, sample: str) -> None:
        """
        Run the single-cell RNA sequencing simulation based on the configuration settings.

        This method generates synthetic single-cell RNA sequences, simulates PCR errors,
        and prepares the data for FASTQ generation. The ground truth and FASTQ data are saved to files.
        """
        sc_simulator = GenerateBarcodeSequencesSC()
        output_path = self.config["general"]["output_path"]
        csv_output_path = os.path.join(output_path, f"{sample}_ground_truth.csv")
        barcodes_csv_path = os.path.join(output_path, f"{sample}_whitelist.tsv")
        fastq_prefix = os.path.join(output_path, f"{sample}")
        distribution_params = self._get_distribution_params()

        ground_truth, fastq_input = sc_simulator.generate_fastq_generator_input(
            umis_per_cell=self.config["sc"]["umis_per_cell"],
            barcode_length=self.config["common"]["barcode_length"],
            pcr_error_rate=self.config["common"]["pcr_error_rate"],
            pcr_error_fraction=self.config["common"]["pcr_error_fraction"],
            max_pcr_hamming_dist=self.config["common"]["max_pcr_hamming_dist"],
            num_pcr_chimeras=self.config["common"]["num_pcr_chimeras"],
            sequence_length=self.config["common"]["sequence_length"],
            flank_left=self.config["common"]["flank_left"],
            flank_right=self.config["common"]["flank_right"],
            num_cells=self.config["sc"]["num_cells"],
            num_clones=self.config["sc"]["num_clones"],
            cell_bc_length=self.config["sc"]["cell_bc_length"],
            umi_length=self.config["sc"]["umi_length"],
            distribution_params=distribution_params,
        )

        self._save_ground_truth(ground_truth, csv_output_path, is_dna=False, barcodes_csv_path=barcodes_csv_path)
        writer = FastqWriter(fastq_prefix)
        writer.write_fastq(fastq_input, is_dna=False)

    def _get_distribution_params(self) -> Dict[str, Any]:
        """
        Extract distribution parameters from the configuration file based on the distribution type.

        Returns:
            Dict[str, Any]: A dictionary of distribution-specific parameters.

        Raises:
            ValueError: If the distribution type is not supported.
        """
        distribution_type = self.config["distribution_params"]["distribution_type"]
        params = {"distribution_type": distribution_type}

        if distribution_type == "uniform":
            params["counts"] = self.config["distribution_params"]["counts"]
        elif distribution_type == "random":
            params["max_count"] = self.config["distribution_params"]["max_count"]
        elif distribution_type == "normal":
            params["mean"] = self.config["distribution_params"]["mean"]
            params["std_dev"] = self.config["distribution_params"]["std_dev"]
        elif distribution_type == "powerlaw":
            params["exponent"] = self.config["distribution_params"]["exponent"]
        else:
            e = f"Unsupported distribution type: {distribution_type}"
            raise ValueError(e)

        if distribution_type in ["random", "normal", "powerlaw"]:
            params["min_frequency"] = self.config["distribution_params"]["min_frequency"]

        return params

    def _save_ground_truth(
        self,
        ground_truth: pd.DataFrame,
        output_file: str,
        is_dna: bool,
        barcodes_csv_path: Optional[str],
        use_art: bool = False,
    ) -> None:
        """
        Save the ground truth DataFrame to a CSV file.

        Parameters:
            ground_truth (pd.DataFrame): The DataFrame containing the ground truth barcodes.
            output_file (str): Path to the output CSV file.
        """
        ground_truth.to_csv(output_file, index=False)
        if not is_dna:
            if not barcodes_csv_path:
                e = "Missing barcodes CSV path for single-cell simulation"
                raise ValueError(e)
            cell_barcodes = pd.Series([f"{barcode}" for barcode in ground_truth["cell_bc"].unique()])
            cell_barcodes.to_csv(barcodes_csv_path, index=False, sep="\t", header=False)
