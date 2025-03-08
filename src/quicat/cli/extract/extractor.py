import gc
import multiprocessing
import os
from typing import Any, Dict, Union

import numpy as np
import polars as pl
import yaml
from joblib import Parallel, delayed
from rich.console import Console

from quicat import read_dna, read_sc
from quicat.cli.utils import ValidateInputs, ValidateReference

from .barcodes_processor import BarcodesProcessor
from .reads_processor_bam import ReadsProcessorBAM
from .reads_processor_fastq import ReadsProcessorFASTQ


class Extractor:
    """
    The Extractor class orchestrates the processing of sequencing reads and the extraction of barcodes
    from different input formats (BAM or FASTQ) based on the provided configuration.

    Attributes:
        config_dict (Dict[str, Any]): Configuration parameters loaded from the YAML file.
        input_csv (str): Path to the input CSV file specified in the configuration.
        output_path (str): Path to save the output files.
        input_dict (Dict[str, Any]): Dictionary containing processed input configuration data.
        console (Console): Rich Console object for printing formatted output to the console.
    """

    def __init__(self, config_path: str) -> None:
        """
        Initializes the Extractor object by loading and validating the configuration and inputs.

        Parameters:
            config_path (str): Path to the YAML configuration file.
        """
        # Load the YAML configuration file
        config = self._load_config(config_path)
        self.console = Console()

        # Extract sections from the YAML config
        self.config_dict = config.get("config", {})
        self.input_csv = config.get("input_csv", "")
        self.output_path = config.get("output_path", "")
        validator_reference = ValidateReference(self.config_dict["reference"])
        self.config_dict["reference"] = validator_reference.get_validated_data()
        validator_input = ValidateInputs(self.input_csv, self.config_dict)
        self.input_dict = validator_input.build_config()

    @staticmethod
    def _load_config(config_path: str) -> Dict[str, Any]:
        """
        Loads the YAML configuration file.

        Parameters:
            config_path (str): Path to the YAML configuration file.

        Returns:
            Dict[str, Any]: Loaded configuration data as a dictionary.
        """
        with open(config_path) as file:
            return yaml.safe_load(file) or {}

    @staticmethod
    def _process_sample(sample: str, config_dict: dict, input_dict: dict, output_path: str) -> pl.DataFrame:
        """
        Processes a single sample based on the provided configuration and input data.

        Parameters:
            sample (str): Name of the sample to process.
            config_dict (dict): Configuration dictionary for the Extractor run.
            input_dict (dict): Input dictionary containing sample-specific parameters.
            output_path (str): Path to save the intermediate output files.

        Returns:
            pl.DataFrame: DataFrame containing processed reads data for the sample.
        """
        _console = Console()
        try:
            reads_processor: Union[ReadsProcessorBAM, ReadsProcessorFASTQ]
            replicate = input_dict[sample].get("replicate", None)
            condition = input_dict[sample].get("condition", None)
            if config_dict["input"] == "fastq":
                reads_processor = ReadsProcessorFASTQ(
                    paired_end=input_dict[sample]["paired_end"],
                    read_qc_threshold=input_dict[sample]["read_qc_threshold"],
                    read_qc_percentage=input_dict[sample]["read_qc_percentage"],
                    phred33=input_dict[sample]["phred33"],
                    fastq_path_r1=input_dict[sample].get("fastq_path_r1", None),
                    fastq_path_r2=input_dict[sample].get("fastq_path_r2", None),
                    reference=input_dict[sample]["reference"],
                    flanked_pattern=input_dict[sample]["flanked_pattern"],
                    min_read_length=input_dict[sample].get("min_read_length", 0),
                    max_read_length=input_dict[sample].get("max_read_length", 0),
                    read_length=input_dict[sample].get("read_length", 0),
                    sequencing_technology=input_dict[sample]["sequencing_technology"],
                    barcode_path=input_dict[sample].get("barcodes_path", None),
                    aln_mismatches=input_dict[sample].get("aln_mismatches", None),
                    flanking_mismatches=input_dict[sample].get("flanking_mismatches", None),
                    left_flanking_coverage=input_dict[sample].get("left_flanking_coverage", None),
                    right_flanking_coverage=input_dict[sample].get("right_flanking_coverage", None),
                )
            elif config_dict["input"] == "bam":
                reads_processor = ReadsProcessorBAM(
                    read_qc_threshold=input_dict[sample]["read_qc_threshold"],
                    read_qc_percentage=input_dict[sample]["read_qc_percentage"],
                    phred33=input_dict[sample]["phred33"],
                    bam_path=input_dict[sample]["bam_path"],
                    reference=input_dict[sample]["reference"],
                    flanked_pattern=input_dict[sample]["flanked_pattern"],
                    min_read_length=input_dict[sample].get("min_read_length", 0),
                    max_read_length=input_dict[sample].get("max_read_length", 0),
                    read_length=input_dict[sample].get("read_length", 0),
                    sequencing_technology=input_dict[sample]["sequencing_technology"],
                    barcode_path=input_dict[sample].get("barcodes_path", None),
                    contig=input_dict[sample].get("contig", "*"),  # default to unmapped reads
                    aln_mismatches=input_dict[sample].get("aln_mismatches", None),
                    flanking_mismatches=input_dict[sample].get("flanking_mismatches", None),
                    left_flanking_coverage=input_dict[sample].get("left_flanking_coverage", None),
                    right_flanking_coverage=input_dict[sample].get("right_flanking_coverage", None),
                )
            reads = reads_processor.extract_barcodes()

            if input_dict[sample]["sequencing_technology"] == "DNA" and isinstance(reads, list):
                unique_reads, counts = np.unique(reads, return_counts=True)
                df_sample_reads = pl.DataFrame(
                    {
                        "sample": sample,
                        "barcode": [str(read) for read in unique_reads],
                        "counts": counts,
                    }
                )
            elif isinstance(reads, dict):
                df_sample_reads = pl.DataFrame(
                    {
                        "sample": sample,
                        "cell": list(reads),
                        "barcode": [
                            list(barcodes.keys()) if isinstance(barcodes, dict) else [] for barcodes in reads.values()
                        ],
                        "counts": [
                            list(barcodes.values()) if isinstance(barcodes, dict) else [] for barcodes in reads.values()
                        ],
                    }
                )

            # Add replicate and condition columns if available
            if replicate is not None:
                df_sample_reads = df_sample_reads.with_columns(pl.lit(replicate).alias("replicate"))
            if condition is not None:
                df_sample_reads = df_sample_reads.with_columns(pl.lit(condition).alias("condition"))

            if input_dict[sample]["sequencing_technology"] != "DNA":
                # Save the intermediate output to a CSV file
                output_filename = os.path.join(output_path, f"{sample}_output.csv")
                df_sample_reads = df_sample_reads.explode(["barcode", "counts"])
                df_sample_reads.write_csv(output_filename)

                # Add sample suffix to 'cell' column if present
                df_sample_reads = df_sample_reads.with_columns((pl.col("cell") + f"_{sample}").alias("cell"))

            _console.print(f"[green]Sample {sample} processed successfully[/green] ✅")
            _console.print(f"[yellow]Number of barcodes extracted: {df_sample_reads.shape[0]}[/yellow]")
            gc.collect()
        except Exception as e:
            _console.print(f"[red]Sample {sample} failed: {e}[/red] 🚨\n")
            return pl.DataFrame()

        return df_sample_reads

    def run(self) -> None:
        """
        Runs the Extractor pipeline to process all samples, extract barcodes, filter them, and save the results.

        This method manages the entire workflow from loading the samples, processing them in parallel,
        filtering barcodes, and outputting the results to multiple CSV files.
        """
        n_threads: int = self.config_dict["n_threads"]
        available_threads = multiprocessing.cpu_count()
        if n_threads > available_threads:  # safety check
            n_threads = available_threads
        self.console.print(f"\n[bold green]Starting Extractor run with {len(self.input_dict)} samples[/bold green]🚀\n")
        self.console.print("[yellow]Configuration:[/yellow]\n")
        for key, value in self.config_dict.items():
            self.console.print(f"[yellow]{key}:[/yellow]", f"[green]{value}[/green]")
        process_sample = self._process_sample

        if n_threads == 1:
            df_samples_reads = pl.DataFrame()
            for sample in self.input_dict:
                df_sample_reads = process_sample(sample, self.config_dict, self.input_dict, self.output_path)
                df_samples_reads = df_samples_reads.vstack(df_sample_reads)
        else:
            results = Parallel(n_jobs=n_threads)(
                delayed(process_sample)(sample, self.config_dict, self.input_dict, self.output_path)
                for sample in self.input_dict
            )
            df_samples_reads = pl.concat(results)
            results = None  # free memory

        gc.collect()
        self.console.print(
            "\n[bold green]All samples processed.[/bold green]✅\n"
            "[bold green]Filtering barcodes, and building dataframe.[/bold green]\n"
        )

        filter_barcodes_raw_numbers = self.config_dict.get("filter_barcodes_raw_numbers", 0)
        filter_barcodes_relative_abundance = self.config_dict.get("filter_barcodes_relative_abundance", 0)
        distance_threshold = self.config_dict.get("distance_threshold", 0)
        n_threads = self.config_dict.get("n_threads", -1)
        barcode_ratio = self.config_dict.get("barcode_ratio", 1)
        sequencing_technology = self.config_dict["sequencing_technology"]

        barcode_filter = BarcodesProcessor(
            sequencing_technology=sequencing_technology,
            filter_barcodes_raw_numbers=filter_barcodes_raw_numbers,
            filter_barcodes_relative_abundance=filter_barcodes_relative_abundance,
            distance_threshold=distance_threshold,
            n_threads=n_threads,
            barcode_ratio=barcode_ratio,
        )

        df = barcode_filter.build_dataframe(df_samples_reads)
        if not isinstance(df, pl.DataFrame):
            e = "df is not an instance of pl.DataFrame"
            raise TypeError(e)
        gc.collect()

        print(df)
        out_csv = os.path.join(self.output_path, "barcodes_output.csv")
        df.write_csv(out_csv)
        self.console.print(f"\n[bold green]Barcodes filtered.\nCSV Output saved to {out_csv}[/bold green]✅\n")
        out_adata = os.path.join(self.output_path, "barcodes_output.h5ad")
        adata = read_dna(out_csv) if sequencing_technology == "DNA" else read_sc(out_csv)
        adata.write_h5ad(out_adata)
        self.console.print(f"\n[bold green]AnnData object saved to {out_adata}[/bold green]✅\n")
        return None
