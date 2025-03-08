import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd


class ValidationError(Exception):
    """Custom exception for validation errors."""

    pass


class ColumnValidationStrategy:
    """Base class for column validation strategy.

    This class defines the interface for validating columns based on the
    configuration provided.

    Attributes:
        config_dict (Dict[str, Union[str, bool, int, float]]): Configuration dictionary with settings.
    """

    def __init__(self, config_dict: Dict[str, Union[str, bool, int, float]]) -> None:
        """Initializes ColumnValidationStrategy with the given configuration.

        Parameters:
            config_dict (Dict[str, Union[str, bool, int, float]]): Configuration dictionary.
        """
        self.config_dict = config_dict

    def _get_required_columns(self) -> List[str]:
        """Returns a list of required columns based on the configuration.

        This method should be implemented by subclasses.

        Returns:
            List[str]: A list of required column names.

        Raises:
            NotImplementedError: If the method is not implemented by a subclass.
        """
        raise NotImplementedError("Must implement get_required_columns")

    def validate_columns(self, df: pd.DataFrame) -> None:
        """Validates the presence and integrity of required columns in the DataFrame.

        Parameters:
            df (pd.DataFrame): The DataFrame to validate.

        Raises:
            ValidationError: If any required column is missing or contains missing values.
        """
        required_columns = self._get_required_columns()

        for col in required_columns:
            if col not in df.columns:
                e = f"Column {col} not found in CSV file, please follow the template."
                raise ValidationError(e)
            if df[col].isna().any():
                e = f"Column {col} contains missing values, please check your data."
                raise ValidationError(e)


class DNACsvValidation(ColumnValidationStrategy):
    """Validation strategy for DNA sequencing technology."""

    def _get_required_columns(self) -> List[str]:
        """Returns the required columns for DNA sequencing based on the configuration.

        Returns:
            List[str]: A list of required column names.
        """
        if self.config_dict.get("paired_end", False):
            return ["sample", "fastq_path_r1", "fastq_path_r2"]
        return ["sample", "fastq_path_r1"]


class SingleCellFastqCsvValidation(ColumnValidationStrategy):
    """Validation strategy for 10x sequencing technology."""

    def _get_required_columns(self) -> List[str]:
        """Returns the required columns for 10x sequencing based on the configuration.

        Returns:
            List[str]: A list of required column names.
        """
        cols = ["sample", "fastq_path_r1", "fastq_path_r2"]
        if self.config_dict.get("cells_barcodes_whitelist", False):
            cols.append("barcodes_path")
        return cols


class BamCsvValidation(ColumnValidationStrategy):
    """Validation strategy for BAM input files."""

    def _get_required_columns(self) -> List[str]:
        """Returns the required columns for BAM input files.

        Returns:
            List[str]: A list of required column names.
        """
        cols = ["sample", "bam_path"]
        if self.config_dict.get("cells_barcodes_whitelist", False):
            cols.append("barcodes_path")
        return cols


class CsvValidatorFactory:
    """Factory for creating column validation strategies based on the configuration."""

    @staticmethod
    def _create_validator(config_dict: Dict[str, Union[str, bool, int, float]]) -> ColumnValidationStrategy:
        """Creates and returns the appropriate validation strategy based on the configuration.

        Parameters:
            config_dict (Dict[str, Union[str, bool, int, float]]): Configuration dictionary.

        Returns:
            ColumnValidationStrategy: The appropriate validation strategy.

        Raises:
            ValueError: If the configuration is not supported.
        """
        tech = config_dict["sequencing_technology"]
        input_type = config_dict["input"]

        if tech == "DNA":
            return DNACsvValidation(config_dict)
        elif tech == "10x" and input_type == "fastq":
            return SingleCellFastqCsvValidation(config_dict)
        elif input_type == "bam":
            return BamCsvValidation(config_dict)
        else:
            e = f"Unsupported configuration: {tech} with input type {input_type}"
            raise ValueError(e)


class ValidateInputs:
    """Class for validating and processing input CSV files for sequencing data.

    This class handles the validation of input CSV files based on the specified
    configuration and returns a validated DataFrame. It also builds a final
    configuration dictionary for further processing.

    Attributes:
        csv_file (str): Path to the CSV file to validate.
        config_dict (Dict[str, Union[str, bool, int, float]]): Configuration dictionary.
        validator (ColumnValidationStrategy): The validation strategy used for this instance.
    """

    def __init__(self, csv_file: str, config_dict: Dict[str, Union[str, bool, int, float]]) -> None:
        """Initializes the ValidateInputs class with the CSV file and configuration.

        Parameters:
            csv_file (str): Path to the CSV file to validate.
            config_dict (Dict[str, Union[str, bool, int, float]]): Configuration dictionary.
        """
        self.csv_file = csv_file
        self.config_dict = config_dict
        self.validator = CsvValidatorFactory._create_validator(config_dict)

    @staticmethod
    def _validate_path(path: str) -> None:
        """Validates if a given path exists.

        Parameters:
            path (str): The path to validate.

        Raises:
            ValidationError: If the path does not exist.
        """
        if pd.notna(path) and not os.path.exists(path):
            e = f"Path does not exist: {path}"
            raise ValidationError(e)

    @staticmethod
    def _validate_sample(sample: str) -> None:
        """Validates that the sample is a non-empty string.

        Parameters:
            sample (str): The sample value to validate.

        Raises:
            ValidationError: If the sample is not a non-empty string.
        """
        if pd.isna(sample) or not isinstance(sample, str) or sample.strip() == "":
            e = "Sample must be a non-empty string"
            raise ValidationError(e)

    def _validate(self) -> Tuple[pd.DataFrame, bool, bool]:
        """Validates the input CSV file against the expected columns and checks if paths exist.

        Returns:
            pd.DataFrame: The validated DataFrame.

        Raises:
            ValidationError: If any validation check fails.
        """
        df = pd.read_csv(self.csv_file, keep_default_na=False, na_values=[None])
        # Validate columns based on the strategy
        self.validator.validate_columns(df)

        # Validate sample names
        df["sample"].apply(self._validate_sample)

        use_condition = df["condition"].notna().any()
        use_replicate = df["replicate"].notna().any()

        if use_condition or use_replicate:  # Check if condition and replicate are partially filled
            for col in ["condition", "replicate"]:
                if df[col].notna().any() and df[col].isna().any():
                    e = f"Column {col} is partially filled. If any value is filled, all values must be filled."
                    raise ValidationError(e)

        # Validate paths for required columns
        for col in self.validator._get_required_columns():
            if col not in ["sample", "condition", "replicate"]:
                df[col].apply(self._validate_path)

        return df, use_condition, use_replicate

    def _add_metadata(self, row: pd.Series) -> Dict[str, Optional[Union[str, int]]]:
        """Checks if the 'condition' and 'replicate' columns are filled and adds them to the final config.

        Parameters:
            row (pd.Series): A row of the DataFrame.

        Returns:
            Dict[str, Optional[Union[str, int]]]: The updated sample dictionary.
        """
        # Initialize the dictionary with the expected type
        sample_dict: Dict[str, Optional[Union[str, int]]] = {}

        for key in row.index:
            value = row[key]
            sample_dict[key] = value if pd.notna(value) else None

        if pd.notna(row["condition"]) and row["condition"] != "":
            sample_dict["condition"] = str(row["condition"])
        if pd.notna(row["replicate"]) and row["replicate"] != "":
            sample_dict["replicate"] = str(row["replicate"])

        return sample_dict

    def build_config(self) -> Dict[str, Dict[str, Optional[Union[str, int]]]]:
        """Builds a final configuration dictionary by merging input CSV data with the configuration settings.

        Returns:
            Dict[str, Dict[str, Optional[Union[str, int]]]]: A dictionary where keys are sample names and values are settings for each sample.
        """
        final_config: defaultdict[str, Dict[str, Optional[Union[str, int]]]] = defaultdict(dict)
        input_df, use_condition, use_replicate = self._validate()

        for _, row in input_df.iterrows():
            sample = row["sample"]
            sample_dict = (
                self._add_metadata(row)
                if (use_condition or use_replicate)
                else {key: value for key, value in row.items() if key not in ["condition", "replicate"]}
            )
            config_dict = {
                key: value
                for key, value in self.config_dict.items()
                if isinstance(value, (str, int, dict, list, type(None)))  # For type checking
            }
            keys_to_check = ["fastq_path_r1", "fastq_path_r2", "barcodes_path", "bam_path"]
            for key in keys_to_check:
                sample_dict.setdefault(key, None)
            sample_dict.update(config_dict)
            final_config[sample] = sample_dict

        return dict(final_config)
