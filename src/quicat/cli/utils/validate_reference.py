import json
import os
from pathlib import Path
from typing import Dict, List, Union


class ValidateReference:
    """
    A class to validate biological sequences provided as strings or in various file formats (TXT, JSON, FASTA).

    Attributes:
        input_data (str): The input data, either a sequence string or a file path.
        validated_data (Union[str, List[str], Dict[str, str]]): The validated sequence(s) in the format determined by the input.
    """

    def __init__(self, input_data: str):
        """
        Initialize the ValidateReference class with the input data.

        Parameters:
            input_data (str): A sequence string or a file path.

        Raises:
            ValueError: If the input sequence or file content is invalid.
            FileNotFoundError: If the provided file path does not exist.
        """
        self.input_data = input_data
        self.validated_data = self._validate_input()

    def _validate_input(self) -> Union[str, List[str], Dict[str, str]]:
        """
        Validate the input data, either as a string sequence or a file.

        Returns:
            Union[str, List[str], Dict[str, str]]: Validated sequence(s) based on the input type.

        Raises:
            ValueError: If the sequence contains invalid characters or improper format.
            FileNotFoundError: If the file does not exist.
        """
        if self._is_path(self.input_data):
            return self._validate_file(Path(self.input_data))
        elif "." in self.input_data:
            e = f"Path does not exist: {self.input_data}"
            raise FileNotFoundError(e)
        else:
            return self._validate_string(self.input_data)

    def _is_path(self, input_data: str) -> bool:
        """
        Check if the input data is a file path.

        Parameters:
            input_data (str): The input data to check.

        Returns:
            bool: True if the input is a valid file path, False otherwise.
        """
        return os.path.isfile(input_data)

    def _validate_string(self, sequence: str) -> str:
        """
        Validate a sequence string.

        Parameters:
            sequence (str): The sequence string to validate.

        Returns:
            str: The validated sequence string.

        Raises:
            ValueError: If the sequence contains invalid characters or improper format.
        """
        valid_chars = {"A", "C", "G", "T", "N", "*"}
        if not all(c in valid_chars for c in sequence):
            e = "Sequence contains invalid characters"
            raise ValueError(e)

        if "N" in sequence and "*" in sequence:
            e = "Sequence cannot contain both 'N' and '*'"
            raise ValueError(e)

        if "*" in sequence:
            flanking_sequences = sequence.split("*")

            # Ensure there are at most 2 flanking sequences
            if len(flanking_sequences) > 2:
                e = "Sequence with '*' must be in the format 'flanking_left*flanking_right', '*flanking_right', or 'flanking_left*'"
                raise ValueError(e)

            if len(flanking_sequences) == 2:
                left, right = flanking_sequences
                if not left and not right:
                    e = "Invalid sequence: '*' cannot appear on its own"
                    raise ValueError(e)

        return sequence

    def _validate_file(self, file_path: Path) -> Union[List[str], Dict[str, str]]:
        """
        Validate a file containing sequences in TXT, JSON, or FASTA format.

        Parameters:
            file_path (Path): The path to the file to validate.

        Returns:
            Union[List[str], Dict[str, str]]: Validated sequences from the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file format is unsupported or the sequences are invalid.
        """
        if file_path.suffix == ".txt":
            return self._validate_txt(file_path)
        elif file_path.suffix == ".json":
            return self._validate_json(file_path)
        elif file_path.suffix == ".fasta" or file_path.suffix == ".fa":
            return self._validate_fasta(file_path)
        else:
            e = "Unsupported file type"
            raise ValueError(e)

    def _validate_txt(self, file_path: Path) -> List[str]:
        """
        Validate sequences from a TXT file.

        Parameters:
            file_path (Path): The path to the TXT file.

        Returns:
            List[str]: A list of validated sequences.

        Raises:
            ValueError: If the sequences contain 'N' or '*' or if the sequence format is invalid.
        """
        sequences = []
        with file_path.open("r") as file:
            for line in file:
                sequence = line.strip()
                self._validate_string(sequence)
                if "N" in sequence or "*" in sequence:
                    e = "TXT file sequences cannot contain 'N' or '*'"
                    raise ValueError(e)
                sequences.append(sequence)
        return sequences

    def _validate_json(self, file_path: Path) -> Dict[str, str]:
        """
        Validate sequences from a JSON file.

        Parameters:
            file_path (Path): The path to the JSON file.

        Returns:
            Dict[str, str]: A dictionary of validated sequences mapped to their names.

        Raises:
            ValueError: If the sequences contain 'N' or '*' or if the sequence format is invalid.
        """
        with file_path.open("r") as file:
            data = json.load(file)
            validated_data = {}
            for sequence, name in data.items():
                self._validate_string(sequence)
                if "N" in sequence or "*" in sequence:
                    e = "JSON file sequences cannot contain 'N' or '*'"
                    raise ValueError(e)
                validated_data[sequence] = name
            return validated_data

    def _validate_fasta(self, file_path: Path) -> Dict[str, str]:
        """
        Validate sequences from a FASTA file.

        Parameters:
            file_path (Path): The path to the FASTA file.

        Returns:
            Dict[str, str]: A dictionary of validated sequences mapped to their headers.

        Raises:
            ValueError: If the sequences contain 'N' or '*' or if the sequence format is invalid.
        """
        sequences = {}
        with file_path.open("r") as file:
            header = None
            sequence = ""
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        self._validate_string(sequence)
                        if "N" in sequence or "*" in sequence:
                            e = "FASTA file sequences cannot contain 'N' or '*'"
                            raise ValueError(e)
                        sequences[sequence] = header
                    header = line[1:]
                    sequence = ""
                else:
                    sequence += line
            if header and sequence:
                self._validate_string(sequence)
                if "N" in sequence or "*" in sequence:
                    e = "FASTA file sequences cannot contain 'N' or '*'"
                    raise ValueError(e)
                sequences[sequence] = header
        return sequences

    def get_validated_data(self) -> Union[str, List[str], Dict[str, str]]:
        """
        Retrieve the validated reference data.

        Returns:
            Union[str, List[str], Dict[str, str]]: The validated sequence(s).
        """
        return self.validated_data
