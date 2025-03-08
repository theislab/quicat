import json

import pytest
from quicat.cli.utils import ValidateReference


@pytest.fixture
def valid_txt_file(tmpdir):
    file_path = tmpdir.join("valid_sequences.txt")
    sequences = ["ACGT", "GTCGAT", "TGCATGCA"]
    with open(file_path, "w") as f:
        for seq in sequences:
            f.write(f"{seq}\n")
    return str(file_path)


@pytest.fixture
def invalid_txt_file_with_n(tmpdir):
    file_path = tmpdir.join("invalid_sequences.txt")
    sequences = ["ACGTN", "GTCGAT", "TG*CA"]
    with open(file_path, "w") as f:
        for seq in sequences:
            f.write(f"{seq}\n")
    return str(file_path)


@pytest.fixture
def valid_json_file(tmpdir):
    file_path = tmpdir.join("valid_sequences.json")
    sequences = {"ACGT": "sequence_1", "GTCGAT": "sequence_2", "TGCATGCA": "sequence_3"}
    with open(file_path, "w") as f:
        json.dump(sequences, f)
    return str(file_path)


@pytest.fixture
def invalid_json_file_with_n(tmpdir):
    file_path = tmpdir.join("invalid_sequences.json")
    sequences = {"ACGTN": "sequence_1", "GTCGAT": "sequence_2", "TG*CA": "sequence_3"}
    with open(file_path, "w") as f:
        json.dump(sequences, f)
    return str(file_path)


@pytest.fixture
def valid_fasta_file(tmpdir):
    file_path = tmpdir.join("valid_sequences.fasta")
    sequences = {"sequence_1": "ACGT", "sequence_2": "GTCGAT", "sequence_3": "TGCATGCA"}
    with open(file_path, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")
    return str(file_path)


@pytest.fixture
def invalid_fasta_file_with_n(tmpdir):
    file_path = tmpdir.join("invalid_sequences.fasta")
    sequences = {"sequence_1": "ACGTN", "sequence_2": "GTCGAT", "sequence_3": "TG*CA"}
    with open(file_path, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")
    return str(file_path)


# Test valid string sequence
def test_validate_string_sequence():
    validator = ValidateReference("ACGTGTCAGT")
    assert validator.get_validated_data() == "ACGTGTCAGT"


# Test valid string sequence
def test_validate_string_sequence_star():
    validator = ValidateReference("ACGTG*CAGT")
    assert validator.get_validated_data() == "ACGTG*CAGT"


# Test valid string sequence
def test_validate_string_sequence_N():
    validator = ValidateReference("ACGTGNCAGT")
    assert validator.get_validated_data() == "ACGTGNCAGT"


# Test invalid string sequence with invalid characters
def test_invalid_string_sequence():
    with pytest.raises(ValueError, match="Sequence contains invalid characters"):
        ValidateReference("ACGTBXGT")


# Test invalid string sequence with both 'N' and '*'
def test_invalid_string_with_n_and_star():
    with pytest.raises(ValueError, match="Sequence cannot contain both 'N' and '*'"):
        ValidateReference("ACGT*N")


# Test valid string sequence with '*'
def test_valid_string_with_star():
    validator = ValidateReference("ACGT*GTCA")
    assert validator.get_validated_data() == "ACGT*GTCA"


# Test valid sequences from TXT file
def test_validate_txt_file(valid_txt_file):
    validator = ValidateReference(valid_txt_file)
    assert validator.get_validated_data() == ["ACGT", "GTCGAT", "TGCATGCA"]


# Test invalid sequences from TXT file with 'N'
def test_invalid_txt_file_with_n(invalid_txt_file_with_n):
    with pytest.raises(ValueError, match="TXT file sequences cannot contain 'N' or '*'"):
        ValidateReference(invalid_txt_file_with_n)


# Test valid sequences from JSON file
def test_validate_json_file(valid_json_file):
    validator = ValidateReference(valid_json_file)
    assert validator.get_validated_data() == {"ACGT": "sequence_1", "GTCGAT": "sequence_2", "TGCATGCA": "sequence_3"}


# Test invalid sequences from JSON file with 'N'
def test_invalid_json_file_with_n(invalid_json_file_with_n):
    with pytest.raises(ValueError, match="JSON file sequences cannot contain 'N' or '*'"):
        ValidateReference(invalid_json_file_with_n)


# Test valid sequences from FASTA file
def test_validate_fasta_file(valid_fasta_file):
    validator = ValidateReference(valid_fasta_file)
    assert validator.get_validated_data() == {"ACGT": "sequence_1", "GTCGAT": "sequence_2", "TGCATGCA": "sequence_3"}


# Test invalid sequences from FASTA file with 'N'
def test_invalid_fasta_file_with_n(invalid_fasta_file_with_n):
    with pytest.raises(ValueError, match="FASTA file sequences cannot contain 'N' or '*'"):
        ValidateReference(invalid_fasta_file_with_n)


# Test non-existent file error
def test_nonexistent_file_error():
    with pytest.raises(FileNotFoundError, match="Path does not exist: non_existent_file.txt"):
        ValidateReference("non_existent_file.txt")


# Test unsupported file type error
def test_unsupported_file_type(tmpdir):
    file_path = tmpdir.join("unsupported_file_type.csv")
    with open(file_path, "w") as f:
        f.write("sample data")
    with pytest.raises(ValueError, match="Unsupported file type"):
        ValidateReference(str(file_path))
