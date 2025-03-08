from typing import List, Optional, Tuple, Union

from .pattern_matching import PatternMatching


def _filter(
    phred33: bool, quality: Union[str, List[int]], threshold: int, read_qc_percentage: int, bam: bool = False
) -> bool:
    """
    Filters the sequence based on the quality threshold.

    Parameters:
        quality (str): The quality string of the sequence.
        threshold (int): The quality threshold.
        read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.

    Returns:
        bool: True if the sequence passes the quality filter, False otherwise.
    """
    if threshold == 0 or read_qc_percentage == 0:
        return True
    if bam and isinstance(quality, list):
        above_threshold = sum(1 for q in quality if q >= threshold)
        return (above_threshold / len(quality)) * 100 >= read_qc_percentage
    elif not bam and isinstance(quality, str):
        threshold = threshold + 33 if phred33 else threshold + 64
        above_threshold = sum(1 for q in quality if ord(q) >= threshold)
        return (above_threshold / len(quality)) * 100 >= read_qc_percentage
    return False


def _match_sequence(
    sequence: str,
    pattern_matching: PatternMatching,
    flanked_sequence: bool,
    min_read_length: int,
    max_read_length: int,
    read_length: int,
) -> Optional[str]:
    """
    Matches a sequence using the provided pattern matching object.

    Parameters:
        sequence (str): The sequence to be matched.
        pattern_matching (PatternMatching): An object used to find sequences.
        flanked_sequence (bool): If True, only return matched sequences that are flanked.
        min_read_length (int): The minimum read length required for a match if `flanked_sequence` is True.

    Returns:
        Optional[str]: The matched sequence if found and valid, otherwise None.
    """
    if pattern_matching is None:
        return None
    matched_read = pattern_matching.find_sequence(sequence)
    if matched_read:
        read_len = len(matched_read)
        if flanked_sequence and (
            (read_length > 0 and read_len != read_length)
            or (min_read_length > 0 and read_len < min_read_length)
            or (max_read_length > 0 and read_len > max_read_length)
        ):
            return None

        return matched_read

    return None


def extract_tags_parse(header: str) -> str:
    """
    Extracts the tags from the header of a BAM file.

    Parameters:
        header (str): The header of a BAM file.

    Returns:
        Tuple[str, str, str, str]: The extracted tags.
    """
    cell_bc = header.split("__")[0]
    return cell_bc


def process_reads(
    pattern_matching: PatternMatching,
    r1: Optional[Tuple[str, str, str]] = None,
    r2: Optional[Tuple[str, str, str]] = None,
    seq: Optional[str] = None,
    quality: Optional[List[int]] = None,
    phred33: bool = True,
    threshold: int = 0,
    read_qc_percentage: int = 0,
    flanked_sequence: bool = False,
    min_read_length: int = 0,
    max_read_length: int = 0,
    read_length: int = 0,
    dna: bool = True,
    bam: bool = False,
) -> Union[str, None]:
    """
    Processes a single read pair, single-end read, or individual read (r1 or r2),
    filtering by quality and matching sequences.

    Parameters:
        pattern_matching (PatternMatching): An instance of the PatternMatching class for sequence matching.
        r1 (Optional[Tuple[str, str, str]]): The first read tuple containing name, sequence, and quality.
        r2 (Optional[Tuple[str, str, str]]): The second read tuple (if paired-end).
        seq (Optional[str]): The sequence to be matched if input is bam.
        quality (Optional[List[str]]): The quality extracted from a bam file read.
        phred33 (bool): Whether the quality scores are in Phred33 format.
        threshold (int): The quality threshold for filtering reads.
        read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.
        flanked_sequence (bool): Whether to check for flanked sequences.
        min_read_length (int): The minimum read length required for a match.
        max_read_length (int): The maximum read length required for a match.
        read_length (int): The expected read length for sequence matching.
        dna (bool): Whether the data is from DNA sequencing (False for 10x Genomics).
        bam (bool): Whether the input is coming from a BAM file.

    Returns:
        Union[str, None]: The matched sequence if a match is found and meets quality and length requirements, otherwise None.
    """
    if not bam:
        # override seq and quality if not bam
        if r1 and r2:
            combined_quality = r1[2] + r2[2]
            sequence = r1[1] + r2[1] if dna else r2[1]  # 10x
        elif r1:
            combined_quality = r1[2]
            sequence = r1[1]
        elif r2:
            combined_quality = r2[2]
            sequence = r2[1]
        else:
            e = "At least one read must be provided."
            raise ValueError(e)
        if _filter(phred33, combined_quality, threshold, read_qc_percentage, bam=bam):
            return _match_sequence(
                sequence, pattern_matching, flanked_sequence, min_read_length, max_read_length, read_length
            )
    elif seq and quality and bam:
        if _filter(phred33, quality, threshold, read_qc_percentage, bam=bam):
            return _match_sequence(
                seq, pattern_matching, flanked_sequence, min_read_length, max_read_length, read_length
            )
    return None  # Ensure all paths return a value
