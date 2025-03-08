import gzip
from collections import defaultdict
from typing import Dict, List, Optional, Union

import pyfastx

from .pattern_matching import PatternMatching
from .reads_utils import process_reads


class ReadsProcessorFASTQ:
    """
    A class to process sequencing reads and extract barcodes based on the sequencing technology.

    Attributes:
        paired_end (bool): Whether the sequencing is paired-end.
        read_qc_threshold (int): The quality threshold for filtering reads.
        read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.
        phred33 (bool): Whether the quality scores are in Phred33 format.
        fastq1 (pyfastx.Fastx): The first FASTQ file.
        fastq2 (Optional[pyfastx.Fastx]): The second FASTQ file (if paired-end).
        pattern_matching (PatternMatching): An instance of the PatternMatching class.
        flanked_sequence (bool): Whether to check for flanked sequences.
        min_read_length (int): The minimum read length required for a match.
        max_read_length (int): The maximum read length required for a match.
        read_length (int): The read length required for a match.
        sequencing_technology (str): The sequencing technology being used.
        barcode_path (Optional[str]): Path to the barcode file for 10x Genomics sequencing.
    """

    def __init__(
        self,
        paired_end: bool,
        read_qc_threshold: int,
        read_qc_percentage: int,
        phred33: bool,
        fastq_path_r1: str,
        fastq_path_r2: Optional[str],
        reference: Union[str, List[str], Dict[str, str]],
        flanked_pattern: bool,
        min_read_length: int,
        max_read_length: int,
        read_length: int,
        sequencing_technology: str,
        barcode_path: Optional[str] = None,
        aln_mismatches: Optional[int] = None,
        flanking_mismatches: Optional[float] = None,
        left_flanking_coverage: Optional[int] = None,
        right_flanking_coverage: Optional[int] = None,
    ) -> None:
        """
        Initializes the ReadsProcessor with the given parameters.

        Parameters:
            paired_end (bool): Whether the sequencing is paired-end.
            read_qc_threshold (int): The quality threshold for filtering reads.
            read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.
            phred33 (bool): Whether the quality scores are in Phred33 format.
            fastq_path_r1 (str): Path to the first FASTQ file.
            fastq_path_r2 (Optional[str]): Path to the second FASTQ file (if paired-end).
            reference (Union[str, List[str], Dict[str, str]]): The reference pattern(s) for sequence matching.
            flanked_pattern (bool): Whether to check for flanked sequences.
            min_read_length (int): The minimum read length required for a match.
            sequencing_technology (str): The sequencing technology being used.
            barcode_path (Optional[str]): Path to the barcode file for 10x Genomics sequencing.
            aln_mismatches (Optional[int]): The maximum number of mismatches allowed in the alignment. If not None, edlib is used, otherwise Aho-Corasick is used.
            flanking_mismatches (Optional[float]): The maximum fraction of mismatches allowed in flanking sequences.
            left_flanking_coverage (Optional[int]): The minimum overlap required for the left flanking sequence.
            right_flanking_coverage (Optional[int]): The minimum overlap required for the right flanking sequence.
        """
        self.paired_end = paired_end
        self.read_qc_threshold = read_qc_threshold
        self.read_qc_percentage = read_qc_percentage
        self.phred33 = phred33
        if paired_end:
            self.fastq1 = pyfastx.Fastx(fastq_path_r1)
            self.fastq2 = pyfastx.Fastx(fastq_path_r2)
        else:
            self.fastq1 = pyfastx.Fastx(fastq_path_r1)
            self.fastq2 = None
        self.flanked_sequence = flanked_pattern
        self.min_read_length = min_read_length
        self.max_read_length = max_read_length
        self.read_length = read_length
        self.sequencing_technology = sequencing_technology
        self.barcode_path = barcode_path
        self.pattern_matching = PatternMatching(
            pattern=reference,
            flanked_pattern=flanked_pattern,
            read_length=read_length if read_length and read_length > 0 else None,
            max_read_length=max_read_length if max_read_length and max_read_length > 0 else None,
            flanking_mismatches=flanking_mismatches,
            left_flanking_coverage=left_flanking_coverage,
            right_flanking_coverage=right_flanking_coverage,
            aln_mismatches=aln_mismatches,
        )

    def _process_dna(self) -> List[str]:
        """
        Processes single-end DNA sequencing reads, filtering by quality and matching sequences.

        Returns:
            List[str]: A list of matched sequences.
        """
        matched_reads = []
        # Local variable assignments to avoid repeated lookups
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        phred33 = self.phred33
        threshold = self.read_qc_threshold
        read_qc_percentage = self.read_qc_percentage
        pattern_matching = self.pattern_matching
        flanked_sequence = self.flanked_sequence
        min_read_length = self.min_read_length
        max_read_length = self.max_read_length
        read_length = self.read_length
        if not self.fastq2:
            for read in fastq1:
                matched_read = process_reads(
                    pattern_matching=pattern_matching,
                    r1=read,
                    r2=None,
                    phred33=phred33,
                    threshold=threshold,
                    read_qc_percentage=read_qc_percentage,
                    flanked_sequence=flanked_sequence,
                    min_read_length=min_read_length,
                    max_read_length=max_read_length,
                    read_length=read_length,
                )
                if matched_read:
                    matched_reads.append(matched_read)
        else:
            for r1, r2 in zip(fastq1, fastq2):
                matched_read = process_reads(
                    pattern_matching=pattern_matching,
                    r1=r1,
                    r2=r2,
                    phred33=phred33,
                    threshold=threshold,
                    read_qc_percentage=read_qc_percentage,
                    flanked_sequence=flanked_sequence,
                    min_read_length=min_read_length,
                    max_read_length=max_read_length,
                    read_length=read_length,
                )
                if matched_read:
                    matched_reads.append(matched_read)
        return matched_reads

    def _process_10x(self) -> Dict[str, Dict[str, int]]:
        """
        Processes 10x Genomics sequencing reads, filtering by quality and matching sequences.
        Returns a nested dictionary with cell barcodes, matched sequences, and their corresponding UMI counts.

        Returns:
            Dict[str, Dict[str, int]]: A nested dictionary where keys are cell barcodes, and values are dictionaries
            with matched sequences as keys and their UMI counts as values.
        """
        counts: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        filtered_cell_barcodes = set()

        # Local variable assignments to avoid repeated lookups
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        phred33 = self.phred33
        threshold = self.read_qc_threshold
        read_qc_percentage = self.read_qc_percentage
        pattern_matching = self.pattern_matching
        flanked_sequence = self.flanked_sequence
        min_read_length = self.min_read_length
        max_read_length = self.max_read_length
        read_length = self.read_length
        barcode_path = self.barcode_path
        suffix = "-1"
        if barcode_path:
            with gzip.open(barcode_path, "rt") as barcode_file:
                filtered_cell_barcodes = set(barcode_file.read().splitlines())
                suffix = "-" + next(iter(filtered_cell_barcodes)).split("-")[1]
        for r1, r2 in zip(fastq1, fastq2):
            cell_bc = r1[1][:16] + suffix
            if not filtered_cell_barcodes or cell_bc in filtered_cell_barcodes:
                umi: str = r1[1][16:24]
                matched_read = process_reads(
                    r1=r1,
                    r2=r2,
                    phred33=phred33,
                    threshold=threshold,
                    read_qc_percentage=read_qc_percentage,
                    pattern_matching=pattern_matching,
                    flanked_sequence=flanked_sequence,
                    min_read_length=min_read_length,
                    max_read_length=max_read_length,
                    read_length=read_length,
                    dna=False,  # 10x specific
                )
                if matched_read:
                    counts[cell_bc][matched_read][umi] += 1

        final_counts = {
            cell_bc: {barcode: len(umis) for barcode, umis in barcodes.items()} for cell_bc, barcodes in counts.items()
        }
        return final_counts

    def extract_barcodes(self) -> Union[List[str], Dict[str, Dict[str, int]]]:
        """
        Extracts barcodes based on the sequencing technology and read type.

        Returns:
            Union[List[str], Dict[str, Dict[str, int]]]: The extracted barcodes as a list of matched sequences
            (for DNA) or a nested dictionary (for 10x Genomics).
        """
        if self.sequencing_technology == "DNA":
            return self._process_dna()
        elif self.sequencing_technology == "10x":
            return self._process_10x()
        else:
            e = f"Unsupported sequencing technology for fastq files: {self.sequencing_technology}\nSelect either 'DNA' or '10x'."
            raise ValueError(e)
