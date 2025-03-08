import gzip
from collections import defaultdict
from typing import Dict, List, Optional, Union

import pysam

from .pattern_matching import PatternMatching
from .reads_processor_bam_conf import BAM_TAGS
from .reads_utils import extract_tags_parse, process_reads


class ReadsProcessorBAM:
    """
    Processes sequencing reads from BAM files and extracts barcodes based on the sequencing technology.

    Attributes:
        read_qc_threshold (int): The quality threshold for filtering reads.
        read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.
        bam_file (pysam.AlignmentFile): The BAM file object for reading sequencing data.
        pattern_matching (PatternMatching): An instance of the PatternMatching class for sequence matching.
        flanked_sequence (bool): Whether to check for flanked sequences.
        min_read_length (int): The minimum read length required for a match.
        max_read_length (int): The maximum read length required for a match.
        read_length (int): The expected read length for sequence matching.
        sequencing_technology (str): The sequencing technology being used.
        phred33 (bool): Whether the quality scores are in Phred33 format.
        tags (Dict[str, str]): Dictionary of BAM tags specific to the sequencing technology.
        coordinate_based (bool): Whether the sequencing technology uses coordinate-based tags.
        contig (str): The contig to fetch reads from.
        barcode_path (Optional[str]): Path to a file containing cell barcodes to filter reads.
    """

    def __init__(
        self,
        read_qc_threshold: int,
        read_qc_percentage: int,
        phred33: bool,
        bam_path: str,
        reference: Union[str, List[str], Dict[str, str]],
        flanked_pattern: bool,
        min_read_length: int,
        max_read_length: int,
        read_length: int,
        sequencing_technology: str,
        contig: str,
        barcode_path: Optional[str] = None,
        aln_mismatches: Optional[int] = None,
        flanking_mismatches: Optional[float] = None,
        left_flanking_coverage: Optional[int] = None,
        right_flanking_coverage: Optional[int] = None,
    ) -> None:
        """
        Initializes the ReadsProcessorBAM with the provided parameters.

        Parameters:
            read_qc_threshold (int): The quality threshold for filtering reads.
            read_qc_percentage (int): The minimum percentage of bases that must meet the quality threshold.
            phred33 (bool): Whether the quality scores are in Phred33 format.
            bam_path (str): Path to the BAM file.
            reference (Union[str, List[str], Dict[str, str]]): The reference pattern(s) for sequence matching.
            flanked_pattern (bool): Whether to check for flanked sequences.
            min_read_length (int): The minimum read length required for a match.
            max_read_length (int): The maximum read length required for a match.
            read_length (int): The expected read length for sequence matching.
            sequencing_technology (str): The sequencing technology being used.
            contig (str): The contig to fetch reads from.
            barcode_path (Optional[str]): Path to a file containing cell barcodes to filter reads.
            aln_mismatches (Optional[int]): The maximum number of mismatches allowed in the alignment. If not None, edlib is used, otherwise Aho-Corasick is used.
            flanking_mismatches (Optional[float]): The maximum fraction of mismatches allowed in flanking sequences.
            left_flanking_coverage (Optional[int]): The minimum overlap required for the left flanking sequence.
            right_flanking_coverage (Optional[int]): The minimum overlap required for the right flanking sequence.

        Raises:
            ValueError: If the sequencing technology is unsupported or tags are not defined.
        """
        self.read_qc_threshold = read_qc_threshold
        self.read_qc_percentage = read_qc_percentage
        self.bam_file = pysam.AlignmentFile(bam_path, "rb")
        self.flanked_sequence = flanked_pattern
        self.min_read_length = min_read_length
        self.max_read_length = max_read_length
        self.read_length = read_length
        self.sequencing_technology = sequencing_technology
        self.phred33 = phred33
        self.contig = contig
        tags = BAM_TAGS.get(sequencing_technology)
        if tags is None:
            e = f"Tags for the sequencing technology '{sequencing_technology}' are not defined."
            raise ValueError(e)
        self.tags: Dict[str, str] = tags
        self.barcode_path = barcode_path
        self.coordinate_based = self._is_coordinate_based(self.tags)
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

    @staticmethod
    def _is_coordinate_based(tags: Dict[str, str]) -> bool:
        """
        Determines if the tags dictionary indicates a coordinate-based system
        by checking for the presence of both X_TAG and Y_TAG.

        Parameters:
            tags (Optional[Dict[str, str]]): A dictionary of tags for the sequencing technology.

        Returns:
            bool: True if both X_TAG and Y_TAG are present, indicating a coordinate-based system.
        """
        return "X_TAG" in tags and "Y_TAG" in tags

    def _process_umi_based(self) -> Dict[str, Dict[str, int]]:
        """
        Processes UMI-based sequencing reads, filtering by quality and matching sequences.

        Returns:
            Dict[str, Dict[str, int]]: A dictionary where keys are cell barcodes, and values are dictionaries with matched sequences as keys and UMI counts as values.
        """
        counts: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        min_read_length = self.min_read_length
        max_read_length = self.max_read_length
        read_length = self.read_length
        pattern_matching = self.pattern_matching
        flanked_sequence = self.flanked_sequence
        bam_file = self.bam_file
        phred33 = self.phred33
        read_qc_percentage = self.read_qc_percentage
        threshold = self.read_qc_threshold
        barcode_path = self.barcode_path
        if barcode_path:
            with gzip.open(barcode_path, "rt") as barcode_file:
                filtered_cell_barcodes = set(barcode_file.read().splitlines())
        for read in bam_file.fetch(contig=self.contig, until_eof=True):
            try:
                cell_bc = (
                    f"{read.get_tag(self.tags['X_TAG'])!s}_{read.get_tag(self.tags['Y_TAG'])!s}"
                    if self.coordinate_based
                    else str(read.get_tag(self.tags["CELL_BARCODE_TAG"]))
                )
                umi = str(read.get_tag(self.tags["UMI_TAG"]))
                sequence = read.query_sequence
                quality = list(read.query_qualities) if read.query_qualities else None
                if sequence is None or quality is None:
                    continue
                if barcode_path and cell_bc not in filtered_cell_barcodes:
                    continue
                matched_read = process_reads(
                    pattern_matching=pattern_matching,
                    seq=sequence,
                    quality=quality,
                    phred33=phred33,
                    threshold=threshold,
                    read_qc_percentage=read_qc_percentage,
                    flanked_sequence=flanked_sequence,
                    min_read_length=min_read_length,
                    max_read_length=max_read_length,
                    read_length=read_length,
                    bam=True,
                )
                if matched_read:
                    counts[cell_bc][matched_read][umi] += 1
            except KeyError:
                continue

        final_counts: Dict[str, Dict[str, int]] = {
            cell_bc: {barcode: len(umis) for barcode, umis in barcodes.items()} for cell_bc, barcodes in counts.items()
        }
        counts.clear()
        return final_counts

    def _process_non_umi_based(self, parse: bool = True) -> Dict[str, Dict[str, int]]:
        """
        Processes non-UMI-based sequencing reads, filtering by quality and matching sequences.

        Returns:
            Dict[str, Dict[str, int]]: A dictionary where keys are cell barcodes, and values are dictionaries with matched sequences as keys and their counts as values.
        """
        counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        min_read_length = self.min_read_length
        max_read_length = self.max_read_length
        read_length = self.read_length
        pattern_matching = self.pattern_matching
        flanked_sequence = self.flanked_sequence
        bam_file = self.bam_file
        read_qc_percentage = self.read_qc_percentage
        threshold = self.read_qc_threshold
        phred33 = self.phred33
        barcode_path = self.barcode_path
        if barcode_path:
            with gzip.open(barcode_path, "rt") as barcode_file:
                filtered_cell_barcodes = set(barcode_file.read().splitlines())
        for read in bam_file.fetch(contig=self.contig if not parse else None, until_eof=True):
            try:
                if self.coordinate_based:
                    cell_bc = f"{read.get_tag(self.tags['X_TAG'])!s}_{read.get_tag(self.tags['Y_TAG'])!s}"
                elif parse:
                    cell_bc = extract_tags_parse(str(read.query_name))
                else:
                    cell_bc = str(read.get_tag(self.tags["CELL_BARCODE_TAG"]))
                sequence = read.query_sequence
                quality = list(read.query_qualities) if read.query_qualities else None
                if sequence is None or quality is None:
                    continue
                if barcode_path and cell_bc not in filtered_cell_barcodes:
                    continue
                matched_read = process_reads(
                    pattern_matching=pattern_matching,
                    seq=sequence,
                    quality=quality,
                    phred33=phred33,
                    threshold=threshold,
                    read_qc_percentage=read_qc_percentage,
                    flanked_sequence=flanked_sequence,
                    min_read_length=min_read_length,
                    max_read_length=max_read_length,
                    read_length=read_length,
                    bam=True,
                )
                if matched_read:
                    counts[cell_bc][matched_read] += 1
            except KeyError:
                continue

        return {k: dict(v) for k, v in counts.items()}

    def extract_barcodes(self) -> Union[Dict[str, Dict[str, int]], Dict[str, Dict[str, Dict[str, int]]]]:
        """
        Extracts barcodes from the BAM file based on the sequencing technology.

        Returns:
            Union[Dict[str, Dict[str, int]], Dict[str, Dict[str, Dict[str, int]]]]: The extracted barcodes, either as a nested dictionary for UMI-based or Parse-based technologies, or a simpler dictionary for non-UMI-based technologies.
        """
        if "UMI_TAG" in BAM_TAGS[self.sequencing_technology]:
            return self._process_umi_based()
        else:
            return self._process_non_umi_based(parse=self.sequencing_technology == "Parse")
