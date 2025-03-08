from typing import Dict, List, Optional, Union

import edlib
import flpc
from ahocorasick_rs import AhoCorasick, MatchKind
from cutadapt.adapters import BackAdapter, FrontAdapter, LinkedAdapter, Match
from cutadapt.modifiers import AdapterCutter
from dnaio import Sequence


class PatternMatching:
    """
    Utility class for sequence operations, supporting regex and Aho-Corasick pattern matching.

    Attributes:
        flanked_pattern (bool): Whether to check for flanked sequences when using regex.
        regex_pattern (Optional[flpc.Pattern]): The compiled regex pattern for sequence matching.
        sequence_pattern (Optional[Union[List[str], Dict[str, str]]]): The list or dictionary for Aho-Corasick sequence matching.
        ac (Optional[AhoCorasick]): The Aho-Corasick automaton for list or dictionary pattern matching.
        pattern_type (Optional[str]): The type of pattern ('N', 'exact', 'left', 'right', 'both').
        read_length (Optional[int]): The length of the read sequence to be extracted when pattern is 'flanking*' or '*flanking'.
        max_read_length (Optional[int]): The maximum length of the read sequence to be extracted.
        aln_mismatches (Optional[int]): The maximum number of mismatches allowed in the alignment. If not None, edlib is used, otherwise Aho-Corasick is used.
        cutadapt_adapter (Optional[Union[BackAdapter, FrontAdapter, LinkedAdapter]]): The initialized Cutadapt adapter.
        adapter_cutter (Optional[AdapterCutter]): The initialized Cutadapt adapter cutter.
        left_flanking_coverage (Optional[int]): The minimum bp overlap of the left flanking region for cutadapt.
        right_flanking_coverage (Optional[int]): The minimum bp overlap of the right flanking region for cutadapt.
        references (Union[List[str], Dict[str, str]]): The list or dictionary of reference sequences for alignment.
    """

    def __init__(
        self,
        pattern: Union[str, List[str], Dict[str, str]],
        flanked_pattern: bool = False,
        read_length: Optional[int] = None,
        max_read_length: Optional[int] = None,
        aln_mismatches: Optional[int] = None,
        flanking_mismatches: Optional[float] = None,
        left_flanking_coverage: Optional[int] = None,
        right_flanking_coverage: Optional[int] = None,
    ) -> None:
        """
        Initialize PatternMatching with a pattern and optional flanked sequence detection.

        Args:
            pattern (Union[str, List[str], Dict[str, str]]): The pattern for sequence matching. It can be a regex pattern,
                a list of sequences, or a dictionary with sequences as keys and associated values.
            flanked_pattern (bool, optional): Whether to check for flanked sequences when using regex patterns. Defaults to False.
            read_length (Optional[int], optional): The length of the read sequence to be extracted when pattern is 'flanking*' or '*flanking'.
            max_read_length (Optional[int], optional): The maximum length of the read sequence to be extracted.
            aln_mismatches (Optional[int]): The maximum number of mismatches allowed in the alignment. If not None, edlib is used, otherwise Aho-Corasick is used.
            flanking_mismatches (Optional[float]): The number of mismatches allowed in the flanking regions. If None, regex matching is used.
            left_flanking_coverage (Optional[int]): The minimum bp overlap of the left flanking region for cutadapt.
            right_flanking_coverage (Optional[int]): The minimum bp overlap of the right flanking region for cutadapt.

        Raises:
            TypeError: If the pattern is not a string, list, or dictionary.
            ValueError: If the pattern contains both 'N' and '*', or invalid '*' usage.
        """
        self.flanked_pattern: bool = flanked_pattern
        self.read_length: Optional[int] = read_length
        self.max_read_length: Optional[int] = max_read_length
        self.regex_pattern: Optional[flpc.Pattern] = None
        self.sequence_pattern: Optional[Union[List[str], Dict[str, str]]] = None
        self.ac: Optional[AhoCorasick] = None
        self.pattern_type: Optional[str] = None
        self.flanking_mismatches: Optional[float] = flanking_mismatches
        self.aln_mismatches: Optional[int] = aln_mismatches
        self.cutadapt_adapter: Optional[Union[BackAdapter, FrontAdapter, LinkedAdapter]] = None
        self.adapter_cutter: Optional[AdapterCutter] = None
        self.left_flanking_coverage: Optional[int] = left_flanking_coverage
        self.right_flanking_coverage: Optional[int] = right_flanking_coverage
        self.references_list: Optional[List[str]] = None
        self.references_dict: Optional[Dict[str, str]] = None
        self.find_sequence = self._default_find_sequence
        if isinstance(pattern, str):
            self._process_string_pattern(pattern)
            if self.cutadapt_adapter and self.flanking_mismatches:
                self.find_sequence = self._find_sequence_with_cutadapt
            if self.regex_pattern and self.pattern_type in {"N", "exact", "both", "left", "right"}:
                self.find_sequence = self._find_sequence_from_pattern
        elif isinstance(pattern, dict):
            self.sequence_pattern = pattern
            if not aln_mismatches:
                self.ac = AhoCorasick(list(self.sequence_pattern.keys()), matchkind=MatchKind.LeftmostLongest)
                self.find_sequence = self._find_sequence_from_dict_aho_corasick
            else:
                self.references_dict = pattern
                self.find_sequence = self._find_sequence_from_dict_edlib
        elif isinstance(pattern, list):
            self.sequence_pattern = pattern
            if not aln_mismatches:
                self.ac = AhoCorasick(self.sequence_pattern, matchkind=MatchKind.LeftmostLongest)
                self.find_sequence = self._find_sequence_from_list_aho_corasick
            else:
                self.references_list = pattern
                self.find_sequence = self._find_sequence_from_list_edlib
        else:
            e = "Pattern must be a string, list, or dictionary"
            raise TypeError(e)

    def _process_string_pattern(self, pattern: str) -> None:
        """
        Process the string pattern to determine the pattern type and flanking sequences.

        Args:
            pattern (str): The input pattern containing 'N' or '*'.

        Raises:
            ValueError: If the pattern contains both 'N' and '*', or invalid '*' usage.
        """
        if "N" in pattern:
            regex_str = self._convert_n_to_regex(pattern)
            self.regex_pattern = flpc.compile(regex_str)
            self.pattern_type = "N"
        elif "*" in pattern and not self.flanking_mismatches:  # regex
            if pattern.startswith("*"):
                # Pattern is '*flanking'
                flanking = pattern[1:]
                self.pattern_type = "left"
                self.regex_pattern = flpc.compile(
                    f"(.{{{self.read_length}}}){flanking}" if self.read_length else f"(.*){flanking}"
                )
            elif pattern.endswith("*"):
                # Pattern is 'flanking*'
                flanking = pattern[:-1]
                self.pattern_type = "right"
                self.regex_pattern = flpc.compile(
                    f"{flanking}(.{{{self.read_length}}})" if self.read_length else f"{flanking}(.*)"
                )
            else:
                # Pattern is 'flanking*flanking'
                self.pattern_type = "both"
                self.regex_pattern = flpc.compile(pattern.replace("*", "(.*)"))
        elif self.flanking_mismatches:  # cutadapt
            self._initialize_cutadapt_adapters(pattern)
        else:
            self.regex_pattern = flpc.compile(pattern)
            self.pattern_type = "exact"

    def _convert_n_to_regex(self, pattern: str) -> str:
        """
        Convert 'N' in the pattern to regex equivalents that match sequences of fixed length.

        Args:
            pattern (str): The input pattern containing 'N'.

        Returns:
            str: The regex string with 'N' converted to match sequences of fixed length.
        """
        # Convert sequences of 'N's into .{count}
        parts = []
        count = 0
        for char in pattern:
            if char == "N":
                count += 1
            else:
                if count > 0:
                    parts.append(f".{{{count}}}")
                    count = 0
                parts.append(char)
        if count > 0:
            parts.append(f".{{{count}}}")
        return "".join(parts)

    def _initialize_cutadapt_adapters(self, pattern: str) -> None:
        """
        Initialize Cutadapt adapters based on the pattern structure and flanking_mismatches.

        Args:
            pattern (str): The input pattern containing wildcards (*).

        Raises:
            ValueError: If the pattern starts and ends with '*', or if the pattern format is invalid for flanking_mismatches.

        Behavior:
            - If the pattern starts with '*', a BackAdapter is initialized.
            - If the pattern ends with '*', a FrontAdapter is initialized.
            - If the pattern contains a single '*', a LinkedAdapter is initialized.
            - The adapters share the `flanking_mismatches` value for controlling mismatches.
        """
        if pattern.startswith("*"):
            flanking = pattern[1:]
            self.pattern_type = "left"
            self.cutadapt_adapter = BackAdapter(
                sequence=flanking,
                max_errors=self.flanking_mismatches,
                min_overlap=len(flanking) if not self.right_flanking_coverage else self.right_flanking_coverage,
            )
        elif pattern.endswith("*"):
            flanking = pattern[:-1]
            self.pattern_type = "right"
            self.cutadapt_adapter = FrontAdapter(
                sequence=flanking,
                max_errors=self.flanking_mismatches,
                min_overlap=len(flanking) if not self.left_flanking_coverage else self.left_flanking_coverage,
            )
        else:
            parts = pattern.split("*")
            if len(parts) != 2:
                e = "Invalid pattern format for cutadapt."
                raise ValueError(e)
            front_adapter = FrontAdapter(
                sequence=parts[0],
                max_errors=self.flanking_mismatches,
                min_overlap=len(parts[0]) if not self.left_flanking_coverage else self.left_flanking_coverage,
            )
            back_adapter = BackAdapter(
                sequence=parts[1],
                max_errors=self.flanking_mismatches,
                min_overlap=len(parts[1]) if not self.right_flanking_coverage else self.right_flanking_coverage,
            )
            self.pattern_type = "both"
            self.cutadapt_adapter = LinkedAdapter(
                front_adapter=front_adapter,
                back_adapter=back_adapter,
                front_required=True,
                back_required=True,
                name="linked_adapter",
            )
        self.adapter_cutter = AdapterCutter([self.cutadapt_adapter])

    def _find_sequence_from_dict_aho_corasick(self, read: str) -> Optional[str]:
        """
        Find a valid sequence in a read using the provided reference dictionary.

        Parameters:
            read (str): The input sequence read.

        Returns:
            Optional[str]: The matched sequence if found, otherwise None.
        """
        if "N" in read or not self.ac or not isinstance(self.sequence_pattern, dict):
            return None
        matches = self.ac.find_matches_as_strings(read)
        if matches:
            return (
                self.sequence_pattern[matches[0]] if len(matches) == 1 else self.sequence_pattern[max(matches, key=len)]
            )
        return None

    def _find_sequence_from_list_aho_corasick(self, read: str) -> Optional[str]:
        """
        Find a valid sequence in a read using the provided reference list.

        Parameters:
            read (str): The input sequence read.

        Returns:
            Optional[str]: The matched sequence if found, otherwise None.
        """
        if "N" in read or not self.ac or not isinstance(self.sequence_pattern, list):
            return None
        matches = self.ac.find_matches_as_strings(read)
        if matches:
            return matches[0] if len(matches) == 1 else max(matches, key=len)
        return None

    def _find_sequence_from_dict_edlib(self, read: str) -> Optional[str]:
        """
        Aligns a query sequence to a list or dictionary of references and stops when the best alignment is found.

        Parameters:
            read (str): The sequence to align.

        Returns:
            (Optional[str]): Aligned reference with the best match.
        """
        best_alignment = None
        if not self.references_dict:
            return None
        for sequence, value in self.references_dict.items():
            result = edlib.align(sequence, read, mode="HW", task="distance", k=self.aln_mismatches)
            edit_distance = result["editDistance"]
            print(result)
            if edit_distance != -1 and (best_alignment is None or edit_distance < best_alignment["edit_distance"]):
                best_alignment = {
                    "aligned_sequence": sequence,
                    "name": value,
                    "edit_distance": edit_distance,
                }
                if edit_distance == 0:
                    break

        return best_alignment["name"] if best_alignment else None

    def _find_sequence_from_list_edlib(self, read: str) -> Optional[str]:
        """
        Aligns a query sequence to a list or dictionary of references and stops when the best alignment is found.

        Parameters:
            read (str): The sequence to align.

        Returns:
            (Optional[str]): Aligned reference with the best match.
        """
        best_alignment = None
        if not self.references_list:
            return None
        for sequence in self.references_list:
            result = edlib.align(sequence, read, mode="HW", task="distance", k=self.aln_mismatches)
            edit_distance = result["editDistance"]
            if edit_distance != -1 and (best_alignment is None or edit_distance < best_alignment["edit_distance"]):
                best_alignment = {
                    "aligned_sequence": sequence,
                    "edit_distance": edit_distance,
                }
                if edit_distance == 0:
                    break

        return best_alignment["aligned_sequence"] if best_alignment else None

    def _find_sequence_from_pattern(self, read: str) -> Optional[str]:
        """
        Find a valid sequence in a read using the compiled regex pattern based on the pattern type.

        Parameters:
            read (str): The input sequence read.

        Returns:
            Optional[str]: The matched sequence if found, otherwise None.
        """
        if not self.regex_pattern:
            e = "No regex pattern found"
            raise ValueError(e)

        if self.pattern_type == "N" or self.pattern_type == "exact":
            _match = flpc.search(self.regex_pattern, read)
            if _match:
                match_group = _match.group(0)
                if match_group and "N" not in match_group:
                    return match_group

        else:
            _match = flpc.search(self.regex_pattern, read)
            if _match:
                match_group = _match.group(1)
                if match_group and "N" not in match_group:
                    if self.max_read_length:
                        # trim the sequence to the max_read_length
                        if self.pattern_type == "left" and len(match_group) > self.max_read_length:
                            return match_group[-self.max_read_length :]
                        elif self.pattern_type == "right" and len(match_group) > self.max_read_length:
                            return match_group[: self.max_read_length]
                    return match_group

        return None

    def _find_sequence_with_cutadapt(self, read: str) -> Optional[str]:
        """
        Use Cutadapt to find and extract the matching sequence from the input.

        Args:
            read (str): The input sequence to search within.

        Returns:
            Optional[str]: The matched sequence if found, otherwise None.
        """
        if self.adapter_cutter is None:
            return None
        _read = Sequence(name="read_placeholder", sequence=read)
        info = TrimInfo()
        trimmed_read = self.adapter_cutter(_read, info=info)
        if len(trimmed_read.sequence) != len(read):
            _match: str = trimmed_read.sequence
            if self.max_read_length:
                # trim the sequence to the max_read_length
                if self.pattern_type == "left" and len(_match) > self.max_read_length:
                    return _match[-self.max_read_length :]
                elif self.pattern_type == "right" and len(_match) > self.max_read_length:
                    return _match[: self.max_read_length]
            elif self.read_length and len(_match) != self.read_length:
                return None
            return _match
        return None

    def _default_find_sequence(self, read: str) -> Optional[str]:
        """Default implementation of find_sequence."""
        raise NotImplementedError("No find_sequence implementation assigned.")


class TrimInfo:
    def __init__(self) -> None:
        self.matches: List[Match] = []
