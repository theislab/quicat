from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple, Union

import numpy as np
from joblib import Parallel, delayed
from rapidfuzz import distance, process
from scipy.sparse import csr_matrix
from sklearn.feature_extraction.text import CountVectorizer


class ForestBuilder:
    """
    Builds a forest structure to group similar sequences based on a distance metric.

    Attributes:
        root_nodes (defaultdict[str, List[str]]): The root nodes of the forest, grouping similar sequences.
        scorer (Callable): The distance metric function used for grouping sequences.
    """

    def __init__(self, distance_metric: str = "Levenshtein") -> None:
        """
        Initializes ForestBuilder with the specified distance metric.

        Parameters:
            distance_metric (str): The distance metric to use ('Levenshtein' or 'Hamming').
        """
        self.root_nodes: defaultdict[str, List[str]] = defaultdict()
        self.scorer = distance.Levenshtein.distance if distance_metric == "Levenshtein" else distance.Hamming.distance

    @staticmethod
    def _add_to_forest(root_nodes: defaultdict[str, List[str]], sequence: str, threshold: int, scorer: Any) -> None:
        """
        Adds a sequence to the forest structure based on the distance threshold.

        Parameters:
            root_nodes (defaultdict[str, List[str]]): The root nodes of the forest.
            sequence (str): The sequence to add.
            threshold (int): The distance threshold for grouping.
            scorer (Callable): The distance metric function.
        """
        match = process.extractOne(sequence, list(root_nodes.keys()), scorer=scorer, score_cutoff=threshold)
        if match:
            root_nodes[match[0]].append(sequence)
        else:
            root_nodes[sequence] = []

    def _build_forest(self, sequences: List[str], threshold: int) -> defaultdict[str, List[str]]:
        """
        Builds the forest structure grouping sequences based on the distance threshold.

        Parameters:
            sequences (List[str]): The list of sequences to group.
            threshold (int): The distance threshold for grouping.

        Returns:
            defaultdict[str, List[str]]: The forest structure grouping similar sequences.
        """
        if not sequences:
            return defaultdict(list)
        add_to_forest = self._add_to_forest
        scorer = self.scorer
        root_nodes = self.root_nodes
        root_nodes[sequences[0]] = []

        for sequence in sequences[1:]:
            add_to_forest(root_nodes, sequence, threshold, scorer)

        return root_nodes

    def extract_mapping(self, sequences: List[str], threshold: int) -> List[Set[str]]:
        """
        Extracts groups of similar sequences from the forest structure.

        Parameters:
            sequences (List[str]): The list of sequences to group.
            threshold (int): The distance threshold for grouping.

        Returns:
            List[Set[str]]: A list of sets, each containing similar sequences.
        """
        forest = self._build_forest(sequences, threshold)
        return [{root, *children} for root, children in forest.items()]


class SequencesCollapser:
    """
    Collapses similar sequences into clusters based on error thresholds.

    Attributes:
        sequences (List[str]): The list of sequences to collapse.
        error_threshold (int): The allowable error threshold for collapsing sequences.
        n_threads (int): Number of threads to use for parallel processing.
        scorer (str): The distance metric used for scoring (Levenshtein or Hamming).
    """

    def __init__(self, sequences: List[str], error_threshold: int = 2, n_threads: int = -1) -> None:
        """
        Initializes SequencesCollapser with sequences, error threshold, and thread count.

        Parameters:
            sequences (List[str]): The list of sequences to collapse.
            error_threshold (int): The allowable error threshold for collapsing sequences.
            n_threads (int): Number of threads to use for parallel processing.
        """
        self.sequences = sequences
        self.error_threshold = error_threshold
        self.n_threads = n_threads

    def _build_ngrams(self) -> csr_matrix:
        """
        Builds an n-gram representation of the sequences for clustering.

        Returns:
            csr_matrix: A sparse matrix of n-grams.
        """
        n, scorer = self._calculate_n_value()
        self.scorer = scorer
        vectorizer = CountVectorizer(analyzer="char", ngram_range=(n, n), lowercase=False, dtype=np.float32)
        matrix = vectorizer.fit_transform(self.sequences)
        return matrix

    def _group_sequences(self) -> List[List[str]]:
        """
        Groups similar sequences based on n-gram similarities.

        Returns:
            List[Set[str]]: A list of sets, each containing similar sequences.
        """
        groups = []
        matrix = self._build_ngrams()
        sequences = self.sequences.copy()

        while True:
            previous_length = len(sequences)
            for i, seq in enumerate(sequences):
                dot_product = matrix[i].dot(matrix.T)
                indices = dot_product.indices
                group = [seq] + [sequences[idx] for idx in indices if sequences[idx] != seq]
                groups.append(group)
                idx = list(set(range(len(sequences))) - {*indices})
                matrix = matrix[idx]
                sequences = [sequences[j] for j in idx]
                break
            if len(sequences) == previous_length:
                break
        return groups

    def cluster(self) -> Dict[str, int]:
        """
        Clusters sequences into groups based on error threshold.

        Returns:
            Dict[str, int]: A dictionary mapping each sequence to its cluster ID.
        """
        initial_groups = self._group_sequences()

        scorer = self.scorer
        error_threshold = self.error_threshold
        final_groups = []
        groups_to_process = [group for group in initial_groups if len(group) > 1]
        single_element_groups = [group for group in initial_groups if len(group) == 1]
        print(f"Initial groups: {len(initial_groups)}")
        print(f"Groups to process: {len(groups_to_process)}")
        if groups_to_process:
            parallel = Parallel(n_jobs=self.n_threads, backend="loky")
            refined_groups_results = parallel(
                delayed(self._process_group)(group, error_threshold, scorer) for group in groups_to_process
            )

            for result in refined_groups_results:
                final_groups.extend(result)
        final_groups.extend(single_element_groups)
        cluster_dict = {}
        for cluster_id, group in enumerate(final_groups):
            for sequence in group:
                cluster_dict[sequence] = cluster_id
        return cluster_dict

    def _calculate_n_value(self) -> Tuple[int, str]:
        """
        Calculates the optimal n value for n-gram creation based on sequence lengths.

        Returns:
            Tuple[int, str]: The n value and the appropriate scorer ('Hamming' or 'Levenshtein').
        """
        lengths = [len(seq) for seq in self.sequences]
        L = lengths[0] if len(set(lengths)) == 1 else max(lengths)
        scorer = "Hamming" if len(set(lengths)) == 1 else "Levenshtein"
        n = L // (self.error_threshold + 1)
        return max(3, n), scorer

    @staticmethod
    def _process_group(group: List[str], error_threshold: int, scorer: str) -> Union[Set[str], List[Set[str]]]:
        """
        Processes a group of sequences to refine clusters.

        Parameters:
            group (Set[str]): The group of sequences to process.
            error_threshold (int): The allowable error threshold for clustering.
            scorer (str): The distance metric to use for scoring.

        Returns:
            List[Set[str]]: A list of refined clusters.
        """
        if len(group) == 1:
            return set(group)
        forest_builder = ForestBuilder(scorer)
        sequences = group
        refined_groups = forest_builder.extract_mapping(sequences, error_threshold)
        return refined_groups
