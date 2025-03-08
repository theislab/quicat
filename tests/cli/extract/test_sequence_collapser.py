from collections import defaultdict

from quicat.cli.extract import SequencesCollapser


def test_sequencescollapser_identical_sequences_hamming():
    sequences = ["ACGT", "ACGG", "ACGA"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=1)
    result = collapser.cluster()
    assert len(result) == len(sequences), "Each sequence should have its own cluster with a low error threshold"


def test_sequencescollapser_varying_lengths_levenshtein():
    sequences = ["ACGT", "ACGTTT", "ACGTT"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=2)
    result = collapser.cluster()
    cluster_values = set(result.values())
    assert len(cluster_values) == 1, "All sequences should be in the same cluster with Levenshtein distance"


def test_sequencescollapser_multithreading():
    sequences = ["ACGT", "ACGTT", "ACGTA", "ACGTC"]
    single_threaded = SequencesCollapser(sequences=sequences, error_threshold=1, n_threads=1).cluster()
    multi_threaded = SequencesCollapser(sequences=sequences, error_threshold=1, n_threads=2).cluster()
    assert single_threaded == multi_threaded, "Single-threaded and multi-threaded results should be consistent"


def test_sequencescollapser_large_error_threshold():
    sequences = ["ATGACGTAACGT", "ATGACGTAATTT", "ATGACGTAAGGG", "ATGACGTAACCC"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=3)
    result = collapser.cluster()
    assert len(set(result.values())) == 1, "With a large error threshold, all sequences should be in the same cluster"


def test_sequencescollapser_grouping_sequences():
    sequences = ["ACTG", "ACTA", "ACGG", "GCTG", "GCTC", "GTTG"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=1)
    result = collapser.cluster()
    groups = defaultdict(list)
    for seq, group_id in result.items():
        groups[group_id].append(seq)

    assert len(groups) > 1, "Sequences should be grouped into multiple clusters"
    assert any(
        len(group) > 1 for group in groups.values()
    ), "There should be at least one group with more than one sequence"


# Test specific scorer (Hamming and Levenshtein)
def test_sequencescollapser_hamming_scorer():
    sequences = ["AAAA", "AAAT", "AAAG", "AATT"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=1)
    result = collapser.cluster()
    groups = defaultdict(list)
    for seq, group_id in result.items():
        groups[group_id].append(seq)

    assert len(groups) > 1, "Sequences should be grouped into multiple clusters with Hamming distance"
    assert any(len(group) > 1 for group in groups.values()), "At least one group should have more than one sequence"


def test_sequencescollapser_levenshtein_scorer():
    sequences = ["ATGACGTACTG", "ATGACGTACTAA", "ATGACGTACGG", "ATGACGTGCTG"]
    collapser = SequencesCollapser(sequences=sequences, error_threshold=3)
    result = collapser.cluster()
    groups = defaultdict(list)
    for seq, group_id in result.items():
        groups[group_id].append(seq)

    assert len(groups) == 1, "All sequences should be grouped into one cluster with Levenshtein distance"
