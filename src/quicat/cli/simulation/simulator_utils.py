import random
from typing import List

from rapidfuzz import distance, process

ALPHABET: str = "ATCG"
PHRED: str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"


def generate_unique_sequences(total_sequences: int, sequences_length: int) -> List[str]:
    unique_umis: set[str] = set()
    while len(unique_umis) < total_sequences:
        umi = "".join(random.choices(ALPHABET, k=sequences_length))
        unique_umis.add(umi)
    return list(unique_umis)


def generate_random_barcodes_with_threshold(num_barcodes: int, barcode_length: int, min_hamming_dist: int) -> List[str]:
    """
    Generates a list of random DNA barcodes ensuring that each barcode is above a specified Hamming distance
    threshold from all others using rapidfuzz's cdist with a score cutoff.

    Parameters:
        num_barcodes (int): Number of barcodes to generate.
        barcode_length (int): Length of each barcode.
        min_hamming_dist (int): Minimum Hamming distance between any two barcodes.

    Returns:
        List[str]: List of random barcode strings that meet the Hamming distance criteria.
    """
    barcodes: List[str] = []

    while len(barcodes) < num_barcodes:
        candidate: str = "".join(random.choices(ALPHABET, k=barcode_length))
        if not barcodes:
            barcodes.append(candidate)
        else:
            distances = process.extract(
                [candidate], barcodes, scorer=distance.Hamming.distance, score_cutoff=min_hamming_dist
            )
            if len(distances) == 0:
                barcodes.append(candidate)
    return barcodes


def generate_pcr_replicates(
    sequence: str, num_replicates: int, max_hamming_dist: int = 2, error_rate: float = 0.1
) -> List[str]:
    """
    Generates PCR error-prone replicates of a given barcode.

    Parameters:
        sequence (str): Original DNA barcode sequence.
        num_replicates (int): Number of PCR replicates to generate.
        max_hamming_dist (int): Maximum Hamming distance for the error-prone replicates. Default is 2.
        error_rate (float): Probability of a mutation at each base. Default is 0.1.

    Returns:
        List[str]: List of PCR replicates with errors.
    """
    replicates: List[str] = []
    seq_len = len(sequence)

    while len(replicates) < num_replicates:
        num_mutations = max_hamming_dist
        mutation_positions = random.sample(range(seq_len), num_mutations)
        seq_list: List[str] = list(sequence)
        for pos in mutation_positions:
            seq_list[pos] = random.choice([x for x in ALPHABET if x != seq_list[pos]])
        mutated_seq: str = "".join(seq_list)
        if distance.Hamming.distance(sequence, mutated_seq) <= max_hamming_dist:
            replicates.append(mutated_seq)
    return replicates


def generate_quality_string(length: int) -> str:
    """
    Generates a random quality string of a given length using PHRED characters.

    Parameters:
        length (int): The length of the quality string.

    Returns:
        str: A randomly generated quality string.
    """
    return "".join(random.choices(PHRED, k=length))
