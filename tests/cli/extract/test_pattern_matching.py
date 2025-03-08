from quicat.cli.extract import PatternMatching


# Test cases for exact regex pattern matching
def test_regex_exact_pattern():
    pm = PatternMatching(pattern="ACGT")
    result = pm.find_sequence("ACGTACGT")
    assert result == "ACGT", "Should find exact match for 'ACGT'"


# Test cases for 'N' pattern
def test_regex_N_pattern():
    pm = PatternMatching(pattern="ACGTN")
    result = pm.find_sequence("ACGTGACGT")
    assert result == "ACGTG", "Should match 'ACGTG' for 'ACGTN' pattern"


# Test case for left-flanking pattern '*flanking'
def test_left_flanking_pattern():
    pm = PatternMatching(pattern="*CGT", flanked_pattern=True)
    result = pm.find_sequence("ACGAACGT")
    assert result == "ACGAA", "Should find sequence matching '*CGT' with left flanking"


# Test case for right-flanking pattern 'flanking*'
def test_right_flanking_pattern():
    pm = PatternMatching(pattern="ACGT*", flanked_pattern=True)
    result = pm.find_sequence("ACGTTCGT")
    assert result == "TCGT", "Should find sequence matching 'ACGT*' with right flanking"


# Test case for both-flanking pattern 'flanking*flanking'
def test_both_flanking_pattern():
    pm = PatternMatching(pattern="A*GT", flanked_pattern=True)
    result = pm.find_sequence("ACGTACGT")
    assert result == "CGTAC", "Should find sequence matching 'A*GT' with both flanking"


# Test case for dictionary pattern matching
def test_dict_pattern_match():
    pm = PatternMatching(pattern={"ACGT": "MATCH1", "GCTA": "MATCH2"})
    result = pm.find_sequence("ACGTGCTAGCAT")
    assert result == "MATCH1", "Should match the longest dictionary key and return the associated value"


# Test case for list pattern matching
def test_list_pattern_match():
    pm = PatternMatching(pattern=["ACGT", "GCTA", "CGTA"])
    result = pm.find_sequence("CGTAGCTAGCAT")
    assert result == "CGTA", "Should match the longest sequence in the list"


# Test case for dictionary pattern matching
def test_dict_pattern_match_with_edlilib():
    pm = PatternMatching(pattern={"ACGC": "MATCH1", "GCTAAA": "MATCH2"}, aln_mismatches=1)
    result = pm.find_sequence("ACGTGCTAGCAT")
    assert result == "MATCH1"
    result = pm.find_sequence("TTTTTTTTTTTTTT")
    assert result is None, "Should not find any match here"


# Test case for list pattern matching
def test_list_pattern_match_with_mismatch():
    pm = PatternMatching(pattern=["ACGC", "GCTAAA"], aln_mismatches=1)
    result = pm.find_sequence("ACGTGCTAGCAT")
    assert result == "ACGC"
    result = pm.find_sequence("TTTTTTTTTTTTTT")
    assert result is None, "Should not find any match here"


# Test case where no match is found
def test_no_match_for_invalid_sequence():
    pm = PatternMatching(pattern="ACGT")
    result = pm.find_sequence("TTTTTTTT")
    assert result is None, "Should return None when no match is found"


# Test case for variable length 'N' matching
def test_N_pattern_with_different_length():
    pm = PatternMatching(pattern="ACGTNN")
    result = pm.find_sequence("ACGTGG")
    assert result == "ACGTGG", "Should match two 'N' with 'GG'"


# Test case for leftmost longest match in dictionary
def test_leftmost_longest_match():
    pm = PatternMatching(pattern={"ACGT": "MATCH1", "ACG": "MATCH2"})
    result = pm.find_sequence("ACGTGCTAGCAT")
    assert result == "MATCH1", "Should return 'MATCH1' because it matches the longest pattern"


# Test case for flanked pattern matching
def test_match_with_flanked_pattern():
    pm = PatternMatching(pattern="AA*CGT", flanked_pattern=True)
    result = pm.find_sequence("AAACTTGCTACGT")
    assert result == "ACTTGCTA", "Should find the match even if flanked_pattern is True"


# Test case for Aho-Corasick automaton with invalid dictionary pattern
def test_aho_corasick_invalid_pattern():
    pm = PatternMatching(pattern={"ACGT": "MATCH1", "GCTN": "MATCH2"})
    result = pm.find_sequence("ACGTGCTAGCAT")
    assert result == "MATCH1", "Should ignore invalid patterns containing 'N'"


def test_cutadapt_front_adapter_without_mismatches():
    pm = PatternMatching(pattern="ACGT*", flanking_mismatches=0.1)
    result = pm.find_sequence("ACGTTCGT")
    assert result == "TCGT", "Should find 'TCGT' even with small mismatches for front adapter"


def test_cutadapt_front_adapter_with_mismatches():
    pm = PatternMatching(pattern="ACGTACGTAT*", flanking_mismatches=0.1)
    result = pm.find_sequence("ACGTCCGTATTCGT")
    assert result == "TCGT", "Should find 'TCGT' even with small mismatches for front adapter"


# Test case for Cutadapt back adapter with mismatches
def test_cutadapt_back_adapter_without_mismatches():
    pm = PatternMatching(pattern="*GCTAGCTATT", flanking_mismatches=0.1)
    result = pm.find_sequence("TTTTGCTAGCTATA")
    assert result == "TTTT", "Should find 'TTTT' even with small mismatches for back adapter"


# Test case for Cutadapt linked adapter with mismatches
def test_cutadapt_linked_adapter_without_mismatches():
    pm = PatternMatching(pattern="ACGT*GCTA", flanking_mismatches=0.1)
    result = pm.find_sequence("ACGTATTAGCTA")
    assert result == "ATTA", "Should find 'ATTA' matching linked adapter with small mismatches"


def test_cutadapt_linked_adapter_with_mismatches():
    pm = PatternMatching(pattern="ACACGT*GCTAAT", flanking_mismatches=0.1)
    result = pm.find_sequence("ACACGTATTAGCTAAT")
    assert result == "ATTA", "Should find 'ATTA' matching linked adapter with small mismatches"
