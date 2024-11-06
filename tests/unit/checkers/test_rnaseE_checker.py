import pytest
from genedesign.checkers.rnaseE_checker import RNaseE_Checker

@pytest.fixture
def rnase_checker():
    """Fixture to initialize the RNaseE_Checker."""
    checker = RNaseE_Checker()
    checker.initiate()
    return checker

def test_no_cleavage_site(rnase_checker):
    """Test case where no cleavage site is present."""
    rna_sequence = "AUGGCCAUUGGCGCUAG"  # No RNase E motifs
    result, motif = rnase_checker.run(rna_sequence)
    assert result is True
    assert motif is None

def test_single_cleavage_site(rnase_checker):
    """Test case where a single cleavage site is present."""
    rna_sequence = "AUGGCCAUUUAUGCGCUAG"  # Contains "AUUUA"
    result, motif = rnase_checker.run(rna_sequence)
    assert result is False
    assert motif == "AUUUA"

def test_multiple_cleavage_sites(rnase_checker):
    """Test case where multiple cleavage sites are present."""
    rna_sequence = "AAUUAAUUUAUUUAA"  # Contains multiple motifs: "AUUUA", "AAUU", "UUAA"
    result, motif = rnase_checker.run(rna_sequence)
    assert result is False
    assert motif in ["AUUUA", "AAUU", "UUAA"]  # Any of these could be detected first

def test_case_insensitivity(rnase_checker):
    """Test case to ensure case insensitivity."""
    rna_sequence = "auggccauuuaugcgcuag"  # Lowercase version of a sequence with "AUUUA"
    result, motif = rnase_checker.run(rna_sequence)
    assert result is False
    assert motif == "AUUUA"
