from tagbrewer.tag import checkers, generators
import pytest
import pathlib

def test_collapse_alleles():
    alleles_dict = {"gene1": {"01": "abc", "02": "aab"},
                    "gene2": {"01": "def", "02": "ddf"}}
    assert checkers.collapse_alleles(alleles_dict) == ["gene1", "gene2"]
    
class TestGeneDifferences:
    @pytest.fixture
    def resource_location(self):
        return pathlib.Path("tests/resources").resolve()

    def test_extract_gene_list_version_parity(self, resource_location):
        original = checkers.extract_gene_list(resource_location, "human", "original", "TRAV")
        extended = checkers.extract_gene_list(resource_location, "human", "extended", "TRAV")
        assert original == extended

    def test_extract_gene_list_extended(self, resource_location):
        extended = checkers.extract_gene_list(resource_location, "human", "extended", "TRAV")
        assert extended == set(["TRAV1-1", "TRAV14/DV4"])

    def test_find_gene_differences(self, resource_location, monkeypatch):
        reference_gene_list = {"TRAV1-1": {"01": "abc", "02": "aab"},
                               "TRAV14/DV4": {"01": "def", "02": "ddf"},
                               "TRAV-2": {"01": "ghi", "02": "ggi"}}
        monkeypatch.setattr(generators, "get_tr_alleles_for_gene_group_for_species", lambda x, y: (None, reference_gene_list))
        differences = checkers.find_gene_differences(resource_location, "human", "extended", "TRAV")
        assert differences == set(["TRAV-2"])

