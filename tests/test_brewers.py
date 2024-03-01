from tagbrewer.tag import brewers
from typing import Dict, List
import pytest

@pytest.fixture(scope="module")
def vbrewer() -> brewers.VBrewer:
    chain = "B"
    species = "Homo sapiens"
    return brewers.VBrewer(chain, species)

@pytest.fixture(scope="module")
def jbrewer() -> brewers.JBrewer:
    chain = "B"
    species = "Homo sapiens"
    return brewers.JBrewer(chain, species)

class TestConstructLengthMethods:

    def test_get_max_gene_length(self, vbrewer: brewers.VBrewer) -> None:
        max_gene_length = vbrewer.get_max_gene_length("D")
        assert max_gene_length == len("gggactagcggggggg") # TRBD2

    def test_conservative_v_gene_start_index(self, vbrewer: brewers.VBrewer) -> None:
        start_index = vbrewer.conservative_v_gene_start_index()
        assert start_index == -57

class TestBrewAllTags:

    def test_jbrewer(self, jbrewer: brewers.JBrewer) -> None:
        all_tags = jbrewer.brew_all_tags()
        assert all_tags['TRBJ1-1'][0] == "agctttctttggacaaggca"

    def test_vbrewer(self, vbrewer: brewers.VBrewer) -> None:
        all_tags = vbrewer.brew_all_tags()
        assert all_tags['TRBV1'][0] == "tgtggtcgcactgcagcaag"

class TestBrewTags:

    def test_jbrewer(self, jbrewer: brewers.JBrewer) -> None:
        tags = jbrewer.brew_tags()
        print(tags.gene_layers)
        assert tags.gene_layers['TRBJ1-1'][0][0] == "agctttctttggacaaggca"
        assert len(tags.gene_layers['TRBJ1-1'].keys()) == 1

    def test_vbrewer(self, vbrewer: brewers.VBrewer) -> None:
        tags = vbrewer.brew_tags()
        print(tags.gene_layers)
        assert tags.gene_layers['TRBV1'][0][0] == "tgtggtcgcactgcagcaag"
        assert tags.gene_layers["TRBV7-3"][1] == "Cannot find unique tag set."
        
class TestBrewerLoopMethods:
    
    @pytest.fixture
    def check_genes(self) -> Dict[str, List[str]]:
        return {"gene1": ["abc"],
                "gene2": ["def"],
                "gene3": ["def", "jkl"]}

    def test_create_check_list(self, jbrewer: brewers.JBrewer, check_genes: Dict[str, List[str]]) -> None:
        gene = "gene1"
        check_list = jbrewer.create_check_list(gene, check_genes)
        assert check_list == ["def", "def", "jkl"]

    def test_find_clashes(self, jbrewer: brewers.JBrewer, check_genes: Dict[str, List[str]]) -> None:
        best_tag = "def"
        clashes = jbrewer.find_clashes(check_genes, best_tag)
        assert clashes == ["gene2", "gene3"]

    def test_check_over_alleles(self, jbrewer: brewers.JBrewer) -> None:
        possible_tags = ["abc", "def"]
        alleles = {"01": "abcdef", "02": "defghi"}
        valid_tags = jbrewer.check_over_alleles(possible_tags, alleles)
        assert valid_tags == ["def"]

class TestFindUndecombinable:

    def test_find_undecombinable_j(self, jbrewer: brewers.JBrewer) -> None:
        undecombinable = jbrewer.find_undecombinable()
        assert undecombinable == set()

    def test_find_undecombinable_v(self, vbrewer: brewers.VBrewer) -> None:
        undecombinable = vbrewer.find_undecombinable()
        assert undecombinable == set(('TRBV11-3',
                                      'TRBV7-3',
                                      'TRBV23/OR9-2',
                                      'TRBV6-6',
                                      'TRBV3-2',
                                      'TRBV6-2',
                                      'TRBV12-4',
                                      'TRBV24/OR9-2',
                                      'TRBV3-1',
                                      'TRBV20/OR9-2',
                                      'TRBV6-3',
                                      'TRBV12-3'))
        
class TestDropNonFunctional:

    def test_drop_functional(self):
        chain = "A"
        species = "Homo sapiens"
        brewer = brewers.JBrewer(chain, species, functional=True)
        tags = brewer.brew_tags()
        assert not ("TRAJ51" in tags.gene_layers.keys())
        assert "TRAJ10" in tags.gene_layers.keys()
