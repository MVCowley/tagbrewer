from tagbrewer.tag import brewers

class TestMaxGeneLengths:

    def test_get_max_gene_length(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.VBrewer(chain, species)
        max_gene_length = brewer.get_max_gene_length("D")
        assert max_gene_length == len("gggactagcggggggg") # TRBD2

class TestVStartIndex:

    def test_conservative_v_gene_start_index(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.VBrewer(chain, species)
        start_index = brewer.conservative_v_gene_start_index()
        assert start_index == -57

class TestBrewAllTags:

    def test_jbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.JBrewer(chain, species)
        all_tags = brewer.brew_all_tags()
        assert all_tags['TRBJ1-1'][0] == "agctttctttggacaaggca"

    def test_vbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.VBrewer(chain, species)
        all_tags = brewer.brew_all_tags()
        assert all_tags['TRBV1'][0] == "tgtggtcgcactgcagcaag"

class TestBrewTags:

    # TODO: alter test to specifically test for a tag which overlaps
    # perhaps negatative assertion

    def test_jbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.JBrewer(chain, species)
        tags = brewer.brew_tags()
        assert tags['TRBJ1-1'][0] == "agctttctttggacaaggca"

    def test_vbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.VBrewer(chain, species)
        tags = brewer.brew_tags()
        assert tags['TRBV1'][0] == "tgtggtcgcactgcagcaag"

class TestFindUndecombinable:

    def test_find_undecombinable_j(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.JBrewer(chain, species)
        undecombinable = brewer.find_undecombinable()
        assert undecombinable == set()

    def test_find_undecombinable_v(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.VBrewer(chain, species)
        undecombinable = brewer.find_undecombinable()
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
        assert not ("TRAJ51" in tags.keys())
        assert "TRAJ10" in tags.keys()

class TestCheckOverAlleles:

    def test_check_over_alleles(self):
        chain = "B"
        species = "Homo sapiens"
        brewer = brewers.JBrewer(chain, species)
        possible_tags = ["abc", "def"]
        alleles = {"01": "abcdef", "02": "defghi"}
        valid_tags = brewer.check_over_alleles(possible_tags, alleles)
        assert valid_tags == ["def"]

