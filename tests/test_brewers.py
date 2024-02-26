from tagbrewer.tag import brewers

TAG_LEN = 20

class TestMaxGeneLengths:

    def test_get_max_gene_length(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.VBrewer(chain, species, tag_len)
        max_gene_length = brewer.get_max_gene_length("D")
        assert max_gene_length == len("gggactagcggggggg") # TRBD2

class TestVStartIndex:

    def test_conservative_v_gene_start_index(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.VBrewer(chain, species, tag_len)
        start_index = brewer.conservative_v_gene_start_index()
        assert start_index == -57

class TestBrewAllTags:

    def test_jbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.JBrewer(chain, species, tag_len)
        all_tags = brewer.brew_all_tags()
        assert all_tags['TRBJ1-1'][0] == "tgaacactgaagctttcttt"

    def test_vbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.VBrewer(chain, species, tag_len)
        all_tags = brewer.brew_all_tags()
        print(all_tags)
        assert all_tags['TRBV1'][0] == "tgtggtcgcactgcagcaag"

class TestBrewTags:

    def test_jbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.JBrewer(chain, species, tag_len)
        tags = brewer.brew_tags()
        assert tags['TRBJ1-1'][0] == "tgaacactgaagctttcttt"

    def test_vbrewer(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.VBrewer(chain, species, tag_len)
        tags = brewer.brew_tags()
        assert tags['TRBV1'][0] == "tgtggtcgcactgcagcaag"

class TestFindUndecombinable:

    def test_find_undecombinable_j(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.JBrewer(chain, species, tag_len)
        undecombinable = brewer.find_undecombinable()
        assert undecombinable == set()

    def test_find_undecombinable_v(self):
        chain = "B"
        species = "Homo sapiens"
        tag_len = TAG_LEN
        brewer = brewers.VBrewer(chain, species, tag_len)
        undecombinable = brewer.find_undecombinable()
        assert undecombinable == set(('TRBV3-2',
                                      'TRBV6-2',
                                      'TRBV12-4',
                                      'TRBV24/OR9-2',
                                      'TRBV3-1',
                                      'TRBV20/OR9-2',
                                      'TRBV24-1',
                                      'TRBV6-3',
                                      'TRBV20-1',
                                      'TRBV12-3'))