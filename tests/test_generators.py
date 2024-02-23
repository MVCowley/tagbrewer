from tagbrewer.tag import generators
import collections

TAG_LEN = 20

def test_get_tr_alleles_for_gene_group_for_species_functionality():
    species = "Homo sapiens"
    gene_group = "TRBD"
    alleles_functionality, _ = \
        generators.get_tr_alleles_for_gene_group_for_species(gene_group, species)
    assert alleles_functionality['TRBD1']["01"] == "F"

def test_get_tr_alleles_for_gene_group_for_species_nucleotides():
    species = "Homo sapiens"
    gene_group = "TRBD"
    _, alleles_fastas = \
        generators.get_tr_alleles_for_gene_group_for_species(gene_group, species)
    assert alleles_fastas['TRBD1']["01"] == "gggacagggggc"

def test_get_max_gene_length():
    region = "D"
    chain = "B"
    species = "Homo sapiens"
    max_gene_length = generators.get_max_gene_length(region, chain, species)
    assert max_gene_length == len("gggactagcggggggg") # TRBD2

def test_conservative_v_gene_start_index():
    chain = "B"
    species = "Homo sapiens"
    start_index = generators.conservative_v_gene_start_index(chain, species)
    assert start_index == -57

def test_gen_tags_tags_found():
    """ Implicitly tests `sliceIterator` """
    alleles_fastas = collections.defaultdict(dict)
    alleles_fastas["gene1"]["01"] = "abcdeabcdeabcdeabcdeabcde"
    gene_group_tags = generators.gen_tags(alleles_fastas, TAG_LEN)
    tags_found = len(gene_group_tags["gene1"])
    assert tags_found == 6, f"value was {tags_found}, should be 6"

def test_gen_tags_tag_len():
    alleles_fastas = collections.defaultdict(dict)
    alleles_fastas["gene1"]["01"] = "abcdeabcdeabcdeabcdeabcde"
    gene_group_tags = generators.gen_tags(alleles_fastas, TAG_LEN)
    tag_len = len(gene_group_tags["gene1"][0])
    assert tag_len == 20, f"value was {tag_len}, should be 6"

def test_gen_tags_allele():
    alleles_fastas = collections.defaultdict(dict)
    alleles_fastas["gene1"]["01"] = "abc"
    alleles_fastas["gene1"]["02"] = "fghijfghijfghijfghijfghij"
    gene_group_tags = generators.gen_tags(alleles_fastas, TAG_LEN)
    assert not gene_group_tags["gene1"], f"gene group tags found, should return none"

def test_find_unique_tags_fail():
    gene_group_tags = {"gene1": ["abc"], "gene2": ["abc"]}
    unique_tags = generators.find_unique_tags(gene_group_tags)
    assert (not unique_tags["gene1"]) and (not unique_tags["gene2"])

def test_find_unique_tags_part_fail():
    gene_group_tags = {"gene1": ["abc", "def"], "gene2": ["abc"]}
    unique_tags = generators.find_unique_tags(gene_group_tags)
    assert (unique_tags["gene1"] == ["def"]) and (not unique_tags["gene2"])

def test_find_unique_tags_pass():
    gene_group_tags = {"gene1": ["abc"], "gene2": ["def"]}
    unique_tags = generators.find_unique_tags(gene_group_tags)
    assert (unique_tags["gene1"] == ["abc"]) and (unique_tags["gene2"] == ["def"])

def test_undecombinable_genes_none():
    alleles_fastas = {"gene1": ["abcdef"], "gene2": ["defghi"]}
    unique_tags = {"gene1": ["abc"], "gene2": ["def"]}
    undecombinable_genes = generators.find_undecombinable_genes(alleles_fastas, unique_tags)
    assert len(undecombinable_genes) == 0

def test_undecombinable_genes_two():
    alleles_fastas = {"gene1": ["abcdef"], "gene2": ["defghi"], "gene3": ["defghi"]}
    unique_tags = {"gene1": ["abc"]}
    undecombinable_genes = generators.find_undecombinable_genes(alleles_fastas, unique_tags)
    assert len(undecombinable_genes) == 2