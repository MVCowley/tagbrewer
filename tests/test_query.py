from tagbrewer.tag import query

TAG_LEN = 20

def test_get_tr_alleles_for_gene_group_for_species_functionality():
    species = "Homo sapiens"
    chain = "B"
    region = "D"
    alleles_functionality, _ = \
        query.get_tr_alleles_for_gene_group_for_species(chain, region, species)
    assert alleles_functionality['TRBD1']["01"] == "F"

def test_get_tr_alleles_for_gene_group_for_species_nucleotides():
    species = "Homo sapiens"
    chain = "B"
    region = "D"
    _, alleles_fastas = \
        query.get_tr_alleles_for_gene_group_for_species(chain, region, species)
    assert alleles_fastas['TRBD1']["01"] == "gggacagggggc"