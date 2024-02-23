from tagbrewer.utils import sequences

def test_subtract_primer_hsa():
    gene_nucleotides = (
        "natatccagaaccctgaccctgccgtgtaccagctgagagactctaaatccagtgacaag"
        "tctgtctgcctattcaccgattttgattctcaaacaaatgtgtcacaaagtaaggattct"
        "gatgtgtatatcacagacaaaactgtgctagacatgaggtctatggacttcaagagcaac"
        "agtgctgtggcctggagcaacaaatctgactttgc"
        ) # TRAC*01
    primer = "gagtctctcagctggtacacg"
    formatted_gene_nucleotides = sequences.subtract_primer(gene_nucleotides, primer)
    assert formatted_gene_nucleotides == "atatccagaaccctgaccctgc"

def test_subtract_primer_hsb():
    gene_nucleotides = (
        "gaggacctgaacaaggtgttcccacccgaggtcgctgtgtttgagccatcagaagcagag"
        "atctcccacacccaaaaggccacactggtgtgcctggccacaggcttcttccccgaccac"
        "gtggagctgagctggtgggtgaatgggaaggaggtgcacagtggggtcagcacggacccg"
        "cagcccctcaaggagcagcccgccctcaatgactccagatactgcctgagcagccgcctg"
        "agggtctcggccaccttctggcagaacccccgcaaccacttccgctgtcaagtccagttc"
        "tacgggctctcggagaatgacgagtggacccaggatagggccaaacccgtcacccagatc"
        "gtcagcgccgaggcctggggtagagca"
        ) # TRBC*01
    primer = "acacagcgacctcgggtgggaa"
    formatted_gene_nucleotides = sequences.subtract_primer(gene_nucleotides, primer)
    assert formatted_gene_nucleotides == "gaggacctgaacaaggtg"

def test_subtract_primer_mma():
    gene_nucleotides = (
        "nacatccagaacccagaacctgctgtgtaccagttaaaagatcctcggtctcaggacagc"
        "accctctgcctgttcaccgactttgactcccaaatcaatgtgccgaaaaccatggaatct"
        "ggaacgttcatcactgacaaaactgtgctggacatgaaagctatggattccaagagcaat"
        "ggggccattgcctggagcaaccagacaagcttcacctgccaagatatcttcaaagagacc"
        "aacgccacctaccccagttca"
        ) # TRAC*01|Mus musculus_B10.D2-H2dm1
    primer = "gagaccgaggatcttttaactgg"
    formatted_gene_nucleotides = sequences.subtract_primer(gene_nucleotides, primer)
    assert formatted_gene_nucleotides == "acatccagaacccagaacctgctgtgta"


def test_subtract_primer_mmb():
    gene_nucleotides = (
        "naggatctgagaaatgtgactccacccaaggtctccttgtttgagccatcaaaagcagag"
        "attgcaaacaaacaaaaggctaccctcgtgtgcttggccaggggcttcttccctgaccac"
        "gtggagctgagctggtgggtgaatggcaaggaggtccacagtggggtcagcacggaccct"
        "caggcctacaaggagagcaattatagctactgcctgagcagccgcctgagggtctctgct"
        "accttctggcacaatcctcgcaaccacttccgctgccaagtgcagttccatgggctttca"
        "gaggaggacaagtggccagagggctcacccaaacctgtcacacagaacatcagtgcagag"
        "gcctggggccgagca"
        ) # TRBC1*01|Mus musculus_B10.A
    primer = "gcttttgatggctcaaacaagg"
    formatted_gene_nucleotides = sequences.subtract_primer(gene_nucleotides, primer)
    assert formatted_gene_nucleotides == "aggatctgagaaatgtgactccacccaaggtct"

def test_get_c_region_post_rt():
    species = "Homo sapiens"
    chain = "B"
    c_region_post_rt = sequences.get_c_region_post_rt(species, chain)
    assert c_region_post_rt == "gaggacctgaacaaggtg"