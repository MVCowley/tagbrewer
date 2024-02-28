from Bio.Seq import Seq

# TODO: Alter get_c_region_post_rt to use IMGT lookup

def subtract_primer(gene_nucleotides, primer):
    rc_primer = Seq(primer).reverse_complement()
    primer_start = gene_nucleotides.index(str(rc_primer))
    post_rt_gene = gene_nucleotides[:primer_start]
    return post_rt_gene.replace("n", "")

def get_c_region_post_rt(species: str, chain: str) -> str:
    if "Homo sapiens":
        if chain == "A":
            return "atatccagaaccctgaccctgc" # TRAC*01 (first n base dropped)
        elif chain == "B":
            trbc1 = (
                "gaggacctgaacaaggtgttcccacccgaggtcgctgtgtttgagccatcagaagcagag"
                "atctcccacacccaaaaggccacactggtgtgcctggccacaggcttcttccccgaccac"
                "gtggagctgagctggtgggtgaatgggaaggaggtgcacagtggggtcagcacggacccg"
                "cagcccctcaaggagcagcccgccctcaatgactccagatactgcctgagcagccgcctg"
                "agggtctcggccaccttctggcagaacccccgcaaccacttccgctgtcaagtccagttc"
                "tacgggctctcggagaatgacgagtggacccaggatagggccaaacccgtcacccagatc"
                "gtcagcgccgaggcctggggtagagca"
            )
            trbc2 = (
                "gaggacctgaaaaacgtgttcccacccgaggtcgctgtgtttgagccatcagaagcagag"
                "atctcccacacccaaaaggccacactggtgtgcctggccacaggcttctaccccgaccac"
                "gtggagctgagctggtgggtgaatgggaaggaggtgcacagtggggtcagcacagacccg"
                "cagcccctcaaggagcagcccgccctcaatgactccagatactgcctgagcagccgcctg"
                "agggtctcggccaccttctggcagaacccccgcaaccacttccgctgtcaagtccagttc"
                "tacgggctctcggagaatgacgagtggacccaggatagggccaaacctgtcacccagatc"
                "gtcagcgccgaggcctggggtagagca"
            )
            primer = "acacagcgacctcgggtgggaa"
            hsb1 = subtract_primer(trbc1, primer) # TTCCCACCCGAGGTCGCTGTGT
            hsb2 = subtract_primer(trbc2, primer)

            if len(hsb1) >= len(hsb2):
                return hsb1
            elif len(hsb2) > len(hsb1):
                return hsb2

        elif any(chain in s for s in ["D", "G"]):
            ValueError("D and G chains have no primer at present.")
    elif "Mus musculus":
        if chain == "A":
            gene_nucleotides = (
                "nacatccagaacccagaacctgctgtgtaccagttaaaagatcctcggtctcaggacagc"
                "accctctgcctgttcaccgactttgactcccaaatcaatgtgccgaaaaccatggaatct"
                "ggaacgttcatcactgacaaaactgtgctggacatgaaagctatggattccaagagcaat"
                "ggggccattgcctggagcaaccagacaagcttcacctgccaagatatcttcaaagagacc"
                "aacgccacctaccccagttca"
            ) # TRAC*01|Mus musculus_B10.D2-H2dm1
            primer = "gagaccgaggatcttttaactgg"
        if chain == "B":
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
        elif any(chain in s for s in ["D", "G"]):
            ValueError("D and G chains have no primer at present.")
        
        return subtract_primer(gene_nucleotides, primer)
    
def get_index_oligo(number: int) -> str:
    if number == 1:
        return "atcacg"
    elif number == 8:
        return "atcacgac"