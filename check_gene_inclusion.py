# Script to check how many genes on IMGT are not in current Decombinator tag dictionary

# Pseudocode
# 1. Download gene names from IMGT
# 2. Assign taggable and functionality properties to each gene
# 3. Compare set of genes from IMG and Decombinator tag dictionary
# 4. Compare taggable properties
# 5. Compare functionality properties

import tag_gen
import collections

def get_alt_species():
    return ("human", "mouse")

def extract_gene_list(directory: str, species: str, version: str, gene_groups: str) -> set:

    filename = f"{directory}{species}_{version}_{gene_groups}.tags"

    with open(filename, "r") as file:
        lines = [line.rstrip() for line in file]

    gene_names = []
    for line in lines:
        if version == "original":
            gene_names.append(",".join(line.split("|")[1].split("/")))
        elif version == "extended":
            gene_names.append(line.split("|")[1])
        else:
            print("Version must be either original or extended")
            return ValueError
    
    return set(gene_names)
    
def find_gene_differences(directory: str, species: str, version: str, tag_len=20) -> set:
    if species == "human":
        species_fmt1 = "Homo sapiens"
        species_fmt2 = "human"
    elif species == "mouse":
        species_fmt1 = "Mus musculus"
        species_fmt2 = "mouse"

    missing_genes = {}
    for gene_group in tag_gen.get_tr_gene_groups():
            try:
                alleles_functionality, alleles_fastas = tag_gen.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
                gene_group_tags = tag_gen.gen_tags(alleles_fastas, tag_len)
                unique_tags = tag_gen.find_unique_tags(gene_group_tags)
                undecombinable_genes = tag_gen.find_undecombinable_genes(alleles_fastas, unique_tags)
                dcr_gene_list = extract_gene_list(directory, species_fmt2, version, gene_group)
                diff = dcr_gene_list - set(unique_tags)
                if len(diff) == 0:
                    missing_genes[gene_group] = None
                else:
                    missing_genes[gene_group] = diff
            except:
                missing_genes[gene_group] = "Not present in current decombinator"

    return missing_genes

# if __name__ == "__main__":
#     TAG_LEN = 20 # Specify tag lengths
#     for species_fmt1, species_fmt2 in zip(tag_gen.get_species(), get_alt_species()):
#         for gene_group in tag_gen.get_tr_gene_groups():
#             alleles_functionality, alleles_fastas = tag_gen.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
#             gene_group_tags = tag_gen.gen_tags(alleles_fastas, TAG_LEN)
#             unique_tags = tag_gen.find_unique_tags(gene_group_tags)
#             undecombinable_genes = tag_gen.find_undecombinable_genes(alleles_fastas, unique_tags)
#             dcr_gene_list = extract_gene_list("Decombinator-Tags-FASTAs/", species_fmt2, "original", gene_group)
#             print(f"Genes missing from {species_fmt1 + '/' + species_fmt2}, {gene_group}: {dcr_gene_list - set(unique_tags)}")