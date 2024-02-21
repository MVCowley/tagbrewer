# Script to check how many genes on IMGT are not in current Decombinator tag dictionary

# Pseudocode
# 1. Download gene names from IMGT
# 2. Compare set of genes from IMGT and Decombinator tag dictionary
# 3. Compare functionality properties

from tagbrewer.tag import generators
from tagbrewer.utils import categories
import collections

def collapse_alleles(alleles_dict: dict) -> list:
    return [i for i in alleles_dict.keys()]

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

def find_gene_differences(directory: str, species: str, version: str) -> set:
    if species == "human":
        species_fmt1 = "Homo sapiens"
        species_fmt2 = "human"
    elif species == "mouse":
        species_fmt1 = "Mus musculus"
        species_fmt2 = "mouse"

    missing_genes = {}
    for gene_group in categories.get_tr_gene_groups():
            try:
                alleles_functionality, alleles_fastas = generators.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
                imgt_gene_list = set(collapse_alleles(alleles_fastas))
                dcr_gene_list = extract_gene_list(directory, species_fmt2, version, gene_group)
                # TODO: fix diff bug
                print(gene_group)
                print(imgt_gene_list)
                print(dcr_gene_list)
                # diff = dcr_gene_list - set(imgt_gene_list)
                diff = dcr_gene_list.difference(imgt_gene_list)
                print(diff)
                if len(diff) == 0:
                    missing_genes[gene_group] = None
                else:
                    missing_genes[gene_group] = diff
            except:
                missing_genes[gene_group] = "Not present in current decombinator"

    return missing_genes

def find_taggable_gene_differences(directory: str, species: str, version: str, tag_len=20) -> set:
    if species == "human":
        species_fmt1 = "Homo sapiens"
        species_fmt2 = "human"
    elif species == "mouse":
        species_fmt1 = "Mus musculus"
        species_fmt2 = "mouse"

    missing_genes = {}
    for gene_group in categories.get_tr_gene_groups():
            try:
                alleles_functionality, alleles_fastas = generators.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
                gene_group_tags = generators.gen_tags(alleles_fastas, tag_len)
                unique_tags = generators.find_unique_tags(gene_group_tags)
                undecombinable_genes = generators.find_undecombinable_genes(alleles_fastas, unique_tags)
                dcr_gene_list = extract_gene_list(directory, species_fmt2, version, gene_group)
                diff = dcr_gene_list - set(unique_tags)
                if len(diff) == 0:
                    missing_genes[gene_group] = None
                else:
                    missing_genes[gene_group] = diff
            except:
                missing_genes[gene_group] = "Not present in current decombinator"

    return missing_genes