# Script to check how many genes on IMGT are not in current Decombinator tag dictionary

# Pseudocode
# 1. Download gene names from IMGT
# 2. Compare set of genes from IMGT and Decombinator tag dictionary
# 3. Compare functionality properties

from tagbrewer.tag import generators
from tagbrewer.utils import categories
import collections
import pathlib

def collapse_alleles(alleles_dict: dict) -> list:
    return [i for i in alleles_dict.keys()]

def extract_gene_list(directory: str, species: str, version: str, gene_groups: str) -> set:

    filename = pathlib.Path(f"{directory}/{species}_{version}_{gene_groups}.tags").resolve()

    with open(filename, "r") as file:
        lines = [line.rstrip() for line in file]

    gene_names = []
    for line in lines:
        if version == "original":
            gene_names.append("/".join(line.split("|")[1].split(",")))
        elif version == "extended":
            gene_names.append(line.split("|")[1])
        else:
            print("Version must be either original or extended")
            return ValueError
    
    return set(gene_names)

def find_gene_differences(directory: str, species: str, version: str, gene_group: str) -> set:
    """ Takes an """
    if species == "human":
        species_fmt1 = "Homo sapiens"
        species_fmt2 = "human"
    elif species == "mouse":
        species_fmt1 = "Mus musculus"
        species_fmt2 = "mouse"

    try:
        alleles_functionality, alleles_fastas = generators.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
    except FileNotFoundError:
        return f"{gene_group} not present in current decombinator."
    
    imgt_gene_list = set(collapse_alleles(alleles_fastas))
    dcr_gene_list = extract_gene_list(directory, species_fmt2, version, gene_group)
    diff = dcr_gene_list.symmetric_difference(imgt_gene_list)
    if len(diff) == 0:
        return f"{gene_group} up to date."
    else:
        return diff