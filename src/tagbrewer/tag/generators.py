# Script to generate tags for Decombinator assignment

# Pseudocode
# 1. Pull all possible sequences for the given chain and gene_groups from tidytcells
# 2. For each sequence generate list of all possible 20bp substrings
# 3. For each substring, search all other sequence substrings for a match,
#   saving each unique substring.
# 4. Save a csv of genes + their unique tag
# 5. Save a csv of genes with no unique tags

# TODO: add functionality which creates files in Decombinator expected format.

import requests
import collections
from tagbrewer.utils import strings
import pandas as pd
from bs4 import BeautifulSoup
import itertools
from typing import Dict, FrozenSet, DefaultDict, Tuple, List

def parse_fasta_header(line: str) -> Tuple[str, str, str]:
    """
    Code from: https://github.com/yutanagano/tidytcells/blob/50af17ff1230cd3312caf14bded48987754528ef/scripts/script_utility.py#L2
    """
    fields = line.split("|")
    allele_name = fields[1]
    gene, allele_designation = allele_name.split("*")
    functionality = fields[3].strip("()[]")

    return gene, allele_designation, functionality

def get_tr_alleles_for_gene_group_for_species(gene_group: str, species: str) -> Tuple[DefaultDict, DefaultDict]:
    """
    Code from: https://github.com/yutanagano/tidytcells/blob/50af17ff1230cd3312caf14bded48987754528ef/scripts/script_utility.py#L2
    """
    alleles = collections.defaultdict(dict)

    response = requests.get(
        f"https://www.imgt.org/genedb/GENElect?query=7.2+{gene_group}&species={species}"
    )

    parser = BeautifulSoup(response.text, features="html.parser")
    fasta = parser.find_all("pre")[1].string
    header_lines = filter(lambda line: line.startswith(">"), fasta.splitlines())

    for line in header_lines:
        gene, allele_designation, functionality = parse_fasta_header(line)
        alleles[gene][allele_designation] = functionality

    # NEW: Now create code which parses fasta data
    gene_fastas = collections.defaultdict(dict)
    fasta_lines = fasta.split(">")
    for line in fasta_lines[1:]:
        fields = line.split("|")
        allele_name = fields[1]
        gene, allele_designation = allele_name.split("*")
        line_fasta = fields[-1].replace("\n", "")
        gene_fastas[gene][allele_designation] = line_fasta

    return alleles, gene_fastas

def gen_tags(fasta_dicts: DefaultDict, tag_len: int=20) -> Dict[str, List[str]]:
    """
    Generate 20bp tags from the prototypical allelel sequence for each gene
    """
    gene_group_tags = {}
    for gene, alleles in fasta_dicts.items():
        prototypical_fasta = alleles['01']
        possible_tags = [i for i in strings.sliceIterator(prototypical_fasta, tag_len)]
        gene_group_tags[gene] = possible_tags
    return gene_group_tags

def find_unique_tags(gene_group_tags: Dict[str, List]) -> Dict[str, List[str]]:

    unique_tags = collections.defaultdict(list)
    for gene, possible_tags in gene_group_tags.items():
        check_list = []
        for check_gene, check_possible_tags in gene_group_tags.items():
            if gene not in check_gene:
                check_list.extend(check_possible_tags)
        for test_tag in possible_tags:
            if test_tag in check_list:
                continue
            else:
                unique_tags[gene].append(test_tag)

    return unique_tags

def find_undecombinable_genes(alleles_fastas: DefaultDict,
                              unique_tags: DefaultDict[str, List[str]]
                              ) -> List[str]:
    
    return set(alleles_fastas) - set(unique_tags)

# if __name__ == "__main__":
#     TAG_LEN = 20 # Specify tag lengths
#     for species in strings.get_species():
#         for gene_group in strings.get_tr_gene_groups():
#             alleles_functionality, alleles_fastas = get_tr_alleles_for_gene_group_for_species(gene_group, species)
#             gene_group_tags = gen_tags(alleles_fastas, TAG_LEN)
#             unique_tags = find_unique_tags(gene_group_tags)
#             undecombinable_genes = find_undecombinable_genes(alleles_fastas, unique_tags)
#             if len(undecombinable_genes) == 0:
#                 decombineable_genes = len([i for i in unique_tags.keys()])
#                 tag_counts = [len(i) for i in unique_tags.values()]
#                 sum_tag_counts = sum(tag_counts)
#                 mean_gene_tags = sum_tag_counts / decombineable_genes
#                 std_dev_gene_tags = (sum([(i - mean_gene_tags)**2 for i in tag_counts]) / decombineable_genes)**(1/2)
#                 print(f"Mean tag counts = {mean_gene_tags} (SD) {std_dev_gene_tags} for {species}, {gene_group} across {decombineable_genes} genes.")
#             else:
#                 print(f"Genes with no valid decombinator tag for {species}, {gene_group}: {undecombinable_genes}. \
#                         {len(set(unique_tags))}/{len(set(alleles_fastas))} decombinable.")
            
            

