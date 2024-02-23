# TODO: add functionality which creates files in Decombinator expected format.

import requests
import collections
from tagbrewer.utils import strings, sequences
import pandas as pd
from bs4 import BeautifulSoup
import itertools
from typing import Dict, FrozenSet, DefaultDict, Tuple, List

READ_1_LENGTH = 150
V_REGION_DELS = 10

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

    # NEW: Code which returns FASTA data
    gene_fastas = collections.defaultdict(dict)
    fasta_lines = fasta.split(">")
    for line in fasta_lines[1:]:
        fields = line.split("|")
        allele_name = fields[1]
        gene, allele_designation = allele_name.split("*")
        line_fasta = fields[-1].replace("\n", "")
        gene_fastas[gene][allele_designation] = line_fasta

    return alleles, gene_fastas

def get_max_gene_length(region: str, chain: str, species: str):
    _, region_fastas = get_tr_alleles_for_gene_group_for_species(f"TR{chain}{region}", species)
    region_lengths = {gene: len(alleles["01"]) for gene, alleles in region_fastas.items()}
    return max(region_lengths.values())

def conservative_v_gene_start_index(chain: str, species: str):
    """ Returns a negative number that indexes the V gene """
    i1_length = len(sequences.get_index_oligo(1))
    c_length = len(sequences.get_c_region_post_rt(chain=chain, species=species))
    j_length = get_max_gene_length("J", chain, species)

    if chain == "B":
        d_length = get_max_gene_length("D", chain, species)
        not_v_read1 = i1_length + c_length + j_length + d_length
    else:
        not_v_read1 = i1_length + c_length + j_length

    return not_v_read1 - READ_1_LENGTH


# TODO: create logic in brewers to make j genes and seperate filtered v genes

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
    
    # TODO: change to symmetric set difference
    
    return set(alleles_fastas) - set(unique_tags)