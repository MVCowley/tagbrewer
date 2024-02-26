# TODO: add functionality which creates files in Decombinator expected format.

import requests
import collections
from tagbrewer.utils import strings, sequences
import pandas as pd
from bs4 import BeautifulSoup
import itertools
from typing import Dict, FrozenSet, DefaultDict, Tuple, List

# TODO: create logic in brewers to make j genes and seperate filtered v genes

def gen_tags(fasta_dicts: DefaultDict, tag_len: int=20) -> Dict[str, List[str]]:
    """
    Generate 20bp tags from the prototypical allelel sequence for each gene
    """
    gene_group_tags = {}
    for gene, alleles in fasta_dicts.items():
        prototypical_fasta = alleles['01']
        possible_tags = [i for i in strings.slice_iterator(prototypical_fasta, tag_len)]
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