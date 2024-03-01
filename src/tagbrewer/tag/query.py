# TODO: add functionality which creates files in Decombinator expected format.

import requests
import collections
from tagbrewer.utils import strings, sequences
import pandas as pd
from bs4 import BeautifulSoup
import itertools
import pathlib
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
    functionality = fields[3]

    return gene, allele_designation, functionality

def get_tr_alleles(
    chain: str, region: str, species: str
) -> Tuple[
    DefaultDict[str, DefaultDict[str, dict]], DefaultDict[str, DefaultDict[str, dict]]
]:
    
    alleles = collections.defaultdict(dict)

    with open(pathlib.Path(f"src/tagbrewer/resources/imgt_{species}_"
              f"TR{chain}{region}.html").resolve()) as file:
        parser = BeautifulSoup(file, features="html.parser")
        
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

def get_tr_alleles_remote(
    chain: str, region: str, species: str
) -> Tuple[
    DefaultDict[str, DefaultDict[str, dict]], DefaultDict[str, DefaultDict[str, dict]]
]:
    alleles = collections.defaultdict(dict)

    response = requests.get(
        f"https://www.imgt.org/genedb/GENElect?query=7.2+TR{chain}{region}&species={species}"
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