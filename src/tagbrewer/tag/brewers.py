from abc import ABC, abstractmethod
from tagbrewer.utils import strings, sequences
from tagbrewer.tag import query
import collections
from typing import List, Dict, DefaultDict

# TODO: Add type hinting

READ_1_LENGTH = 150
REGION_DELS = 10

class Brewer(ABC):
    """ Parent class for Brewers """

    def __init__(self, chain: str, region: str, species: str, tag_len: int=20, functional: bool=False):
        self.functionality, self.fasta_dicts = query.get_tr_alleles(chain, region, species)
        self.chain = chain
        self.region = region
        self.species = species
        self.tag_len = tag_len

        if functional:
            self.fasta_dicts = self.drop_non_functional()

    @abstractmethod
    def slice_fasta(self, sequence: str) -> str:
        pass

    def brew_all_tags(self) -> Dict[str, List[str]]:
        """
        Generate tags for each gene
        """
        gene_group_tags = {}
        for gene, alleles in self.fasta_dicts.items():
            prototypical_fasta = alleles['01']
            sliced_fasta = self.slice_fasta(prototypical_fasta)

            possible_tags = [i for i in strings.slice_iterator(sliced_fasta, self.tag_len)]
            gene_group_tags[gene] = self.check_over_alleles(possible_tags, alleles)

        return gene_group_tags

    def brew_tags(self) -> Dict[str, List[str]]:
        all_tags = self.brew_all_tags()
        unique_tags = {gene: [] for gene in all_tags.keys()}
        for gene, possible_tags in all_tags.items():
            check_list = []
            for check_gene, check_possible_tags in all_tags.items():
                if gene != check_gene:
                    check_list.extend(check_possible_tags)
            check_set = set(check_list)
            for test_tag in possible_tags:
                if test_tag in check_set:
                    continue
                unique_tags[gene].append(test_tag)
        return unique_tags
    
    def find_undecombinable(self) -> List[str]:
        unique_tags = self.brew_tags()
        undecombinable = set([gene for gene in unique_tags.keys()
                              if len(unique_tags[gene]) == 0])
        return undecombinable
    
    def drop_non_functional(self) -> DefaultDict[str, DefaultDict[str, str]]:
        functional = collections.defaultdict(dict)
        for gene, alleles in self.fasta_dicts.items():
            for allele in alleles.keys():
                if self.functionality[gene][allele] == 'F':
                    functional[gene][allele] = self.fasta_dicts[gene][allele]
                else:
                    continue
        return functional
    
    def check_over_alleles(self, possible_tags: List[str], alleles: Dict[str, str]) -> List[str]:
        valid_tags = [possible_tag for possible_tag in possible_tags if
                      all(possible_tag in allele for allele in alleles.values())]
        return valid_tags


class VBrewer(Brewer):
    """ Class which creates V tags """

    def __init__(self, chain: str, species: str, tag_len: int=20, functional: bool=False):
        super().__init__(chain, "V", species, tag_len, functional)

    def slice_fasta(self, sequence: str) -> str:
        start_index = self.conservative_v_gene_start_index()
        sliced_fasta = sequence[start_index: -REGION_DELS]
        return sliced_fasta

    def get_max_gene_length(self, other_region: str) -> int:
        _, region_fastas = query.get_tr_alleles(self.chain, other_region, self.species)
        region_lengths = {gene: len(alleles["01"]) for gene, alleles in region_fastas.items()}
        return max(region_lengths.values())

    def conservative_v_gene_start_index(self) -> int:
        """ Returns a negative number that indexes the V gene """
        ix_length = len(sequences.get_index_oligo(1))
        c_length = len(sequences.get_c_region_post_rt(chain=self.chain, species=self.species))
        j_length = self.get_max_gene_length("J")

        if self.chain == "B":
            d_length = self.get_max_gene_length("D")
            not_v_read1 = ix_length + c_length + j_length + d_length
        else:
            not_v_read1 = ix_length + c_length + j_length

        return not_v_read1 - READ_1_LENGTH
    
class JBrewer(Brewer):
    """ Class which creates J tags """

    def __init__(self, chain: str, species: str, tag_len: int=20, functional: bool=False):
        super().__init__(chain, "J", species, tag_len, functional)

    def slice_fasta(self, sequence: str) -> str:
        sliced_fasta = sequence[REGION_DELS:]
        return sliced_fasta