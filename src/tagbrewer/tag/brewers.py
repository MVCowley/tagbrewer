from abc import ABC, abstractmethod
from tagbrewer.utils import strings, sequences
from tagbrewer.tag import query
import collections

READ_1_LENGTH = 150
V_REGION_DELS = 10

class Brewer(ABC):
    """ Parent class for Brewers """

    def __init__(self, chain: str, region: str, species: str, tag_len):
        _, self.fasta_dicts = query.get_tr_alleles_for_gene_group_for_species(chain, region, species)
        self.chain = chain
        self.region = region
        self.species = species
        self.tag_len = tag_len

    @abstractmethod
    def brew_tags(self):
        print("Please create either a VBrewer or JBrewer object to use this method.")
        pass

    def find_unique_tags(self):

        unique_tags = collections.defaultdict(list)
        for gene, possible_tags in self.all_tags.items():
            check_list = []
            for check_gene, check_possible_tags in self.all_tags.items():
                if gene not in check_gene:
                    check_list.extend(check_possible_tags)
            for test_tag in possible_tags:
                if test_tag in check_list:
                    continue
                else:
                    unique_tags[gene].append(test_tag)

        self.unique_tags = unique_tags

class VBrewer(Brewer):
    """ Class which creates V tags """

    def get_max_gene_length(self, region: str, chain: str, species: str):
        _, region_fastas = query.get_tr_alleles_for_gene_group_for_species(f"TR{chain}{region}", species)
        region_lengths = {gene: len(alleles["01"]) for gene, alleles in region_fastas.items()}
        return max(region_lengths.values())

    def conservative_v_gene_start_index(self):
        """ Returns a negative number that indexes the V gene """
        i1_length = len(sequences.get_index_oligo(1))
        c_length = len(sequences.get_c_region_post_rt(chain=self.chain, species=self.species))
        j_length = self.get_max_gene_length("J", self.chain, self.species)

        if self.chain == "B":
            d_length = self.get_max_gene_length("D", self.chain, self.species)
            not_v_read1 = i1_length + c_length + j_length + d_length
        else:
            not_v_read1 = i1_length + c_length + j_length

        return not_v_read1 - READ_1_LENGTH

    def brew_tags(self):
        """
        Generate 20bp tags from the prototypical allelel sequence for each gene
        """
        gene_group_tags = {}
        for gene, alleles in self.fasta_dicts.items():
            prototypical_fasta = alleles['01']
            start_index = self.conservative_v_gene_start_index()
            sliced_fasta = prototypical_fasta[start_index: -V_REGION_DELS]
            possible_tags = [i for i in strings.sliceIterator(sliced_fasta, self.tag_len)]
            gene_group_tags[gene] = possible_tags
        self.all_tags = gene_group_tags
    
class JBrewer(Brewer):
    """ Class which creates V tags """

    def brew_tags(self):
        """
        Generate 20bp tags from the prototypical allelel sequence for each gene
        """
        gene_group_tags = {}
        for gene, alleles in self.fasta_dicts.items():
            prototypical_fasta = alleles['01']
            
            possible_tags = [i for i in strings.sliceIterator(prototypical_fasta, self.tag_len)]
            gene_group_tags[gene] = possible_tags
        self.all_tags = gene_group_tags