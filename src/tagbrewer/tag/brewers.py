from abc import ABC, abstractmethod
from tagbrewer.utils import strings, sequences
from tagbrewer.tag import query

# TODO: Add type hinting

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
    def brew_all_tags(self) -> dict:
        pass

    def brew_tags(self) -> dict:
        all_tags = self.brew_all_tags()
        unique_tags = {gene: [] for gene in all_tags.keys()}
        for gene, possible_tags in all_tags.items():
            check_list = []
            for check_gene, check_possible_tags in all_tags.items():
                if gene != check_gene:
                    check_list.extend(check_possible_tags)
            for test_tag in possible_tags:
                if test_tag in check_list:
                    continue
                else:
                    unique_tags[gene].append(test_tag)
        return unique_tags
    
    def find_undecombinable(self):
        unique_tags = self.brew_tags()
        undecombinable = set([gene for gene in unique_tags.keys()
                              if len(unique_tags[gene]) == 0])
        return undecombinable

class VBrewer(Brewer):
    """ Class which creates V tags """

    def __init__(self, chain: str, species: str, tag_len):
        super().__init__(chain, "V", species, tag_len)

    def get_max_gene_length(self, other_region):
        _, region_fastas = query.get_tr_alleles_for_gene_group_for_species(self.chain, other_region, self.species)
        region_lengths = {gene: len(alleles["01"]) for gene, alleles in region_fastas.items()}
        return max(region_lengths.values())

    def conservative_v_gene_start_index(self):
        """ Returns a negative number that indexes the V gene """
        i1_length = len(sequences.get_index_oligo(1))
        c_length = len(sequences.get_c_region_post_rt(chain=self.chain, species=self.species))
        j_length = self.get_max_gene_length("J")

        if self.chain == "B":
            d_length = self.get_max_gene_length("D")
            not_v_read1 = i1_length + c_length + j_length + d_length
        else:
            not_v_read1 = i1_length + c_length + j_length

        return not_v_read1 - READ_1_LENGTH

    def brew_all_tags(self):
        """
        Generate 20bp tags from the prototypical allelel sequence for each gene
        """
        gene_group_tags = {}
        for gene, alleles in self.fasta_dicts.items():
            prototypical_fasta = alleles['01']
            start_index = self.conservative_v_gene_start_index()
            sliced_fasta = prototypical_fasta[start_index: -V_REGION_DELS]
            possible_tags = [i for i in strings.slice_iterator(sliced_fasta, self.tag_len)]
            gene_group_tags[gene] = possible_tags
        return gene_group_tags
    
class JBrewer(Brewer):
    """ Class which creates V tags """

    def __init__(self, chain: str, species: str, tag_len):
        super().__init__(chain, "J", species, tag_len)

    def brew_all_tags(self):
        """
        Generate 20bp tags from the prototypical allelel sequence for each gene
        """
        gene_group_tags = {}
        for gene, alleles in self.fasta_dicts.items():
            prototypical_fasta = alleles['01']
            possible_tags = [i for i in strings.slice_iterator(prototypical_fasta, self.tag_len)]
            gene_group_tags[gene] = possible_tags
        return gene_group_tags