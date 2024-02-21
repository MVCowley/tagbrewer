from abc import ABC, abstractmethod
from tagbrewer.utils import strings

class Brewer(ABC):

    @abstractmethod
    def brew_tags(self):
        pass

class VBrewer(Brewer):

    def brew_tags(self, fasta_dicts, tag_len):
        """
        Generate 20bp tags from the prototypical allelel sequence for each gene
        """
        gene_group_tags = {}
        for gene, alleles in fasta_dicts.items():
            prototypical_fasta = alleles['01']
            possible_tags = [i for i in strings.sliceIterator(prototypical_fasta, tag_len)]
            gene_group_tags[gene] = possible_tags
        return gene_group_tags