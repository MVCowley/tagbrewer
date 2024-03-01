from abc import ABC, abstractmethod
from tagbrewer.utils import strings, sequences
from tagbrewer.tag import containers, query
import collections
from typing import List, Dict, DefaultDict, Set

READ_1_LENGTH = 150
REGION_DELS = 10

class Brewer(ABC):
    """ Parent class for Brewers """

    def __init__(self, chain: str, region: str, species: str, tag_len: int=20, functional: bool=False) -> None:
        self.functionality, self.fasta_dicts = query.get_tr_alleles(chain, region, species)
        self.fasta_dicts = self.drop_incomplete()

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
    
    def brew_tags(self) -> containers.Tags:
        all_tags = self.brew_all_tags()
        gene_layers = {gene: [] for gene in all_tags.keys()}
        for gene, possible_tags in all_tags.items():
            min_hits = float('inf')
            check_genes = {key: all_tags[key] for key in all_tags if key != gene}
            unique_tags = []
            tag_layers = {}
            count = 0
            best_tag = ""
            if not possible_tags:
                tag_layers[count] = "No tags cover all confirmed alleles"
                gene_layers[gene] = tag_layers
                continue
            while not unique_tags:
                check_list = self.create_check_list(gene, check_genes)
                for test_tag in possible_tags:
                    if test_tag not in check_list:
                        unique_tags.append(test_tag)
                    else:
                        if not unique_tags:
                            tag_hits = check_list.count(test_tag)
                            if tag_hits < min_hits:
                                best_tag = test_tag
                                min_hits = tag_hits
                if unique_tags:
                    tag_layers[count] = unique_tags
                else:
                    tag_layers[count] = [best_tag]
                    clashes = self.find_clashes(check_genes, best_tag)
                    if set(clashes) == set(check_genes):
                        tag_layers[count+1] = f"Clashes with {clashes}"
                        break
                    check_genes = {key: check_genes[key] for key in clashes}
                count =+ 1
            gene_layers[gene] = tag_layers

        brewed_tags = containers.Tags(gene_layers)
        return brewed_tags

    def create_check_list(self, gene: str, check_genes: Dict[str, List[str]]) -> List[str]:
        check_list = []
        for check_gene, check_possible_tags in check_genes.items():
            if gene != check_gene:
                check_list.extend(check_possible_tags)
        return check_list

    def find_clashes(self, check_genes: Dict[str, List[str]], best_tag: str) -> List[str]:
        clashes = []
        for check_gene, check_gene_tags in check_genes.items():
            if best_tag in check_gene_tags:
                clashes.append(check_gene)
        return clashes
    
    def find_undecombinable(self) -> Set[str]:
        tags = self.brew_tags()
        undecombinable = {}
        for gene, layers in tags.gene_layers.items():
            n_layers = len(layers.keys())
            if type(layers[n_layers-1]) == str:
                undecombinable[gene] = layers[n_layers-1]
        return undecombinable
    
    def drop_incomplete(self) -> DefaultDict[str, DefaultDict[str, str]]:
        complete = collections.defaultdict(dict)
        for gene, alleles in self.fasta_dicts.items():
            for allele in alleles.keys():
                if any(self.functionality[gene][allele] == func for func in ["[F]", "[P]", "(F)", "(P)"]):
                    continue
                else:
                    complete[gene][allele] = self.fasta_dicts[gene][allele]
        return complete
    
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
        valid_tags = []
        for possible_tag in possible_tags:
            if all(possible_tag in allele for allele in alleles.values()):
                valid_tags.append(possible_tag)
        return valid_tags


class VBrewer(Brewer):
    """ Class which creates V tags """

    def __init__(self, chain: str, species: str, tag_len: int=20, functional: bool=False) -> None:
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

    def __init__(self, chain: str, species: str, tag_len: int=20, functional: bool=False) -> None:
        super().__init__(chain, "J", species, tag_len, functional)

    def slice_fasta(self, sequence: str) -> str:
        sliced_fasta = sequence[REGION_DELS:]
        return sliced_fasta