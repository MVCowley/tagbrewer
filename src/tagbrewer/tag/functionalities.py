from tagbrewer.tag import checkers, generators
from tagbrewer.utils import categories

def create_db_functionality_obj(alleles_functionality, unique_tags):
    functionality = {}
    for gene in unique_tags.keys():
        functionality[gene] = alleles_functionality[gene]["01"]

    return functionality

def extract_gene_func(directory: str, species: str, version: str, gene_groups: str) -> dict:

# TODO funct not yet finished
    filename = f"{directory}{species}_{version}_{gene_groups}.translate"

    with open(filename, "r") as file:
        translate_lines = [line.rstrip() for line in file]

    functionality = {}
    for translate_line in translate_lines:
        if version == "original":
            split = translate_line.split(",")
            gene_name = split[0]
            functionality[gene_name] = split
        elif version == "extended":
            gene_name = translate_line.split("|")[1]
        else:
            print("Version must be either original or extended")
            return ValueError
    
    return functionality

def find_functionality_differences(directory: str, species: str, version: str, tag_len=20) -> set:
    if species == "human":
        species_fmt1 = "Homo sapiens"
        species_fmt2 = "human"
    elif species == "mouse":
        species_fmt1 = "Mus musculus"
        species_fmt2 = "mouse"

    missing_genes = {}
    for gene_group in categories.get_tr_gene_groups():
            try:
                alleles_functionality, alleles_fastas = generators.get_tr_alleles_for_gene_group_for_species(gene_group, species_fmt1)
                gene_group_tags = generators.gen_tags(alleles_fastas, tag_len)
                unique_tags = generators.find_unique_tags(gene_group_tags)
                dcr_gene_list = checkers.extract_gene_list(directory, species_fmt2, version, gene_group)
                diff = dcr_gene_list - set(unique_tags)
                if len(diff) == 0:
                    missing_genes[gene_group] = None
                else:
                    missing_genes[gene_group] = diff
            except:
                missing_genes[gene_group] = "Not present in current decombinator"

    return missing_genes