from typing import List, Dict

class Tags:
    def __init__(self, gene_layers: Dict[str, Dict[str, List[str]]]) -> None:
        self.gene_layers = gene_layers