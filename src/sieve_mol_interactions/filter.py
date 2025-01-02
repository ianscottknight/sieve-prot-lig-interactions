from __future__ import annotations
from typing import TYPE_CHECKING, List, Optional, Tuple, Dict, Union
if TYPE_CHECKING:
    pass

from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport


class InteractionFilter:
    def __init__(self, interaction_criteria: Dict[str, Dict[str, int]]) -> None:
        """
        Args:
            interaction_criteria: A dictionary where keys are interaction types (e.g., 'hbonds') and values
                                  are dictionaries with optional 'min' and 'max' keys, e.g.,
                                  {'hbonds': {'min': 3, 'max': 10}, 'hydrophobic': {'min': 5}}.
        """
        self.interaction_criteria = interaction_criteria

    def passes(self, pdb_complex_path: str, ligand_id: str = "LIG") -> bool:
        """
        Filter a molecule based on specified interaction criteria.

        Returns:
            A boolean indicating whether the molecule meets all specified interaction criteria.
        """
        return self.ligand_obeys_interaction_criteria(pdb_complex_path, self.interaction_criteria, ligand_id=ligand_id)

    @staticmethod
    def ligand_obeys_interaction_criteria(
            pdb_complex_path: str,
            interaction_criteria: Dict[str, Dict[str, int]],
            ligand_id: str = "LIG",
    ) -> bool:
        """
        Determine whether specified ligand obeys given interaction criteria.

        Returns:
            A boolean indicating whether the ligand obeys the given interaction criteria.
        """
        interactions_dict = InteractionFilter.analyze_protein_ligand_interactions(pdb_complex_path, ligand_id)

        # Check if all criteria are met
        for interaction_type, criterion in interaction_criteria.items():
            count = interactions_dict.get(interaction_type, 0)
            min_count = criterion.get('min', 0)  # Default to zero
            max_count = criterion.get('max', float('inf'))  # Default to infinity
            if not (min_count <= count <= max_count):
                return False

        return True

    @staticmethod
    def analyze_protein_ligand_interactions(pdb_complex_path: str, ligand_id: str = "LIG") -> Optional[Dict[str, int]]:
        """
        Analyzes interactions for a specific ligand within a PDB complex file.
        """
        # Initialize the PDB complex
        complex = PDBComplex()
        complex.load_pdb(pdb_complex_path)
        complex.analyze()

        # Find the specified ligand and analyze interactions
        ligand_id_found = False
        for ligand in complex.ligands:
            if ligand.hetid == ligand_id:
                ligand_id_found = True
                complex.characterize_complex(ligand)
                site = complex.interaction_sets[ligand.mol.title]
                if site:
                    binding_site = BindingSiteReport(site)
                    interactions_dict = {
                        'hydrophobic': len(binding_site.hydrophobic_info),
                        'hbond': len(binding_site.hbond_info),
                        'waterbridge': len(binding_site.waterbridge_info),
                        'saltbridge': len(binding_site.saltbridge_info),
                        'pistacking': len(binding_site.pistacking_info),
                        'pication': len(binding_site.pication_info),
                        'halogen': len(binding_site.halogen_info),
                        'metal': len(binding_site.metal_info),
                    }
                    interactions_dict['total_interactions'] = sum(interactions_dict.values())
                    return interactions_dict

        if not ligand_id_found:
            raise ValueError(f"No interactions found for ligand {ligand_id} in the PDB complex.")

        return None
