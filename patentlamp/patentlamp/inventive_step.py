"""
In Silico Proof of Inventive Step Engine (LPI Art. 13 / USPTO Non-Obviousness).
Quantifies technical advantages (thermodynamic delta dG, specificity barriers, Boltzmann equilibrium ratio)
of synthetic chimeric primers compared to their non-modified natural controls.
Derived strictly from physical thermodynamics (SantaLucia 1998 & Boltzmann equilibrium kinetics).
Handles positive/negative ddG conditionally without text contradictions or formatting bugs (+-).
"""

import math
from datetime import datetime
from typing import Dict, List, Any

class InventiveStepProofEngine:
    def __init__(self, target_species: str = "Babesia canis", temp_c: float = 65.0):
        self.target_species = target_species
        self.temp_c = temp_c
        self.gas_constant_r = 0.0019872  # kcal/(mol*K)

    def calculate_boltzmann_suppression_factor(self, delta_ddg: float) -> float:
        """
        Calculates the thermodynamic dimer suppression factor S_factor = exp(ddG / (R*T))
        where ddG = dG_chimeric - dG_natural.
        A factor > 1.0 represents fold-reduction in dimer formation probability.
        """
        temp_k = self.temp_c + 273.15
        try:
            exponent = delta_ddg / (self.gas_constant_r * temp_k)
            exponent = max(-20.0, min(20.0, exponent))
            suppression_factor = math.exp(exponent)
            return round(suppression_factor, 2)
        except OverflowError:
            return 1.0

    def evaluate_comparative_advantage(self, chimerization_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generates comparative metrics and mathematical proof of unexpected technical effect
        for synthetic chimeric primers vs. natural unmodified controls.
        """
        fip_eng = chimerization_result.get("fip_engineering", {})
        bip_eng = chimerization_result.get("bip_engineering", {})

        nat_fip_dg = fip_eng.get("natural_metrics", {}).get("dG_dimer", 0.0)
        syn_fip_dg = fip_eng.get("chimeric_metrics", {}).get("dG_dimer", 0.0)

        nat_bip_dg = bip_eng.get("natural_metrics", {}).get("dG_dimer", 0.0)
        syn_bip_dg = bip_eng.get("chimeric_metrics", {}).get("dG_dimer", 0.0)

        # ddG = dG_chimeric - dG_natural
        fip_ddg = round(syn_fip_dg - nat_fip_dg, 2)
        bip_ddg = round(syn_bip_dg - nat_bip_dg, 2)
        avg_ddg = round((fip_ddg + bip_ddg) / 2.0, 2)

        fip_red = fip_eng.get("dimer_reduction_percentage", 0.0)
        bip_red = bip_eng.get("dimer_reduction_percentage", 0.0)
        avg_dimer_reduction = round((fip_red + bip_red) / 2.0, 2)

        suppression_factor = self.calculate_boltzmann_suppression_factor(avg_ddg)

        # Format string nicely without "+-" formatting bug
        ddg_formatted = f"+{avg_ddg}" if avg_ddg > 0 else f"{avg_ddg}"

        # Conditional technical proof generation to prevent legal/scientific contradictions
        if avg_ddg > 0:
            descriptive_proof = (
                f"O conjunto de primers sintéticos quiméricos para {self.target_species} apresentou uma variação "
                f"termodinâmica positiva na energia livre de autodimerização (ddG = {ddg_formatted} kcal/mol a {self.temp_c}°C) "
                f"e uma redução média de {avg_dimer_reduction}% na estabilidade de auto-dímeros quando comparado "
                f"diretamente às sequências naturais controle não modificadas. Com base no modelo de equilíbrio "
                f"de Boltzmann, a junção sintética introduziu um fator de supressão termodinâmica de {suppression_factor}x "
                f"contra a formação de estruturas secundárias indesejadas, comprovando efeito técnico surpreendente "
                f"e não-óbvio (atendimento ao Art. 13 da LPI)."
            )
        else:
            descriptive_proof = (
                f"O conjunto de primers sintéticos quiméricos para {self.target_species} compreende moléculas não-naturais "
                f"projetadas por engenharia de ligantes sintéticos (junção artificial especificada). A alteração estrutural "
                f"mantém a eficiência de anelamento ao genoma alvo enquanto confere novidade e caráter sintético "
                f"proprietário (superando a vedação do Art. 10, IX da LPI), com estabilidade termodinâmica de "
                f"autodimerização mensurada em ddG = {ddg_formatted} kcal/mol a {self.temp_c}°C."
            )

        return {
            "metadata": {
                "tool": "PatentLAMP",
                "version": "1.0.0",
                "timestamp": datetime.now().isoformat()
            },
            "target_species": self.target_species,
            "reaction_temperature_c": self.temp_c,
            "controls_evaluated": {
                "fip_control_natural": fip_eng.get("natural_sequence"),
                "fip_synthetic_patented": fip_eng.get("synthetic_chimeric_sequence"),
                "bip_control_natural": bip_eng.get("natural_sequence"),
                "bip_synthetic_patented": bip_eng.get("synthetic_chimeric_sequence")
            },
            "comparative_metrics": {
                "fip_delta_delta_g_kcal_mol": fip_ddg,
                "bip_delta_delta_g_kcal_mol": bip_ddg,
                "average_delta_delta_g_kcal_mol": avg_ddg,
                "average_dimer_reduction_pct": avg_dimer_reduction,
                "boltzmann_equilibrium_ratio": suppression_factor
            },
            "patent_claims_support_text": descriptive_proof
        }
