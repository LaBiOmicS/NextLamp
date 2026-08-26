"""
Unit tests suite for PatentLAMP modules.
"""

import os
import unittest
from patentlamp.chimerization import (
    ChimerizationEngine,
    calculate_gc,
    calculate_nearest_neighbor_dg,
    calculate_full_matrix_homodimer_dg,
    parse_and_split_fip_bip,
    clean_iupac_sequence
)
from patentlamp.inventive_step import InventiveStepProofEngine
from patentlamp.wipo_st26 import WIPOSequenceListingGenerator

class TestPatentLAMP(unittest.TestCase):

    def setUp(self):
        self.sample_primer_set = {
            "F3": "GCTTGTCGTTAACTAGCG",
            "B3": "CGATCGATCGATCGAT",
            "F1c": "AGCTAGCTAGCTAGCTAGCT",
            "F2": "GCTAGCTAGCTAGCTA",
            "B1c": "TCTAGCTAGCTAGCTAGCTA",
            "B2": "ATCGATCGATCGATCG",
            "LoopF": "GCTAGCTAGCTA",
            "LoopB": "TAGCTAGCTAGC"
        }
        self.chimer_engine = ChimerizationEngine(default_linker="TAAA")

    def test_clean_iupac_sequence(self):
        # R -> A, Y -> C, K -> G
        self.assertEqual(clean_iupac_sequence("AGCTRYK"), "AGCTACG")

    def test_calculate_gc(self):
        self.assertEqual(calculate_gc("GGCC"), 100.0)
        self.assertEqual(calculate_gc("AATT"), 0.0)
        self.assertEqual(calculate_gc(""), 0.0)

    def test_nearest_neighbor_dg(self):
        dg = calculate_nearest_neighbor_dg("GCGCGCGC")
        self.assertIsInstance(dg, float)
        self.assertLess(dg, 0.0)

    def test_full_matrix_homodimer_dg(self):
        dg = calculate_full_matrix_homodimer_dg("AGCTAGCTAGCTAGCTAGCTTAAAGCTAGCTAGCTAGCTA")
        self.assertIsInstance(dg, float)

    def test_parse_and_split_fip_bip_fallback(self):
        fused_set = {
            "FIP": "AGCTAGCTAGCTAGCTAGCTGCTAGCTAGCTAGCTA",
            "BIP": "TCTAGCTAGCTAGCTAGCTAATCGATCGATCGATCG"
        }
        f1c, f2, b1c, b2 = parse_and_split_fip_bip(fused_set)
        self.assertTrue(len(f1c) > 0)
        self.assertTrue(len(f2) > 0)
        self.assertEqual(f1c + f2, fused_set["FIP"])

    def test_chimerization_engine(self):
        res = self.chimer_engine.process_primer_set(self.sample_primer_set)
        self.assertIn("synthetic_chimeric_set", res)
        fip_syn = res["synthetic_chimeric_set"]["FIP_synthetic"]
        self.assertIn("TAAA", fip_syn)

    def test_inventive_step_engine(self):
        chimer_res = self.chimer_engine.process_primer_set(self.sample_primer_set)
        proof_engine = InventiveStepProofEngine(target_species="Babesia canis")
        eval_res = proof_engine.evaluate_comparative_advantage(chimer_res)
        self.assertIn("boltzmann_equilibrium_ratio", eval_res["comparative_metrics"])

    def test_wipo_st26_generator(self):
        chimer_res = self.chimer_engine.process_primer_set(self.sample_primer_set)
        proof_engine = InventiveStepProofEngine(target_species="Babesia canis")
        eval_res = proof_engine.evaluate_comparative_advantage(chimer_res)
        gen = WIPOSequenceListingGenerator()
        out_xml = "/tmp/test_st26.xml"
        out_md = "/tmp/test_patent.md"
        gen.generate_st26_xml([chimer_res], out_xml)
        gen.generate_patent_draft_md([eval_res], out_md)
        self.assertTrue(os.path.exists(out_xml))
        self.assertTrue(os.path.exists(out_md))
        if os.path.exists(out_xml): os.remove(out_xml)
        if os.path.exists(out_md): os.remove(out_md)

if __name__ == "__main__":
    unittest.main()
