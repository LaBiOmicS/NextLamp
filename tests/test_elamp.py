import unittest
from nextlamp.elamp import evaluate_primer_set_quality, simulate_elamp_amplicon, calculate_3end_deltag, calculate_3end_gc_clamp

class TestELAMPModule(unittest.TestCase):
    def test_evaluate_primer_set_quality(self):
        sample_set = {
            "f3": "ATCGATCGATCGATCG",
            "b3": "CGATCGATCGATCGAT",
            "f2": "ATCGATCGATCGATCGATCG",
            "b2": "CGATCGATCGATCGATCGAT",
            "f1c": "GCTAGCTAGCTAGCTAGCTA",
            "b1c": "TAGCTAGCTAGCTAGCTAGC",
            "tm_f3": 60.0,
            "tm_b3": 60.0,
            "tm_f2": 60.0,
            "tm_b2": 60.0,
            "tm_f1c": 62.0,
            "tm_b1c": 62.0
        }

        metrics = evaluate_primer_set_quality(sample_set)
        self.assertIn("quality_score", metrics)
        self.assertGreaterEqual(metrics["quality_score"], 80.0)
        self.assertEqual(metrics["tm_balance"], 0.0)
        self.assertEqual(metrics["diff_f1c_b1c"], 0.0)

    def test_simulate_elamp_amplicon(self):
        target_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        sample_set = {
            "f3": "ATCGATCGATCGATCG",
            "b3": "CGATCGATCGATCGAT",
            "f2": "ATCGATCGATCGATCGATCG",
            "b2": "CGATCGATCGATCGATCGAT"
        }

        sim_results = simulate_elamp_amplicon(target_seq, sample_set)
        self.assertIn("quality_metrics", sim_results)
        self.assertIn("outer_amplicon_size", sim_results)

if __name__ == "__main__":
    unittest.main()
