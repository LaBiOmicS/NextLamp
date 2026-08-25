import unittest
from unittest.mock import patch, MagicMock
from nextlamp.alignment import filter_by_specificity
from nextlamp.candidates import CandidatePrimer

class TestMultiIndexAlignment(unittest.TestCase):
    @patch("nextlamp.alignment.subprocess.Popen")
    @patch("nextlamp.alignment.tempfile.NamedTemporaryFile")
    def test_filter_by_specificity_multi_index(self, mock_temp, mock_popen):
        # Setup mock tempfile
        mock_file = MagicMock()
        mock_file.name = "/tmp/test_temp.fa"
        mock_temp.return_value = mock_file

        # Setup candidates
        cand1 = CandidatePrimer(seq="ATCGATCGATCGATCGATCG", start=10, end=30, strand=1, tm=60.0, gc=50.0)
        candidates_dict = {
            "F3_B3": [cand1],
            "F2_B2": [],
            "F1c_B1c": [],
            "Loop": []
        }

        # Setup mock process
        mock_proc = MagicMock()
        mock_stdout = MagicMock()
        mock_stdout.__iter__.return_value = ["Q0\t0\tTarget1\t1\t255\t20M\t*\t0\t0\tATCGATCGATCGATCGATCG\t*\tNM:i:0\n"]
        mock_proc.stdout = mock_stdout
        mock_proc.wait.return_value = 0
        mock_popen.return_value = mock_proc

        indexes = ["data/idx1", "data/idx2"]
        targets_list = ["Target1"]
        background_list = ["Bg1"]

        filtered = filter_by_specificity(
            candidates_dict=candidates_dict,
            bowtie2_path="/usr/bin/bowtie2",
            index_prefix=indexes,
            targets_list=targets_list,
            background_list=background_list
        )

        self.assertIn("F3_B3", filtered)
        self.assertEqual(len(filtered["F3_B3"]), 1)
        self.assertEqual(mock_popen.call_count, 2)

if __name__ == "__main__":
    unittest.main()
