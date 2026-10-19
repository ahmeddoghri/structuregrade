import unittest
from pathlib import Path

from structuregrade.core import analyze

FIXTURES = Path(__file__).parent / "fixtures"


def grade(name):
    return analyze({"pdb": (FIXTURES / name).read_text()})


class RealAlphaFoldModels(unittest.TestCase):
    """Full-atom models, not the C-alpha-only demo fixture."""

    def test_a_confident_full_atom_model_has_no_false_clashes(self):
        r = grade("AF-P69905-F1-model_v6.pdb")
        # Before the peptide-bond fix this scored C with 141 "clashes":
        # one for every C(i)-N(i+1) bond in the backbone.
        self.assertEqual(r["residue_count"], 142)
        self.assertEqual(r["inter_residue_clashes"], 0)
        self.assertEqual(r["grade"], "A")

    def test_disordered_regions_pull_the_grade_down(self):
        r = grade("AF-P04637-F1-model_v6.pdb")
        # p53's transactivation and C-terminal regions are intrinsically
        # disordered, and pLDDT says so.
        self.assertEqual(r["residue_count"], 393)
        self.assertEqual(r["inter_residue_clashes"], 0)
        self.assertGreater(r["low_confidence_residues"], 100)
        self.assertEqual(r["grade"], "C")


if __name__ == "__main__":
    unittest.main()
