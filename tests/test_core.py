import unittest

from structuregrade.core import DEMO, analyze


class StructureGradeTest(unittest.TestCase):
    def test_demo_reports_confidence_bands(self):
        result = analyze(DEMO)
        self.assertEqual(result["residue_count"], 8)
        self.assertEqual(result["low_confidence_residues"], 2)
        self.assertEqual(result["inter_residue_clashes"], 0)
        self.assertEqual(result["grade"], "C")

    def test_rejects_empty_structure(self):
        with self.assertRaisesRegex(ValueError, "no parseable"):
            analyze({"pdb": "HEADER only"})


if __name__ == "__main__":
    unittest.main()
