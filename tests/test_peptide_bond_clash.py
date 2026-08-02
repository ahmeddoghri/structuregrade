import unittest

from structuregrade.core import DEMO, analyze


def _atom(serial, name, residue, chain, number, x, y, z, conf):
    return (
        f"ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{conf:6.2f}          {name[0]:>2s}"
    )


def _two_residue_backbone():
    """A chemically normal two-residue backbone with standard bond lengths,
    including a ~1.33 A peptide bond (C of residue 1 to N of residue 2) --
    a covalent bond, not a steric clash."""
    return "\n".join(
        [
            _atom(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, 95),
            _atom(2, "CA", "ALA", "A", 1, 1.458, 0.000, 0.000, 95),
            _atom(3, "C", "ALA", "A", 1, 2.009, 1.420, 0.000, 95),
            _atom(4, "O", "ALA", "A", 1, 1.251, 2.397, 0.000, 95),
            _atom(5, "N", "GLY", "A", 2, 3.328, 1.470, 0.000, 92),
            _atom(6, "CA", "GLY", "A", 2, 4.008, 2.760, 0.000, 92),
            _atom(7, "C", "GLY", "A", 2, 5.520, 2.640, 0.000, 92),
            _atom(8, "O", "GLY", "A", 2, 6.150, 1.590, 0.000, 92),
        ]
    )


class PeptideBondClashTest(unittest.TestCase):
    """The clash detector excluded same-residue atom pairs but not the
    canonical C(i)-N(i+1) peptide bond spanning adjacent residues -- a
    normal covalent bond (~1.33 A) well under the 1.5 A clash threshold.
    Any real multi-atom (non-CA-only) structure would have every adjacent
    residue pair flagged as a false clash. The bundled demo never exposed
    this because it uses a single CA atom per residue."""

    def test_normal_peptide_bond_is_not_flagged_as_a_clash(self):
        result = analyze({"pdb": _two_residue_backbone()})
        self.assertEqual(result["inter_residue_clashes"], 0)
        self.assertEqual(result["grade"], "A")

    def test_genuine_non_bonded_clash_between_distant_residues_is_still_caught(self):
        pdb = _two_residue_backbone() + "\n" + _atom(9, "CA", "VAL", "A", 5, 0.100, 0.050, 0.000, 90)
        result = analyze({"pdb": pdb})
        self.assertGreaterEqual(result["inter_residue_clashes"], 1)

    def test_demo_fixture_unaffected(self):
        result = analyze(DEMO)
        self.assertEqual(result["grade"], "C")
        self.assertEqual(result["inter_residue_clashes"], 0)


if __name__ == "__main__":
    unittest.main()
