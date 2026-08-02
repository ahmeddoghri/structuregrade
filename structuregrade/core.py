from __future__ import annotations

import math
from collections import defaultdict
from typing import Any


def _atom(serial: int, atom: str, residue: str, chain: str, number: int, x: float, y: float, z: float, confidence: float) -> str:
    return f"ATOM  {serial:5d} {atom:^4s} {residue:>3s} {chain}{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{confidence:6.2f}          {atom[0]:>2s}"


DEMO = {
    "pdb": "\n".join([
        _atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, 96),
        _atom(2, "CA", "GLY", "A", 2, 3.8, 0.0, 0.0, 92),
        _atom(3, "CA", "SER", "A", 3, 7.6, 0.2, 0.0, 88),
        _atom(4, "CA", "LEU", "A", 4, 11.3, 0.4, 0.1, 82),
        _atom(5, "CA", "GLU", "A", 5, 15.1, 0.3, 0.2, 74),
        _atom(6, "CA", "LYS", "A", 6, 18.8, 0.1, 0.3, 61),
        _atom(7, "CA", "PRO", "A", 7, 22.6, 0.0, 0.2, 44),
        _atom(8, "CA", "GLY", "A", 8, 26.4, 0.1, 0.0, 38),
    ])
}


def analyze(payload: dict[str, Any]) -> dict[str, Any]:
    text = str(payload.get("pdb", ""))
    atoms = []
    residues: dict[tuple[str, int], list[float]] = defaultdict(list)
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 66:
            continue
        try:
            atom_name = line[12:16].strip()
            chain = line[21].strip() or "_"
            residue = int(line[22:26])
            coordinates = tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
            confidence = float(line[60:66])
        except ValueError:
            continue
        atoms.append((chain, residue, atom_name, coordinates))
        residues[(chain, residue)].append(confidence)
    if not residues:
        raise ValueError("input contains no parseable PDB atoms")
    scores = [sum(values) / len(values) for values in residues.values()]
    clashes = 0
    for index, (chain_a, residue_a, atom_a, point_a) in enumerate(atoms):
        for chain_b, residue_b, atom_b, point_b in atoms[index + 1:]:
            if (chain_a, residue_a) == (chain_b, residue_b):
                continue
            # Skip the canonical backbone peptide bond (C of residue i to N
            # of residue i+1 in the same chain) -- a normal covalent bond
            # (~1.33 A) that would otherwise be flagged as a steric clash on
            # every real multi-atom structure, even though nothing is wrong.
            if chain_a == chain_b and abs(residue_a - residue_b) == 1:
                lower, higher = sorted((residue_a, residue_b))
                lower_atom = atom_a if residue_a == lower else atom_b
                higher_atom = atom_b if residue_a == lower else atom_a
                if lower_atom == "C" and higher_atom == "N":
                    continue
            distance = math.dist(point_a, point_b)
            if distance < 1.5:
                clashes += 1
    high = sum(score >= 90 for score in scores)
    low = sum(score < 50 for score in scores)
    mean = sum(scores) / len(scores)
    grade = "A" if mean >= 80 and low == 0 and clashes == 0 else "B" if mean >= 70 and low == 0 and clashes == 0 else "C" if mean >= 50 else "D"
    return {
        "grade": grade,
        "residue_count": len(scores),
        "chain_count": len({chain for chain, _ in residues}),
        "mean_plddt": round(mean, 2),
        "high_confidence_fraction": round(high / len(scores), 3),
        "low_confidence_residues": low,
        "inter_residue_clashes": clashes,
        "bands": {">=90": high, "70-90": sum(70 <= score < 90 for score in scores), "50-70": sum(50 <= score < 70 for score in scores), "<50": low},
        "interpretation": "Local confidence only; use PAE for domain placement and interfaces.",
    }
