# structuregrade

**A mean pLDDT of 72 sounds fine. Two residues stuck at 30 is the part that actually breaks your downstream pipeline.**

![structuregrade cover](demo/cover.png)

AlphaFold writes per-residue confidence into the B-factor column, and almost every pipeline collapses it to one mean number before deciding whether to trust the structure. That throws away exactly the information you need: a model can average 72 pLDDT with every residue reasonably confident, or average 72 with a beautifully modeled core and two residues the model was essentially guessing at. Those are not the same structure. structuregrade reads the actual per-residue distribution and reports a grade, a clash count, and the low-confidence tail the mean was hiding.

![structuregrade workbench](demo/dashboard.png)

## How it works

Parse PDB ATOM/HETATM records, pull pLDDT straight out of the B-factor field per the AlphaFold convention, bucket residues into confidence bands, and flag anything under 50 by name instead of averaging it away. A separate geometry pass checks inter-residue distances for physically impossible clashes, because a locally confident residue can still be placed somewhere it can't actually be.

## What ships

- PDB ATOM/HETATM parser that reads AlphaFold's pLDDT-in-B-factor convention directly
- Confidence banding (>=90 / 70-90 / 50-70 / <50) plus a letter grade, not just a raw mean
- Inter-residue clash detection independent of the confidence score
- CLI, JSON API, browser workbench, Docker, tests, and CI

## Run it end to end

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -e .
structuregrade demo
structuregrade serve
```

Open <http://127.0.0.1:8090>. Analyze your own JSON input with `structuregrade analyze input.json`.

## API

- `GET /api/demo` returns the committed fixture and result.
- `POST /api/analyze` runs the same engine on a JSON body.

## The result

The eight-residue fixture averages 71.88 pLDDT — a number that would sail through most "pLDDT > 70" filters. structuregrade assigns it a C instead, because two of those eight residues are below 50 and the mean was quietly burying them. Zero geometry clashes, for what that's worth on a structure this uninterpretable in places.

## Update: fixed the clash detector flagging every normal peptide bond

The clash pass compared every pair of atoms from different residues
against a 1.5 Å threshold, excluding only atom pairs from the *same*
residue. It never excluded the one bond that legitimately spans a residue
boundary: the C(i)–N(i+1) peptide bond connecting adjacent residues along
the backbone, which sits at roughly 1.33 Å — well under the "clash"
threshold, and completely normal chemistry. On any real, full-atom
AlphaFold structure (N/CA/C/O plus sidechains per residue, not the
single-CA-per-residue trace the bundled demo happens to use), every
adjacent residue pair in the chain would be reported as a clash, silently
downgrading grades on structurally perfect proteins.

Verified directly: a two-residue backbone built with standard, textbook
bond geometry (including the ~1.33 Å peptide bond) was flagged with 1
"clash" and downgraded from what should be an A to a C, despite nothing
being structurally wrong. Fixed by excluding the canonical C(i)–N(i+1)
peptide-bond pair between sequentially adjacent residues in the same
chain from clash counting, while leaving genuine non-bonded steric
clashes between any other atom pair fully detected.
`tests/test_peptide_bond_clash.py` covers the original false positive,
confirms a genuine non-bonded clash between unrelated residues is still
caught, and confirms the bundled CA-only demo fixture (which never
exercised this bond in the first place) is unaffected.

## Scope

pLDDT is local per-residue confidence, not interface confidence or domain-placement confidence — PAE and domain-aware review still matter for multi-domain structures and complexes. This is a triage gate that catches "the mean is lying to you," not a replacement for structural validation.

## Test

```bash
python -m unittest discover -s tests -v
```

## Research basis

- [AlphaFold confidence guidance](https://alphafold.ebi.ac.uk/faq)

MIT licensed.
