# Export mapped crosslinks to AlphaLink2

Requires XMAS 1.1.4. In ChimeraX, install this source checkout using
`devel install /path/to/ChimeraX_XMAS_bundle_v2`, then restart ChimeraX.
Open **Tools → Structure Analysis → XMAS**.

1. Import and map your XL-MS data as usual. Select the crosslink models and
   open their **Export** dialog.
2. Set the molecular-model, link-type, distance, score, and overlap filters.
   **Export AlphaLink2…** uses those same filters. The `.pb` and DisVis
   checkboxes apply only to the original **Export subset** button.
3. Click **Export AlphaLink2…**. Drag the listed chains into the desired
   FASTA order: first is A, second B, and so on. Distinct chains from different
   molecular models remain separate even when their original chain IDs match.
4. Enter the FDR probability justified by your filtered dataset, for example
   `0.01` for 1%. This is applied to every exported restraint. It is not a
   search-engine score or an automatically estimated per-link FDR. XMAS does
   not infer FDR from scores or the default `PeptidePair.QValue` field.
5. Save the file. Three files are created:

| File | Purpose |
| --- | --- |
| `name.csv` | Headerless, space-separated AlphaLink2 restraints |
| `name.fasta` | Matching chain sequences in the chosen order |
| `name.chains.tsv` | AlphaLink2 positions mapped to model, original chain, residue number, and insertion code |

Upload `name.fasta` in the notebook's **Load FASTA file** cell and `name.csv`
in **Load crosslinks**. The tab-separated mapping table is for your records;
do not upload it as restraints.

## Format and residue numbering

The CSV extension follows the notebook's convention, but fields are separated
by spaces, not commas. There is no header. For example:

```text
4 A 17 B 0.01
8 A 22 A 0.01
```

The five fields are `residueFrom chain1 residueTo chain2 FDR`. Positions are
1-based indices in the exported sequence (`chain.characters`), **not PDB
residue numbers**. A chain beginning at PDB residue 41 still begins at sequence
position 1. Missing structural residues represented in the chain sequence keep
their sequence positions; insertion codes are handled through residue identity
and recorded in the mapping table. The original structures are not renumbered.
Reversed pairs and duplicates from replicates or multiple atom pairs collapse
to one unique residue pair. Distinct chain copies remain distinct.

Only chains participating in the filtered links are exported. These are the
sequences available in the loaded ChimeraX model; they are not necessarily
full-length database proteins. Always use the paired FASTA. If you need a
precursor, a differently trimmed mature protein, or a different construct,
load/map against that sequence before exporting. Replacing the FASTA later
requires remapping the restraints.

The exporter rejects nonstandard/unknown amino-acid letters, absent chain or
residue mappings, invalid FDR values, and same-residue self-restraints. It
supports AlphaLink2 chain IDs A–Z, a–z, 0–9 (up to 62 chains). The Colab has its
own sequence-length and compute limits. Homomeric or repeated-sequence mapping
ambiguity is not resolved automatically: alternative mapped pairs are distinct
restraints, so select the intended chain copies and mappings first.

No additional distance cutoff is imposed by the exporter. The existing distance
slider still filters links if you set it; exporting every valid link does not
require that links already satisfy a model distance. Select the appropriate
AlphaLink2 weights and distance settings separately.

## Validation

```bash
python -m unittest discover -s tests -v
python -m compileall -q src/alphalink_export.py src/tool.py
```

To also check round trips using AlphaLink2's actual `load_crosslinks` function,
point this environment variable at an AlphaLink2 checkout:

```bash
ALPHALINK2_DATA_PY=/path/to/AlphaLink2/unifold/colab/data.py python -m unittest discover -s tests -v
```

This integration test needs NumPy. It extracts only that function, avoiding
AlphaLink2's GPU/model dependencies, and checks both single- and multiple-row
files. Without an AlphaLink2 checkout, that integration test is skipped.
The live Qt/ChimeraX dialog still requires manual validation in ChimeraX.

Input format verified against:
- [AlphaLink2 notebook](https://github.com/Rappsilber-Laboratory/AlphaLink2/blob/main/notebooks/alphalink2.ipynb)
- [AlphaLink2 loader](https://github.com/Rappsilber-Laboratory/AlphaLink2/blob/main/unifold/colab/data.py)
