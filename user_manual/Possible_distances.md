# All possible lysine–lysine distances

Implements the request in [upstream issue #6](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle/issues/6):
compare the distribution of all possible K–K distances in a structure with
mapped crosslink distances, and export the underlying distance table.

## Use

1. Load a molecular model in ChimeraX and open XMAS.
2. Choose **Analyze → All possible K–K distances…**. No mapped crosslinks are
   required. The original analysis functions remain under
   **Analyze → Crosslink analysis…**.
3. Select one molecular model. All its chains are included. Choose **CA**
   (default, Cα–Cα) or **NZ** (lysine side-chain Nζ–Nζ), and all pairs,
   intrachain only, or interchain only. Set the histogram bin width in Å.
4. Optionally check crosslink models to compare with the possible pairs.
   Their current check states in XMAS initialize the list. Checked models
   are pooled as unique residue pairs. Uncheck all for background alone.
5. Optionally select `.pb` export. Click **OK**, choose a filename, and the
   calculation writes the table and plots, then displays the plot.
   The calculation can be cancelled from its progress dialog.

The main XMAS toolbar retains **Analyze**, **Export**, **Visualize**, and
**Integrate**; this feature is an action under Analyze, not an extra toolbar
button. AlphaLink2 likewise uses the existing Export subset action.

## Outputs

| File | Contents |
| --- | --- |
| `name.tsv` | Every unique eligible K–K pair, distance in Å, and an `observed_XL` flag (0/1) |
| `name.png` / `name.svg` | Comparison histogram using the same bins and independently normalized populations |
| `name.pb` (optional) | All eligible atom pairs as a ChimeraX pseudobond file |

The TSV records model ID, both original chain IDs, residue numbers, insertion
codes, atom choice, distance, and observed status. It is not restricted to
observed crosslinks. The optional PB file can be opened in ChimeraX; it is
written without automatically adding potentially millions of lines to the
viewer. Output size and runtime grow quadratically with the lysine count.
Pairs are processed one row of the distance matrix at a time, without storing
a complete distance matrix or all output rows in memory.

## Interpretation

- Each unordered pair of distinct resolved `LYS` residues is counted once:
  for N eligible lysines the all-pairs count is N(N−1)/2. There are no
  self-pairs or reversed duplicates. Modified lysines with other residue names
  are not included.
- Lysines missing the requested atom are skipped and counted in the log.
  NZ measurements do not silently fall back to CA. Unresolved residues without
  coordinates cannot contribute a distance.
- Observed crosslinks must map to two eligible lysines in this same model and
  satisfy the selected chain scope. Replicate/reversed/atom-level duplicates
  collapse to one pair. Other-model, non-K–K, missing-atom, and same-residue
  records are excluded; the log reports the excluded record count.
- Observed distances are recomputed using the chosen CA or NZ atoms, even if
  their original pseudobonds connect different atoms. Both distributions
  therefore use the same distance definition and current coordinates.
- Plot values are fractions per bin, each population summing to one, so a large
  theoretical population does not overwhelm a small measured set. Legend counts
  show the respective sample sizes. The observed set is a subset of the possible
  set; no statistical significance test is implied.
- No distance cutoff, SASA/accessibility filter, digestion model, or reactivity
  prediction is imposed. These are geometric possibilities, not a prediction
  that every pair can form a chemical crosslink. All pairs in the chosen model
  are retained, including long distances.
- Current scene coordinates and the active coordinate set/alternate locations
  are used. Independent structures or conformers are not pooled. Use one intended
  assembly/model; select different models in separate runs.
- The comparison uses the checked crosslink models in this dialog. The Export
  dialog's distance/score/overlap filters and display visibility are not applied.
  To compare a curated subset, first export that subset as `.pb`, open it, and
  check only that crosslink model here. Mapping ambiguity is not resolved
  automatically; each distinct eligible mapped residue pair is counted.

Calculations stage output files before publishing them, so cancellation or
calculation errors do not replace previous files. File-system failures while
replacing the completed outputs can still leave a partially updated output set.

## Validation

Tests cover analytically known 3–4–5 Å pairs, exact N(N−1)/2 counts,
insertion codes, both chain scopes, deduplication, missing atoms, matching
CA/NZ definitions, excluded observations, cancellation, invalid settings,
normalized plots, and PNG/SVG saving.

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

The live ChimeraX interface requires a manual check: open a small known model,
run Analyze with and without mapped links, inspect TSV distances/counts, save
and reopen the optional PB file, and check cancellation. Test AlphaLink2 via
Export as well, selecting its format and using the existing Export subset
button. A live ChimeraX installation is not available in the implementation
environment.
