"""All unique structure-derived lysine pairs and observed XL comparison.

No ChimeraX/Qt imports: the adapter uses residue/atom objects supplied by the UI.
Only O(number of lysines + histogram bins + observed links) memory is used.
"""
import csv
from dataclasses import dataclass
import math
import numpy as np


class DistanceCalculationCancelled(Exception):
    pass


@dataclass
class DistanceDistribution:
    edges: np.ndarray
    possible_counts: np.ndarray
    observed_counts: np.ndarray
    lysines: int
    missing_atoms: int
    observed_pairs: int
    ignored_observed: int
    atom_name: str
    scope: str

    @property
    def possible_pairs(self):
        return int(self.possible_counts.sum())


def calculate_distances(model, pseudobonds=(), atom_name='CA', scope='all',
                        bin_width=2.0, tsv=None, pb_file=None, progress=None):
    """Enumerate unordered distinct LYS pairs; compare unique mapped pairs.

    Coordinates are scene coordinates in Angstrom. Recompute observed distances
    with the SAME atom choice as the background, independent of PB endpoints.
    Missing requested atoms are skipped (no NZ-to-CA fallback). No distance
    cutoff is applied. progress(completed_rows, total_rows) may return False
    to cancel; neither the structure nor pseudobonds are modified.
    """
    if atom_name not in ('CA', 'NZ'):
        raise ValueError('Choose CA or NZ distances.')
    if scope not in ('all', 'intra', 'inter'):
        raise ValueError('Choose all, intra, or interchain pairs.')
    bin_width = float(bin_width)
    if not math.isfinite(bin_width) or bin_width <= 0:
        raise ValueError('Histogram bin width must be positive and finite.')
    residues, atoms, missing = [], [], 0
    for residue in model.residues:
        if residue is None or residue.name != 'LYS' or residue.chain is None:
            continue
        atom = residue.find_atom(atom_name)
        if atom is None:
            missing += 1
            continue
        residues.append(residue)
        atoms.append(atom)
    if len(atoms) < 2:
        raise ValueError('The model needs at least two lysines with %s atoms.' % atom_name)
    coords = np.asarray([a.scene_coord for a in atoms], dtype=float)
    if coords.shape != (len(atoms), 3) or not np.isfinite(coords).all():
        raise ValueError('The model contains invalid coordinates.')
    # The bounding-box diagonal bounds every pair, including moved models.
    extent = float(np.linalg.norm(np.ptp(coords, axis=0)))
    bin_count = max(1, int(math.ceil(extent / bin_width)))
    if bin_count > 1000000:
        raise ValueError('Bin width is too small for this model; choose a larger value.')
    edges = np.arange(bin_count + 1, dtype=float) * bin_width
    possible_counts = np.zeros(bin_count, dtype=np.int64)
    observed_counts = np.zeros(bin_count, dtype=np.int64)
    indices = {id(r): i for i, r in enumerate(residues)}
    chains = np.array([id(r.chain) for r in residues], dtype=object)
    observed = set()
    ignored = 0
    for pb in pseudobonds:
        endpoints = [indices.get(id(a.residue)) for a in pb.atoms]
        if len(endpoints) != 2 or None in endpoints or endpoints[0] == endpoints[1]:
            ignored += 1
            continue
        i, j = sorted(endpoints)
        same_chain = chains[i] == chains[j]
        if (scope == 'intra' and not same_chain) or (scope == 'inter' and same_chain):
            ignored += 1
            continue
        observed.add((i, j))
    observed_by_start = {}
    for i, j in observed:
        observed_by_start.setdefault(i, set()).add(j)
    writer = csv.writer(tsv, delimiter='\t', lineterminator='\n') if tsv is not None else None
    if writer:
        writer.writerow(['model_id', 'chain_1', 'residue_1', 'insertion_1',
                         'chain_2', 'residue_2', 'insertion_2', 'atom',
                         'distance_A', 'observed_XL'])
    specs = ([a.string(style='command line', omit_structure=False) for a in atoms]
             if pb_file is not None else None)
    for i in range(len(atoms) - 1):
        if progress is not None and progress(i, len(atoms) - 1) is False:
            raise DistanceCalculationCancelled()
        js = np.arange(i + 1, len(atoms))
        same = chains[js] == chains[i]
        if scope == 'intra':
            js = js[same]
        elif scope == 'inter':
            js = js[~same]
        distances = np.linalg.norm(coords[js] - coords[i], axis=1)
        possible_counts += np.histogram(distances, bins=edges)[0]
        observed_js = observed_by_start.get(i, set())
        mask = np.array([int(j) in observed_js for j in js], dtype=bool)
        observed_counts += np.histogram(distances[mask], bins=edges)[0]
        if writer is not None or pb_file is not None:
            r1 = residues[i]
            for j, distance, is_observed in zip(js, distances, mask):
                r2 = residues[j]
                if writer:
                    writer.writerow([model.id_string, r1.chain_id, r1.number,
                                     getattr(r1, 'insertion_code', '').strip(),
                                     r2.chain_id, r2.number,
                                     getattr(r2, 'insertion_code', '').strip(),
                                     atom_name, format(float(distance), '.10g'),
                                     int(is_observed)])
                if pb_file is not None:
                    pb_file.write('%s %s\n' % (specs[i], specs[j]))
    if progress is not None and progress(len(atoms) - 1, len(atoms) - 1) is False:
        raise DistanceCalculationCancelled()
    if not possible_counts.sum():
        raise ValueError('No lysine pairs match the selected chain scope.')
    return DistanceDistribution(edges, possible_counts, observed_counts,
                                len(atoms), missing, len(observed), ignored,
                                atom_name, scope)


def plot_distribution(result, model_name):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    total = result.possible_pairs
    ax.stairs(result.possible_counts / total, result.edges, fill=True,
              alpha=.35, label='All possible K–K (n=%s)' % total)
    if result.observed_pairs:
        ax.stairs(result.observed_counts / result.observed_pairs, result.edges,
                  linewidth=2, label='Mapped K–K XLs (n=%s)' % result.observed_pairs)
    ax.set(xlabel='%s–%s distance (Å)' % (result.atom_name, result.atom_name),
           ylabel='Fraction of pairs per bin',
           title='%s: K–K distances (%s)' % (model_name, result.scope))
    ax.legend()
    fig.tight_layout()
    return fig
