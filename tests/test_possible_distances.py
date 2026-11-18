import csv
import importlib.util
import io
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location('possible_distances',
    Path(__file__).parents[1] / 'src' / 'possible_distances.py')
distance = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = distance
spec.loader.exec_module(distance)


class Atom:
    def __init__(self, residue, name, xyz):
        self.residue, self.name, self.scene_coord = residue, name, np.array(xyz, dtype=float)

    def string(self, **kwargs):
        r = self.residue
        return '#%s/%s:%s%s@%s' % (r.structure.id_string, r.chain_id, r.number,
                                  r.insertion_code, self.name)


class Residue:
    def __init__(self, model, chain, number, xyz, name='LYS', insertion=''):
        self.structure, self.chain, self.chain_id = model, chain, chain.chain_id
        self.number, self.name, self.insertion_code = number, name, insertion
        self.atoms = {'CA': Atom(self, 'CA', xyz)}

    def find_atom(self, name):
        return self.atoms.get(name)


def bond(r1, r2):
    return NS(atoms=[r1.find_atom('CA'), r2.find_atom('CA')])


class DistanceTests(unittest.TestCase):
    def setUp(self):
        self.model = NS(id_string='7', name='Test')
        a, b = NS(chain_id='A'), NS(chain_id='B')
        self.r = [Residue(self.model, a, 10, (0, 0, 0)),
                  Residue(self.model, a, 10, (3, 0, 0), insertion='A'),
                  Residue(self.model, b, 20, (0, 4, 0))]
        self.model.residues = self.r

    def test_all_pairs_exact_distances_and_table(self):
        table, pb = io.StringIO(), io.StringIO()
        observed = [bond(self.r[0], self.r[2]), bond(self.r[2], self.r[0])]
        result = distance.calculate_distances(self.model, observed, tsv=table, pb_file=pb)
        self.assertEqual(result.possible_pairs, 3)  # n(n-1)/2, without self or reverse pairs
        self.assertEqual(result.observed_pairs, 1)
        self.assertEqual(result.possible_counts.tolist(), [0, 1, 2])
        self.assertEqual(result.observed_counts.tolist(), [0, 0, 1])
        rows = list(csv.DictReader(io.StringIO(table.getvalue()), delimiter='\t'))
        self.assertEqual([float(r['distance_A']) for r in rows], [3., 4., 5.])
        self.assertEqual([r['observed_XL'] for r in rows], ['0', '1', '0'])
        self.assertEqual(rows[0]['insertion_2'], 'A')
        self.assertEqual(len(pb.getvalue().splitlines()), 3)
        self.assertIn('#7/A:10A@CA', pb.getvalue())

    def test_chain_scopes(self):
        self.assertEqual(distance.calculate_distances(self.model, scope='intra').possible_pairs, 1)
        self.assertEqual(distance.calculate_distances(self.model, scope='inter').possible_pairs, 2)
        self.assertEqual(distance.calculate_distances(self.model,
            [bond(self.r[0], self.r[2])], scope='intra').observed_pairs, 0)

    def test_missing_atoms_and_non_lysine(self):
        self.r[1].name = 'ALA'
        result = distance.calculate_distances(self.model)
        self.assertEqual(result.possible_pairs, 1)
        self.r[1].name = 'LYS'
        self.r[1].atoms.clear()
        result = distance.calculate_distances(self.model)
        self.assertEqual(result.missing_atoms, 1)
        self.assertEqual(result.possible_pairs, 1)
        with self.assertRaises(ValueError):
            distance.calculate_distances(self.model, atom_name='NZ')

    def test_observed_uses_same_atom_as_background(self):
        self.r[0].atoms['NZ'] = Atom(self.r[0], 'NZ', (0, 0, 0))
        self.r[1].atoms['NZ'] = Atom(self.r[1], 'NZ', (10, 0, 0))
        self.r[2].atoms['NZ'] = Atom(self.r[2], 'NZ', (20, 0, 0))
        table = io.StringIO()
        result = distance.calculate_distances(self.model, [bond(self.r[0], self.r[1])],
                                              atom_name='NZ', tsv=table)
        observed = [r for r in csv.DictReader(io.StringIO(table.getvalue()), delimiter='\t')
                    if r['observed_XL'] == '1']
        self.assertEqual(float(observed[0]['distance_A']), 10.)
        self.assertEqual(result.observed_counts.sum(), 1)

    def test_other_model_and_self_links_excluded(self):
        other = Residue(NS(id_string='8'), self.r[0].chain, 10, (0, 0, 0))
        result = distance.calculate_distances(self.model,
            [bond(self.r[0], other), bond(self.r[0], self.r[0])])
        self.assertEqual(result.ignored_observed, 2)
        self.assertEqual(result.observed_pairs, 0)

    def test_cancellation_and_validation(self):
        with self.assertRaises(distance.DistanceCalculationCancelled):
            distance.calculate_distances(self.model, progress=lambda *_: False)
        for width in (0, -1, float('nan'), float('inf'), 1e-12):
            with self.subTest(width=width), self.assertRaises(ValueError):
                distance.calculate_distances(self.model, bin_width=width)
        for options in ({'atom_name': 'CB'}, {'scope': 'between'}):
            with self.assertRaises(ValueError):
                distance.calculate_distances(self.model, **options)

    def test_degenerate_geometry_and_no_scope_pairs(self):
        self.r[1].atoms['CA'].scene_coord[:] = 0
        self.r[2].atoms['CA'].scene_coord[:] = 0
        self.assertEqual(distance.calculate_distances(self.model).possible_pairs, 3)
        self.r[2].chain = self.r[0].chain
        with self.assertRaises(ValueError):
            distance.calculate_distances(self.model, scope='inter')

    def test_plot_normalization_and_save(self):
        os.environ['MPLBACKEND'] = 'Agg'
        import matplotlib.pyplot as plt
        result = distance.calculate_distances(self.model, [bond(self.r[0], self.r[1])])
        fig = distance.plot_distribution(result, 'Test')
        self.assertEqual(fig.axes[0].get_ylabel(), 'Fraction of pairs per bin')
        for patch in fig.axes[0].patches:
            self.assertAlmostEqual(patch.get_data().values.sum(), 1.)
        with TemporaryDirectory() as folder:
            for ext in ('png', 'svg'):
                path = Path(folder) / ('plot.' + ext)
                fig.savefig(path)
                self.assertGreater(path.stat().st_size, 1000)
        plt.close(fig)


if __name__ == '__main__':
    unittest.main()
