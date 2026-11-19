"""Run with python -m unittest discover -s tests -v (no ChimeraX required)."""
import ast
import importlib.util
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import unittest

MODULE = Path(__file__).resolve().parents[1] / 'src' / 'alphalink_export.py'
spec = importlib.util.spec_from_file_location('alphalink_export', MODULE)
export = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = export
spec.loader.exec_module(export)


def chain(model, name, numbers, sequence=None):
    c = NS(structure=NS(id_string=model), chain_id=name,
           characters=sequence or 'A' * len(numbers))
    c.residues = [None if n is None else NS(number=n, insertion_code='', chain=c)
                  for n in numbers]
    return c


def bond(a, i, b, j):
    return NS(atoms=[NS(residue=a.residues[i]), NS(residue=b.residues[j])])


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.a = chain('2', 'X', [41, None, 45, 45, 90, 91], 'AKGSTA')
        self.a.residues[3].insertion_code = 'A'
        self.b = chain('1', 'X', [201, 202, 203, 204, 205, 206], 'GASKTA')
        self.pb = bond(self.a, 3, self.b, 2)

    def test_numbering_gaps_insertion_codes_and_fasta_order(self):
        result = export.build_alphalink_export([self.pb], [self.b, self.a], '.01')
        self.assertEqual(result.crosslinks, '3 A 4 B 0.01\n')
        self.assertEqual(result.fasta, '>A model_1 chain_X\nGASKTA\n>B model_2 chain_X\nAKGSTA\n')
        self.assertIn('B\t2\t2\tX\t\t\n', result.mapping)
        self.assertIn('B\t4\t2\tX\t45\tA\n', result.mapping)
        self.assertEqual(self.a.residues[3].number, 45)

    def test_duplicate_reverse_and_intralinks(self):
        reverse = NS(atoms=self.pb.atoms[::-1])
        intra = bond(self.a, 0, self.a, 4)
        result = export.build_alphalink_export([self.pb, reverse, self.pb, intra],
                                              [self.a, self.b], .02)
        self.assertEqual(result.count, 2)
        self.assertEqual(result.crosslinks, '1 A 5 A 0.02\n4 A 3 B 0.02\n')

    def test_stable_distinct_chains_with_same_chain_id(self):
        self.assertEqual([id(c) for c in export.linked_chains([self.pb, self.pb])],
                         [id(self.b), id(self.a)])

    def test_invalid_fdr(self):
        for fdr in ['', None, 'nan', 'inf', -.1, 1.1, '1%']:
            with self.subTest(fdr=fdr), self.assertRaises(ValueError):
                export.build_alphalink_export([self.pb], [self.a, self.b], fdr)

    def test_invalid_mapping_and_sequence(self):
        with self.assertRaises(ValueError):
            export.build_alphalink_export([self.pb], [self.a], .01)
        with self.assertRaises(ValueError):
            export.build_alphalink_export([self.pb], [self.a, self.a], .01)
        self.a.characters = 'AXGSTA'
        with self.assertRaises(ValueError):
            export.build_alphalink_export([self.pb], [self.a, self.b], .01)
        self.a.characters = 'AA'
        with self.assertRaises(ValueError):
            export.build_alphalink_export([self.pb], [self.a, self.b], .01)

    def test_empty_self_and_chain_limit(self):
        for pbs, chains in [([], [self.a]), ([bond(self.a, 0, self.a, 0)], [self.a]),
                            ([self.pb], [chain(str(i), 'A', [1]) for i in range(63)])]:
            with self.assertRaises(ValueError):
                export.build_alphalink_export(pbs, chains, .01)

    def test_upstream_loader(self):
        # Integration harness supplies the exact function fetched from upstream.
        source_path = Path(os.environ.get('ALPHALINK2_DATA_PY', 'AlphaLink2/unifold/colab/data.py'))
        if not source_path.exists():
            self.skipTest('Upstream AlphaLink2 checkout not available')
        import numpy as np
        tree = ast.parse(source_path.read_text())
        func = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'load_crosslinks')
        namespace = {'np': np}
        exec(compile(ast.Module(body=[func], type_ignores=[]), str(source_path), 'exec'), namespace)
        for pbs in [[self.pb], [self.pb, bond(self.a, 0, self.a, 4)]]:
            with TemporaryDirectory() as folder:
                path = Path(folder) / 'links.csv'
                result = export.build_alphalink_export(pbs, [self.a, self.b], .01)
                path.write_text(result.crosslinks)
                links = namespace['load_crosslinks'](str(path))
                self.assertEqual(links['A']['B'], [(3, 2, .01)])
                if len(pbs) > 1:
                    self.assertEqual(links['A']['A'], [(0, 4, .01)])


if __name__ == '__main__':
    unittest.main()
