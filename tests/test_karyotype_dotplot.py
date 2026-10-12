"""Regression tests for the native WGDI inputs and karyotype coverage."""
import contextlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np

from wgdi.karyotype_dotplot import karyotype_dotplot


class KaryotypeDotplotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.gff1 = self.write('species.gff', 's1\ts1a\t1\t10\t+\t1\ns1\ts1b\t21\t30\t+\t2\ns1\ts1c\t41\t50\t+\t3\ns2\ts2a\t1\t10\t+\t1\n')
        self.gff2 = self.write('ancestor.gff', 'a1\ta1a\t1\t10\t+\t1\na1\ta1b\t21\t30\t+\t2\na1\ta1c\t41\t50\t+\t3\n')
        self.lens1 = self.write('species.lens', 's1\t100\t3\ns2\t50\t1\n')
        self.lens2 = self.write('ancestor.lens', 'a1\t60\t3\n')
        self.left = self.write('species.ancestor.txt', 's1\t1\t2\tred\t1\ns1\t3\t3\tblue\t2\ns2\t1\t1\tred\t2\n')
        self.top = self.write('ancestor.ancestor.txt', 'a1\t1\t3\tred\t1\n')
        self.col = self.write('pairs.collinearity', '# Alignment 1: score=100 pvalue=0.01 N=3 s1&a1 minus\ns1a 1 a1c 3 -1\ns1b 2 a1b 2 -1\ns1c 3 a1a 1 -1\n')
        self.options = dict(gff1=self.gff1, gff2=self.gff2, lens1=self.lens1,
                            lens2=self.lens2, ancestor_left=self.left,
                            ancestor_top=self.top, collinearity=self.col,
                            savefig=str(self.path/'plot.svg'), dpi='72')

    def write(self, name, content):
        path = self.path/name
        path.write_text(content, encoding='utf-8')
        return str(path)

    def make_plot(self, **options):
        with contextlib.redirect_stdout(io.StringIO()):
            plot = karyotype_dotplot(dict(self.options, **options).items())
        return plot

    def run_plot(self, **options):
        plot = self.make_plot(**options)
        with contextlib.redirect_stdout(io.StringIO()):
            plot.run()
        return plot

    def test_native_order_coordinates_and_complete_regions(self):
        plot = self.run_plot()
        self.assertEqual(plot.stats['pairs'], 3)
        self.assertEqual(plot.stats['regions'], 3)
        self.assertEqual(plot.regions['color'].tolist(), ['red', 'blue', 'red'])
        self.assertEqual(plot.regions['alpha'].tolist(), [0.5, 1, 1])
        self.assertEqual(plot.regions[['start','end']].values.tolist(), [[0,2],[2,3],[0,1]])
        np.testing.assert_allclose(plot.points['x'], [0.3, 0.1, -0.1])
        np.testing.assert_allclose(plot.points['y'], [1/3, 2/3, 1])
        # s2 has no anchors, but its complete karyotype is still colored.
        self.assertEqual(plot.regions[plot.regions['chr']=='s2']['end'].iloc[0], 1)
        root = ET.parse(plot.savefig).getroot()
        points = root.find('.//*[@id="synteny-points"]')
        self.assertEqual(len(points.findall('.//{*}use')), 3)
        self.assertFalse(plt.get_fignums())

    def test_end_coordinates_and_gapless_physical_intervals(self):
        plot = self.run_plot(position='end')
        self.assertEqual(plot.regions[['start','end']].values.tolist(), [[0,40],[40,100],[0,50]])
        np.testing.assert_allclose(plot.points['x'], [0.2,0,-0.2])
        np.testing.assert_allclose(plot.points['y'], [0.1,0.3,0.5])
        # Native WGDI uses gene END positions, not interval midpoints.
        self.assertEqual(plot.points.iloc[0]['gene2'], 'a1c')

    def test_colors_independent_of_collinearity_density(self):
        full = self.run_plot()
        sparse = self.write('sparse.collinearity', '# Alignment 1: score=100 pvalue=0.01 N=1 s1&a1 plus\ns1b 2 a1b 2 1\n')
        few = self.run_plot(collinearity=sparse)
        self.assertTrue(full.regions.equals(few.regions))
        self.assertEqual(few.stats['pairs'], 1)

    def test_coverage_errors(self):
        fixtures = {
            'gap': 's1\t1\t1\tred\t1\ns1\t3\t3\tblue\t1\ns2\t1\t1\tred\t1\n',
            'overlap': 's1\t1\t2\tred\t1\ns1\t2\t3\tblue\t1\ns2\t1\t1\tred\t1\n',
            'missing': 's1\t1\t3\tred\t1\n',
            'overflow': 's1\t1\t4\tred\t1\ns2\t1\t1\tred\t1\n',
        }
        for name, content in fixtures.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.run_plot(ancestor_left=self.write(name+'.txt', content))

    def test_gene_matching_direction_and_selection(self):
        reverse = self.write('reverse.collinearity', '# Alignment 1: score=100 pvalue=0.01 N=1 a1&s1 plus\na1a 1 s1a 1 1\n')
        with self.assertRaisesRegex(ValueError, 'supplied direction'):
            self.run_plot(collinearity=reverse)
        missing = self.write('missing.collinearity', '# Alignment 1: score=100 pvalue=0.01 N=2 s1&a1 plus\ns1a 1 a1a 1 1\nunknown 2 a1b 2 1\n')
        with self.assertWarnsRegex(UserWarning, '1 collinear pairs'):
            plot = self.run_plot(collinearity=missing, strict='false')
        self.assertEqual(plot.stats['missing'], 1)
        selected = self.write('selected.lens', 's1\t100\t3\n')
        plot = self.run_plot(lens1=selected)
        self.assertEqual(plot.stats['chromosomes'], 1)
        self.assertEqual(plot.stats['regions'], 2)

    def test_formats_alias_and_validation(self):
        for suffix, signature in (('.pdf', b'%PDF'), ('.png', b'\x89PNG')):
            plot = self.run_plot(savefig=str(self.path/('plot'+suffix)))
            self.assertTrue(Path(plot.savefig).read_bytes().startswith(signature))
        values = dict(self.options)
        values['ancestor'] = values.pop('ancestor_left')
        with contextlib.redirect_stdout(io.StringIO()):
            plot = karyotype_dotplot(values.items())
            plot.run()
        self.assertEqual(plot.stats['regions'], 3)
        for invalid in ({'position':'bp'}, {'width':'nan'}, {'markersize':'0'},
                        {'figsize':'4,0'}, {'labels':'maybe'}, {'mode':'rank'},
                        {'ancestor_left':self.top}, {'savefig':self.gff1}):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.run_plot(**invalid)

    def test_missing_boundary_gene_and_duplicate_gene(self):
        sparse = self.write('sparse.gff', 's1\ts1a\t1\t10\t+\t1\ns1\ts1c\t41\t50\t+\t3\n')
        with self.assertRaisesRegex(ValueError, 'missing boundary'):
            self.run_plot(gff1=sparse, position='end')
        duplicate = self.write('duplicate.gff', Path(self.gff1).read_text()+'s1\ts1a\t1\t10\t+\t1\n')
        with self.assertRaisesRegex(ValueError, 'unique'):
            self.run_plot(gff1=duplicate)

    def test_cli_template_and_dispatch(self):
        repo = Path(__file__).resolve().parents[1]
        env = dict(os.environ, MPLBACKEND='Agg')
        command = [sys.executable, '-m', 'wgdi.run', '-kd']
        result = subprocess.run(command+['?'], cwd=repo, env=env, text=True, capture_output=True, check=True)
        self.assertIn('ancestor_left =', result.stdout)
        conf = self.write('plot.conf', '[karyotype_dotplot]\n'+'\n'.join(f'{k} = {v}' for k,v in self.options.items()))
        subprocess.run(command+[conf], cwd=repo, env=env, capture_output=True, check=True)
        self.assertTrue(Path(self.options['savefig']).is_file())


if __name__ == '__main__':
    unittest.main()
