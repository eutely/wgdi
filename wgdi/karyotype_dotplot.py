"""Plot WGDI karyotype intervals with chromosome-internal collinearity dots."""

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_rgba
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

import wgdi.base as base


class karyotype_dotplot:
    def __init__(self, options):
        # Follow the option names and left/top convention of WGDI dotplot.
        self.position = 'order'
        self.width = 0.6
        self.figsize = '10,6.18'
        self.markersize = 0.5
        self.markeralpha = 0.5
        self.fontsize = 10
        self.dpi = 500
        self.block_length = 1
        self.labels = True
        self.legend = True
        self.strict = True
        self.genome1_name = ''
        self.genome2_name = ''
        self.title = ''
        self.ancestor_top = None
        self.savefig = 'wgdi.karyotype_dotplot.svg'

        required = {'gff1', 'gff2', 'lens1', 'lens2', 'collinearity', 'ancestor_left'}
        allowed = set(vars(self)) | required | {'ancestor'}
        values = dict(options)
        for k, v in values.items():
            if k not in allowed:
                raise ValueError(f'Unknown karyotype_dotplot option: {k}; see wgdi -kd ?')
            setattr(self, str(k), v)
            print(f'{k} = {v}')
        # `ancestor` is the input name used by WGDI's standalone -k module.
        if hasattr(self, 'ancestor'):
            if hasattr(self, 'ancestor_left') and self.ancestor_left != self.ancestor:
                raise ValueError('ancestor and ancestor_left specify different files')
            self.ancestor_left = self.ancestor
        for name in required:
            if not getattr(self, name, None):
                raise ValueError(f'Missing required option: {name}')
        if self.position not in ('order', 'end'):
            raise ValueError('position must be order or end, as in WGDI dotplot')
        for name in ('labels', 'legend', 'strict'):
            value = getattr(self, name)
            if str(value).lower() not in ('true', 'false'):
                raise ValueError(f'{name} must be true or false')
            setattr(self, name, base.str_to_bool(value))
        for name in ('width', 'markersize', 'fontsize', 'dpi', 'markeralpha'):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive')
            setattr(self, name, value)
        if self.width > 1 or self.markeralpha > 1:
            raise ValueError('width and markeralpha must not exceed 1')
        self.block_length = int(self.block_length)
        if self.block_length < 1:
            raise ValueError('block_length must be at least 1')
        self.figsize = [float(k) for k in self.figsize.split(',')]
        if len(self.figsize) != 2 or any(not np.isfinite(k) or k <= 0 for k in self.figsize):
            raise ValueError('figsize must contain two positive finite values')
        if self.ancestor_top == 'none':
            self.ancestor_top = None
        if Path(self.savefig).suffix.lower() not in ('.svg', '.pdf', '.png'):
            raise ValueError('savefig must end in .svg, .pdf or .png')
        for name in required | {'ancestor_top'}:
            filename = getattr(self, name, None)
            if filename and Path(filename).resolve() == Path(self.savefig).resolve():
                raise ValueError('savefig must not overwrite an input file')

    def read_genome(self, gff_file, lens_file):
        """Use the same GFF/lens readers as WGDI dotplot and karyotype_mapping."""
        gff = base.newgff(gff_file)
        gff.index = gff.index.astype(str)
        counts = base.newlens(lens_file, 'order')
        lengths = base.newlens(lens_file, 'end')
        if gff.empty or counts.empty:
            raise ValueError('GFF and lens must not be empty')
        if gff.index.has_duplicates or counts.index.has_duplicates:
            raise ValueError('Gene IDs and lens chromosome IDs must be unique')
        for values in (counts, lengths):
            if not np.isfinite(values).all() or (values <= 0).any() or (values % 1 != 0).any():
                raise ValueError('Lens lengths and gene counts must be positive integers')
        if (gff['order'] < 1).any() or (gff['start'] < 0).any() or (gff['end'] < gff['start']).any():
            raise ValueError('GFF contains invalid coordinates or gene orders')
        if gff.duplicated(['chr', 'order']).any():
            raise ValueError('GFF gene orders must be unique within each chromosome')
        for seqid, group in gff.groupby('chr', sort=False):
            if seqid in counts.index:
                if group['order'].max() > counts[seqid] or group['end'].max() > lengths[seqid]:
                    raise ValueError(f'{lens_file}: lens is smaller than GFF for chromosome {seqid}')
        return gff, counts, lengths

    def read_karyotype(self, filename, gff, counts, lengths):
        """Read native five-column intervals and require complete coverage.

        Orders are inclusive: [1, N] covers the entire chromosome. In end mode,
        internal cuts lie halfway between the flanking gene end coordinates;
        terminal intervals extend to 0 and the lens chromosome length.
        """
        ancestor = base.read_classification(filename)
        ancestor = ancestor[ancestor[0].isin(counts.index)].copy()
        if ancestor.empty:
            raise ValueError(f'{filename}: no selected chromosomes')
        if not ancestor[3].map(is_color_like).all():
            raise ValueError(f'{filename}: invalid Matplotlib color')
        if (ancestor[4] < 0).any() or ancestor[4].max() <= 0:
            raise ValueError(f'{filename}: column five must be nonnegative with a positive maximum')
        # Match karyotype.py and align_dotplot.py: normalize column five to alpha.
        ancestor[4] = ancestor[4] / ancestor[4].max()
        regions = []
        for seqid in counts.index:
            group = ancestor[ancestor[0] == seqid].sort_values(1)
            expected = 1
            for _, row in group.iterrows():
                start, end = int(row[1]), int(row[2])
                if start != expected or end < start or end > counts[seqid]:
                    raise ValueError(f'{filename}: gap, overlap or invalid range on chromosome {seqid}; '
                                     f'expected next start {expected}, got {start}..{end}')
                expected = end + 1
            if expected != counts[seqid] + 1:
                raise ValueError(f'{filename}: incomplete coverage on chromosome {seqid}; '
                                 f'expected orders 1..{counts[seqid]}')
            cuts = [0.0]
            if self.position == 'order':
                cuts.extend(float(value) for value in group[2])
            else:
                ends = gff[gff['chr'] == seqid].set_index('order')['end']
                for end in group[2].iloc[:-1]:
                    if end not in ends.index or end + 1 not in ends.index:
                        raise ValueError(f'{filename}: missing boundary genes at {seqid}:{end}/{end+1}')
                    cuts.append(float(ends[end] + ends[end + 1]) / 2)
                cuts.append(float(lengths[seqid]))
            for i, (_, row) in enumerate(group.iterrows()):
                if cuts[i+1] <= cuts[i]:
                    raise ValueError(f'{filename}: non-increasing physical boundaries on {seqid}')
                regions.append([seqid, cuts[i], cuts[i+1], row[3], row[4]])
        return pd.DataFrame(regions, columns=['chr', 'start', 'end', 'color', 'alpha'])

    def pair_positions(self, gff1, gff2, lens1, lens2):
        """First pair gene is left/extant; second is top/ancestor, as in WGDI."""
        positions = {seqid: i for i, seqid in enumerate(lens1.index)}
        data, missing, outside, kept, small = [], 0, 0, 0, 0
        for _, block, _ in base.read_collinearity(self.collinearity):
            pairs = []
            for row in block:
                fields = [field for field in row if field]
                if not fields or fields[0].startswith('#'):
                    continue
                if len(fields) < 4:
                    raise ValueError(f'{self.collinearity}: expected gene1 order1 gene2 order2 [strand]')
                pairs.append((fields[0], fields[2]))
            if len(pairs) < self.block_length:
                small += len(pairs)
                continue
            kept += 1
            for gene1, gene2 in pairs:
                if gene1 not in gff1.index or gene2 not in gff2.index:
                    missing += 1
                    continue
                left, top = gff1.loc[gene1], gff2.loc[gene2]
                if left['chr'] not in lens1.index or top['chr'] not in lens2.index:
                    outside += 1
                    continue
                fraction = float(top[self.position]) / lens2[top['chr']]
                x = positions[left['chr']] - self.width / 2 + fraction * self.width
                y = float(left[self.position]) / lens1.max()
                data.append([gene1, gene2, x, y])
        if missing:
            message = f'{missing} collinear pairs do not match gff1/gff2 in the supplied direction'
            if self.strict:
                raise ValueError(message)
            warnings.warn(message + '; skipped')
        if not data:
            raise ValueError('No drawable collinear pairs remain')
        self.stats = dict(pairs=len(data), blocks=kept, missing=missing,
                          outside_selection=outside, small_block_pairs=small)
        return pd.DataFrame(data, columns=['gene1', 'gene2', 'x', 'y'])

    def legend_handles(self, regions, top):
        """Use native colors and the ancestor lens order for legend labels."""
        handles, seen = [], set()
        order = regions if top is None else pd.concat([top, regions], ignore_index=True)
        region_colors = regions['color'].map(to_rgba)
        for _, entry in order.iterrows():
            key = to_rgba(entry['color'])
            matches = regions[region_colors == key]
            if key in seen or matches.empty:
                continue
            seen.add(key)
            row = matches.iloc[0]
            label = row['color']
            if top is not None:
                seqids = top[top['color'].map(to_rgba) == key]['chr'].drop_duplicates()
                if not seqids.empty:
                    label = ', '.join(seqids)
            handles.append(Patch(facecolor=row['color'], alpha=row['alpha'], label=label))
        return handles

    def plot(self, lens, regions, points, top):
        with plt.rc_context({'svg.fonttype': 'none', 'font.size': self.fontsize}):
            fig, ax = plt.subplots(figsize=self.figsize)
            try:
                positions = {seqid: i for i, seqid in enumerate(lens.index)}
                for _, row in regions.iterrows():
                    loc = [positions[row['chr']] - self.width / 2, row['start'] / lens.max()]
                    height = (row['end'] - row['start']) / lens.max()
                    # Same rectangle helper and alpha convention as WGDI -k.
                    base.Rectangle(ax, loc, height, self.width, row['color'], row['alpha'])
                    ax.patches[-1].set_antialiased(False)
                scatter = ax.scatter(points['x'], points['y'], s=self.markersize,
                                     c='black', alpha=self.markeralpha,
                                     edgecolors=None, linewidths=0, marker='o', zorder=3)
                scatter.set_gid('synteny-points')
                for seqid, length in lens.items():
                    ax.add_patch(Rectangle((positions[seqid] - self.width / 2, 0),
                                           self.width, length / lens.max(), fill=False,
                                           edgecolor='#444444', linewidth=0.5, zorder=4))
                for spine in ax.spines.values():
                    spine.set_visible(False)
                ax.set_yticks([])
                ax.set_xlim(-0.65, len(lens) - 0.35)
                ax.set_ylim(-0.015, 1.025)
                ax.tick_params(axis='x', length=0, labelsize=self.fontsize)
                if self.labels:
                    ax.set_xticks(range(len(lens)))
                    ax.set_xticklabels(lens.index)
                    ax.set_xlabel(self.genome1_name, labelpad=10)
                else:
                    ax.set_xticks([])
                if self.title:
                    ax.set_title(self.title, fontsize=self.fontsize + 2, pad=14)
                if self.legend:
                    ax.legend(handles=self.legend_handles(regions, top),
                              title=self.genome2_name or 'Ancestral karyotype', frameon=False,
                              loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0,
                              handlelength=1, handletextpad=0.5)
                fig.subplots_adjust(left=0.035, right=0.84 if self.legend else 0.98,
                                    bottom=0.12 if self.labels else 0.04,
                                    top=0.91 if self.title else 0.98)
                fig.savefig(self.savefig, dpi=self.dpi, bbox_inches='tight', facecolor='white')
            finally:
                plt.close(fig)

    def run(self):
        gff1, counts1, lengths1 = self.read_genome(self.gff1, self.lens1)
        gff2, counts2, lengths2 = self.read_genome(self.gff2, self.lens2)
        lens1 = counts1 if self.position == 'order' else lengths1
        lens2 = counts2 if self.position == 'order' else lengths2
        self.regions = self.read_karyotype(self.ancestor_left, gff1, counts1, lengths1)
        top = None
        if self.ancestor_top:
            top = self.read_karyotype(self.ancestor_top, gff2, counts2, lengths2)
        self.points = self.pair_positions(gff1, gff2, lens1, lens2)
        self.plot(lens1, self.regions, self.points, top)
        self.stats['regions'] = len(self.regions)
        self.stats['chromosomes'] = len(lens1)
        print(f'Wrote {len(self.points)} collinear pairs and {len(self.regions)} karyotype '
              f'regions to {self.savefig}; all {len(lens1)} displayed chromosomes are covered')
        return self.stats
