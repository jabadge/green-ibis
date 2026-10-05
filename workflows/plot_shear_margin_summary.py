#!/usr/bin/env python3
"""plot_shear_margin_summary.py - map the shear margin softening of a domain.

Reads the plain .mat export written by shear_margin_softening_analysis.m and
maps the damage: the factor rheology_B is divided by in the shear margins
(1 everywhere else). If the export also holds the softening already in the
saved 'Rheology' model (rheology_B before / after), that is mapped next to it
on the same color scale, which stops at the largest damage the method can give.

The map shows the whole domain unless --xlim/--ylim zoom in on a glacier.

Usage:
    python3 plot_shear_margin_summary.py shear_margin_softening_CW.mat
    python3 plot_shear_margin_summary.py shear_margin_softening_CW.mat \\
        --xlim -195000 -140000 --ylim -2295000 -2255000 --output jakobshavn.pdf
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.io import loadmat

# -- Palette (dataviz skill default palette; see references/palette.md) -----
INK_PRIMARY = '#0b0b0b'
INK_MUTED = '#898781'
BASELINE = '#c3c2b7'
SURFACE = '#fcfcfb'
NO_ICE = '#d9d8d2'

# Single-hue sequential ramp (blue, light -> dark): undamaged ice is the
# lightest step, the strongest softening the darkest
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list('seq_blue', [
    '#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
    '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b',
])

PAD_FRAC = 0.03  # margin around the ice-covered extent of the whole-domain view


def make_triangulation(data):
    x = data['x'].ravel()
    y = data['y'].ravel()
    elements = data['elements'].astype(int) - 1  # MATLAB is 1-based
    triang = tri.Triangulation(x, y, elements)

    no_ice = data['ice_levelset'].ravel()[elements] > 0
    triang.set_mask(np.all(no_ice, axis=1))
    return triang


def margin_info(data):
    """The margin_info struct from shear_margin_damage, as a dict of floats."""
    info = data['margin_info'][0, 0]
    return {name: float(np.asarray(info[name]).ravel()[0]) for name in info.dtype.names}


def ice_extent(data):
    x = data['x'].ravel()
    y = data['y'].ravel()
    ice = data['ice_levelset'].ravel() <= 0
    if not ice.any():
        ice = np.ones_like(x, dtype=bool)
    xpad = PAD_FRAC * (x[ice].max() - x[ice].min())
    ypad = PAD_FRAC * (y[ice].max() - y[ice].min())
    return (x[ice].min() - xpad, x[ice].max() + xpad), (y[ice].min() - ypad, y[ice].max() + ypad)


def plot_damage_map(fig, ax, triang, values, title, xlim, ylim, vmax):
    # Rasterized: the margins are only a few elements wide, so they are drawn
    # at full resolution without a file of ~1e5 vector triangles
    tpc = ax.tripcolor(triang, values, shading='gouraud', cmap=SEQUENTIAL_CMAP, vmin=1, vmax=vmax, rasterized=True)
    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title, fontsize=11, color=INK_PRIMARY)
    ax.set_xlabel('x (m)', fontsize=8, color=INK_MUTED)
    ax.set_ylabel('y (m)', fontsize=8, color=INK_MUTED)
    ax.tick_params(colors=INK_MUTED, labelsize=7)
    ax.ticklabel_format(style='plain', useOffset=False)
    for spine in ax.spines.values():
        spine.set_color(BASELINE)
    ax.set_facecolor(NO_ICE)
    # The scale stops at the largest damage the method can give; anything
    # softened more than that in a saved model saturates (arrow on the bar)
    cbar = fig.colorbar(tpc, ax=ax, fraction=0.046, pad=0.04,
                        extend='max' if np.nanmax(values) > vmax * 1.001 else 'neither')
    cbar.set_label('rheology B divided by', fontsize=8, color=INK_MUTED)
    cbar.ax.tick_params(labelsize=7, colors=INK_MUTED)


def build_map(mat_path, out_path, xlim=None, ylim=None, dpi=400):
    data = loadmat(mat_path)
    domain = str(np.asarray(data['domain']).ravel()[0]) if 'domain' in data else ''
    info = margin_info(data)
    triang = make_triangulation(data)

    damage = data['damage'].ravel()
    damage_applied = data['damage_applied'].ravel() if 'damage_applied' in data else np.array([])

    full_xlim, full_ylim = ice_extent(data)
    xlim = xlim or full_xlim
    ylim = ylim or full_ylim

    panels = [(damage, 'Damage from the lateral shear\n{} vertices softened, enhancement on A up to {:g}'.format(
        int((damage > 1).sum()), info['max_enhancement']))]
    vmax = info['max_damage']
    if damage_applied.size == damage.size:
        panels.append((damage_applied, 'Softening in the saved model\nrheology B before / after'))

    aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
    panel_w = 8.0
    panel_h = float(np.clip(panel_w * aspect, 3.0, 14.0))
    fig, axes = plt.subplots(1, len(panels), figsize=(panel_w * len(panels) + 1.0, panel_h + 1.5), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, (values, title) in zip(axes[0], panels):
        plot_damage_map(fig, ax, triang, values, title, xlim, ylim, vmax)

    fig.suptitle('Shear margin softening{}'.format(' — ' + domain if domain else ''),
                 fontsize=14, color=INK_PRIMARY)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)
    print('Wrote {}'.format(out_path))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('export', help='.mat export from shear_margin_softening_analysis.m')
    parser.add_argument('--output', help='Figure to write (default: the export name with _map.pdf)')
    parser.add_argument('--xlim', nargs=2, type=float, metavar=('XMIN', 'XMAX'), help='Zoom: x range (m)')
    parser.add_argument('--ylim', nargs=2, type=float, metavar=('YMIN', 'YMAX'), help='Zoom: y range (m)')
    parser.add_argument('--dpi', type=int, default=400, help='Resolution of the rasterized map (default 400)')
    args = parser.parse_args()

    if not os.path.isfile(args.export):
        raise SystemExit(
            "Export does not exist: {}\n"
            "Run shear_margin_softening_analysis.m first.".format(args.export))

    out_path = args.output or os.path.splitext(args.export)[0] + '_map.pdf'
    build_map(args.export, out_path,
              xlim=tuple(args.xlim) if args.xlim else None,
              ylim=tuple(args.ylim) if args.ylim else None,
              dpi=args.dpi)


if __name__ == '__main__':
    main()
