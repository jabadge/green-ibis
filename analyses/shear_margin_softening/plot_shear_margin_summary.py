#!/usr/bin/env python3
"""plot_shear_margin_summary.py - build a one-page summary PDF per glacier
catchment from the shear margin softening analysis.

Reads the plain .mat exports written by shear_margin_softening_analysis.m
(one file per catchment: catchment_<id>.mat) and, for each, renders:

  - a scatter of strain rate vs. velocity error slope, colored by whether
    the vertex was softened, with the best-fit line, its equation, and
    goodness-of-fit statistics (R^2, RMSE, n)
  - maps of velocity error, velocity error slope, strain rate, rheology_B
    before softening, rheology_B after softening, the before/after
    rheology_B difference, and the resulting stress balance velocity
    difference

The maps are zoomed to the ice-covered part of the catchment that reaches at
least the FAST_VELOCITY_PERCENTILE-th percentile of that catchment's own
velocities (before and/or after softening), and that view is then clipped to
the catchment's ice-covered extent so it doesn't waste space on masked-out
(no ice) territory. Slower/masked-out ice is still plotted, just not shown,
so nothing is masked or hidden from the data itself.

Usage:
    python3 plot_shear_margin_summary.py
    python3 plot_shear_margin_summary.py --input-dir shear_margin_softening_exports \
        --output-dir shear_margin_softening_summaries --catchments 3 4 5
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
import numpy as np
from scipy.io import loadmat

# -- Palette (dataviz skill default palette; see references/palette.md) -----
INK_PRIMARY = '#0b0b0b'
INK_SECONDARY = '#52514e'
INK_MUTED = '#898781'
GRIDLINE = '#e1e0d9'
BASELINE = '#c3c2b7'
SURFACE = '#fcfcfb'

CAT_UNCHANGED = '#2a78d6'  # categorical slot 1 (blue)
CAT_SOFTENED = '#eb6834'   # categorical slot 2 (orange)

# Maps are zoomed to where the ice reaches at least this percentile of the
# catchment's own velocities (pooling before and after softening together) -
# slower ice is still drawn, just outside the view.
FAST_VELOCITY_PERCENTILE = 80.0
ZOOM_PAD_FRAC = 0.15   # extra margin around the fast-flow bounding box
ZOOM_PAD_MIN_M = 2000.0  # minimum margin, so small fast-flow patches aren't cropped edge-to-edge

# For the velocity error map, shrink the color range to this fraction of the
# full +/-max so mid-range values (especially positive ones) show up better;
# values beyond the shrunk range simply saturate to the end colors.
V_ERROR_COLOR_SCALE = 0.5

# Single-hue sequential ramp (blue, light -> dark) for magnitude fields
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list('seq_blue', [
    '#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
    '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b',
])

# Diverging pair (blue <-> red, neutral gray midpoint) for signed fields
DIVERGING_CMAP = LinearSegmentedColormap.from_list('div_blue_red', [
    '#0d366b', '#3987e5', '#f0efec', '#e34948', '#8a1f1f',
])


def make_triangulation(data):
    x = data['x'].ravel()
    y = data['y'].ravel()
    elements = data['elements'].astype(int) - 1  # MATLAB is 1-based
    triang = tri.Triangulation(x, y, elements)

    ice_levelset = data['ice_levelset'].ravel()
    no_ice = ice_levelset[elements] > 0
    triang.set_mask(np.all(no_ice, axis=1))
    return triang, ice_levelset


def compute_zoom_extent(data):
    """Bounding box around the fast-flowing part of the catchment (ice-covered
    vertices at or above the FAST_VELOCITY_PERCENTILE-th percentile of the
    catchment's own before/after velocities), padded a bit so the fast region
    isn't cropped edge to edge. The padded box is then clipped to the
    ice-covered extent, so the view doesn't waste space on masked-out (no
    ice) territory beyond that. Falls back to the full ice-covered extent if
    nothing qualifies."""
    x = data['x'].ravel()
    y = data['y'].ravel()
    vel_before = data['vel_before'].ravel()
    vel_after = data['vel_after'].ravel()
    ice = data['ice_levelset'].ravel() <= 0

    domain_vel = np.concatenate([vel_before[ice], vel_after[ice]]) if ice.any() else np.concatenate([vel_before, vel_after])
    if domain_vel.size:
        cutoff = np.percentile(domain_vel, FAST_VELOCITY_PERCENTILE)
    else:
        cutoff = 0.0

    fast = ice & ((vel_before >= cutoff) | (vel_after >= cutoff))
    if not fast.any():
        print('   No vertices reach the {:g}th percentile ({:g} m/yr) - showing full ice-covered extent.'.format(
            FAST_VELOCITY_PERCENTILE, cutoff))
        fast = ice if ice.any() else np.ones_like(x, dtype=bool)

    xmin, xmax = x[fast].min(), x[fast].max()
    ymin, ymax = y[fast].min(), y[fast].max()

    xpad = max(ZOOM_PAD_FRAC * (xmax - xmin), ZOOM_PAD_MIN_M)
    ypad = max(ZOOM_PAD_FRAC * (ymax - ymin), ZOOM_PAD_MIN_M)

    xlim = [xmin - xpad, xmax + xpad]
    ylim = [ymin - ypad, ymax + ypad]

    # Don't let padding spill into masked-out (no ice) territory beyond the
    # catchment's own ice-covered extent.
    if ice.any():
        xlim[0] = max(xlim[0], x[ice].min())
        xlim[1] = min(xlim[1], x[ice].max())
        ylim[0] = max(ylim[0], y[ice].min())
        ylim[1] = min(ylim[1], y[ice].max())

    return tuple(xlim), tuple(ylim)


def style_map_axis(ax, title, xlim=None, ylim=None):
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10, color=INK_PRIMARY)
    ax.tick_params(colors=INK_MUTED, labelsize=7)
    for spine in ax.spines.values():
        spine.set_color(BASELINE)
    ax.set_facecolor('#d9d8d2')  # subdued no-ice background
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)


def plot_sequential_map(fig, ax, triang, values, title, xlim, ylim, cmap=SEQUENTIAL_CMAP, vmin=None, vmax=None):
    tpc = ax.tripcolor(triang, values, shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
    style_map_axis(ax, title, xlim, ylim)
    cbar = fig.colorbar(tpc, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=7, colors=INK_MUTED)
    return tpc


def plot_diverging_map(fig, ax, triang, values, title, xlim, ylim, absmax_scale=1.0):
    finite = values[np.isfinite(values)]
    absmax = np.nanmax(np.abs(finite)) if finite.size else 1.0
    absmax = absmax if absmax > 0 else 1.0
    absmax = absmax * absmax_scale
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-absmax, vmax=absmax)
    tpc = ax.tripcolor(triang, values, shading='gouraud', cmap=DIVERGING_CMAP, norm=norm)
    style_map_axis(ax, title, xlim, ylim)
    cbar = fig.colorbar(tpc, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=7, colors=INK_MUTED)
    return tpc


def plot_scatter(ax, data):
    valid_fit = data['valid_fit'].ravel().astype(bool)
    softened = data['softened'].ravel().astype(bool)
    v_error_slope = data['v_error_slope'].ravel()
    strain_rate = data['strain_rate'].ravel()

    x_unchanged = v_error_slope[valid_fit & ~softened]
    y_unchanged = strain_rate[valid_fit & ~softened]
    x_softened = v_error_slope[valid_fit & softened]
    y_softened = strain_rate[valid_fit & softened]

    ax.scatter(x_unchanged, y_unchanged, s=18, color=CAT_UNCHANGED, alpha=0.6,
               edgecolors='none', label='Unchanged')
    ax.scatter(x_softened, y_softened, s=18, color=CAT_SOFTENED, alpha=0.7,
               edgecolors='none', label='Softened')

    fit_slope = float(data['fit_slope'].ravel()[0])
    fit_intercept = float(data['fit_intercept'].ravel()[0])
    fit_rsquare = float(data['fit_rsquare'].ravel()[0])
    fit_rmse = float(data['fit_rmse'].ravel()[0])
    fit_n = int(data['fit_n'].ravel()[0])

    if np.isfinite(fit_slope) and valid_fit.any():
        xx = np.linspace(v_error_slope[valid_fit].min(), v_error_slope[valid_fit].max(), 100)
        yy = fit_slope * xx + fit_intercept
        ax.plot(xx, yy, color=INK_PRIMARY, linewidth=1.5, linestyle='--',
                label='Best fit')
        sign = '+' if fit_intercept >= 0 else '-'
        eqn = 'strain rate = {:.4g} × error slope {} {:.4g}'.format(fit_slope, sign, abs(fit_intercept))
        stats = 'R² = {:.3f}   RMSE = {:.3g}   n = {}'.format(fit_rsquare, fit_rmse, fit_n)
        ax.text(0.02, 0.98, eqn + '\n' + stats, transform=ax.transAxes,
                fontsize=9, color=INK_SECONDARY, va='top', ha='left')

    ax.set_xlabel('Velocity error slope', fontsize=9, color=INK_SECONDARY)
    ax.set_ylabel('Strain rate', fontsize=9, color=INK_SECONDARY)
    ax.set_title('Strain rate vs. velocity error slope', fontsize=11, color=INK_PRIMARY)
    ax.tick_params(colors=INK_MUTED, labelsize=8)
    ax.grid(True, color=GRIDLINE, linewidth=0.7)
    for spine in ax.spines.values():
        spine.set_color(BASELINE)
    legend = ax.legend(loc='lower right', fontsize=8, frameon=False)
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)


def plot_info_panel(ax, data):
    ax.axis('off')
    catchment = int(data['catchment'].ravel()[0])
    was_softened = bool(data['was_softened'].ravel()[0])
    fit_rsquare = float(data['fit_rsquare'].ravel()[0])
    rsquare_threshold = float(data['rsquare_threshold'].ravel()[0])
    min_thickness_m = float(data['min_thickness_m'].ravel()[0])
    min_velocity_myr = float(data['min_velocity_myr'].ravel()[0])
    n_softened = int(data['softened'].ravel().astype(bool).sum())
    n_valid = int(data['valid_fit'].ravel().astype(bool).sum())

    lines = [
        'Catchment {}'.format(catchment),
        '',
        'Shear margins softened: {}'.format('yes' if was_softened else 'no'),
        'R² threshold: {:.2f}'.format(rsquare_threshold),
        '',
        'Softened vertices: {} / {} valid ({:.1f}%)'.format(
            n_softened, n_valid, 100.0 * n_softened / n_valid if n_valid else 0.0),
        '',
        'Fit inclusion mask:',
        '  ice-covered',
        '  thickness ≥ {:g} m'.format(min_thickness_m),
        '  velocity ≥ {:g} m/yr'.format(min_velocity_myr),
    ]

    if 'vel_cutoff' in data:
        vel_cutoff = float(data['vel_cutoff'].ravel()[0])
        vel_percentile = float(data['softening_vel_percentile'].ravel()[0])
        lines += [
            '',
            'Max cutoff velocity: {:.1f} m/yr'.format(vel_cutoff),
            '  ({:g}th percentile of catchment vertex'.format(vel_percentile),
            '  velocities, NOT area-weighted)',
        ]

    ax.text(0.0, 1.0, '\n'.join(lines), transform=ax.transAxes, fontsize=9,
            color=INK_SECONDARY, va='top', ha='left', family='monospace')


def build_summary_pdf(mat_path, out_path):
    data = loadmat(mat_path)
    catchment = int(data['catchment'].ravel()[0])
    triang, _ = make_triangulation(data)

    fig = plt.figure(figsize=(16, 11))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.3, 1, 1], hspace=0.45, wspace=0.4)

    was_softened = bool(data['was_softened'].ravel()[0])
    fig.suptitle('Shear margin softening summary — catchment {} ({})'.format(
        catchment, 'softened' if was_softened else 'not softened'),
        fontsize=15, color=INK_PRIMARY, y=0.98)

    ax_scatter = fig.add_subplot(gs[0, :3])
    ax_scatter.set_facecolor(SURFACE)
    plot_scatter(ax_scatter, data)

    ax_info = fig.add_subplot(gs[0, 3])
    ax_info.set_facecolor(SURFACE)
    plot_info_panel(ax_info, data)

    strain_rate = data['strain_rate'].ravel()
    v_error = data['v_error'].ravel()
    v_error_slope = data['v_error_slope'].ravel()
    rheology_before = data['rheology_before'].ravel()
    rheology_after = data['rheology_after'].ravel()
    rheology_diff = data['rheology_diff'].ravel()
    vel_diff = data['vel_diff'].ravel()

    rheo_vmin = np.nanmin([np.nanmin(rheology_before), np.nanmin(rheology_after)])
    rheo_vmax = np.nanmax([np.nanmax(rheology_before), np.nanmax(rheology_after)])

    xlim, ylim = compute_zoom_extent(data)

    ax = fig.add_subplot(gs[1, 0]); ax.set_facecolor(SURFACE)
    plot_diverging_map(fig, ax, triang, v_error, 'Velocity error (m/yr)\nmodel − observed', xlim, ylim,
                        absmax_scale=V_ERROR_COLOR_SCALE)

    ax = fig.add_subplot(gs[1, 1]); ax.set_facecolor(SURFACE)
    plot_sequential_map(fig, ax, triang, v_error_slope, 'Velocity error slope', xlim, ylim)

    ax = fig.add_subplot(gs[1, 2]); ax.set_facecolor(SURFACE)
    plot_sequential_map(fig, ax, triang, strain_rate, 'Strain rate (1/yr)', xlim, ylim)

    ax = fig.add_subplot(gs[1, 3]); ax.set_facecolor(SURFACE)
    ax.axis('off')

    ax = fig.add_subplot(gs[2, 0]); ax.set_facecolor(SURFACE)
    plot_sequential_map(fig, ax, triang, rheology_before, 'Rheology B, before', xlim, ylim,
                         vmin=rheo_vmin, vmax=rheo_vmax)

    ax = fig.add_subplot(gs[2, 1]); ax.set_facecolor(SURFACE)
    plot_sequential_map(fig, ax, triang, rheology_after, 'Rheology B, after', xlim, ylim,
                         vmin=rheo_vmin, vmax=rheo_vmax)

    ax = fig.add_subplot(gs[2, 2]); ax.set_facecolor(SURFACE)
    plot_diverging_map(fig, ax, triang, rheology_diff, 'Rheology B difference\nafter − before', xlim, ylim)

    ax = fig.add_subplot(gs[2, 3]); ax.set_facecolor(SURFACE)
    plot_diverging_map(fig, ax, triang, vel_diff, 'Stress balance velocity\ndifference, after − before (m/yr)', xlim, ylim)

    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print('Wrote {}'.format(out_path))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-dir', default='./shear_margin_softening_exports',
                        help='Directory of catchment_<id>.mat files from shear_margin_softening_analysis.m')
    parser.add_argument('--output-dir', default='./shear_margin_softening_summaries',
                        help='Directory to write catchment_<id>_summary.pdf files')
    parser.add_argument('--catchments', nargs='+', type=int,
                        help='Only process these catchment ids (default: all files found)')
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        raise SystemExit(
            "Input directory does not exist: {}\n"
            "Run shear_margin_softening_analysis.m first, or point --input-dir "
            "at wherever its catchment_<id>.mat exports were copied to.".format(args.input_dir))

    os.makedirs(args.output_dir, exist_ok=True)

    if args.catchments:
        mat_paths = [os.path.join(args.input_dir, 'catchment_{}.mat'.format(c)) for c in args.catchments]
    else:
        mat_paths = sorted(glob.glob(os.path.join(args.input_dir, 'catchment_*.mat')))

    if not mat_paths:
        raise SystemExit('Input directory exists but has no catchment_*.mat files: {}'.format(args.input_dir))

    for mat_path in mat_paths:
        if not os.path.exists(mat_path):
            print('Skipping missing file: {}'.format(mat_path))
            continue
        catchment_id = os.path.splitext(os.path.basename(mat_path))[0].replace('catchment_', '')
        out_path = os.path.join(args.output_dir, 'catchment_{}_summary.pdf'.format(catchment_id))
        build_summary_pdf(mat_path, out_path)


if __name__ == '__main__':
    main()
