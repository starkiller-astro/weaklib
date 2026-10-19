#!/usr/bin/env python3
"""Plot the output of wlCalculateSoundSpeed (SoundSpeed.bin).

Compares Gamma1 and the sound speed from the Separate and Combined
interpolations of BaryonsPlusPhotonsPlusLeptonsEOS.h5 with the reference
values of wl-EOS-SFHo-15-25-50.h5, on the reference nodes.

Figures:
  gamma1_vs_rho_Ye<Ye>.png   Gamma1 vs density at a few temperatures
  gamma1_error_map_Ye<Ye>.png  relative Gamma1 error on the (rho, T) plane
  error_histograms.png       distribution of the relative errors, all nodes

Usage:
  python3 plot_sound_speed.py [SoundSpeed.bin] [--ye 0.3] [--out DIR]
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

# --- Colors: method identity is fixed across all figures ---
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
COLORS = {"Reference": INK, "Separate": "#2a78d6", "Combined": "#eb6834"}
STYLES = {"Reference": "-", "Separate": "--", "Combined": ":"}
METHODS = ["Reference", "Separate", "Combined"]

# Single-hue sequential ramp (light -> dark) for error magnitude
ERROR_CMAP = LinearSegmentedColormap.from_list(
    "error_blue",
    ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
ERROR_CMAP.set_bad("#f3f2ef")

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": INK_SECONDARY,
    "axes.labelcolor": INK,
    "axes.titlesize": 10,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
})


def read_bin(path):
    """Read SoundSpeed.bin (stream access, Fortran order)."""
    with open(path, "rb") as f:
        nD, nT, nYe = np.fromfile(f, dtype=np.int32, count=3)
        D = np.fromfile(f, dtype=np.float64, count=nD)
        T = np.fromfile(f, dtype=np.float64, count=nT)
        Ye = np.fromfile(f, dtype=np.float64, count=nYe)
        n = nD * nT * nYe * 3
        gamma = np.fromfile(f, dtype=np.float64, count=n).reshape((nD, nT, nYe, 3), order="F")
        cs = np.fromfile(f, dtype=np.float64, count=n).reshape((nD, nT, nYe, 3), order="F")
        negp = np.fromfile(f, dtype=np.int32, count=nD * nT * nYe).reshape((nD, nT, nYe), order="F")
    return D, T, Ye, gamma, cs, negp.astype(bool)


def rel_err(x, ref):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(x - ref) / ref


def plot_gamma_vs_rho(D, T, Ye, gamma, negp, iYe, out):
    """Small multiples: one panel per temperature, three lines per panel."""
    logT_targets = [9.5, 10.5, 11.0, 11.5]
    iTs = [int(np.argmin(np.abs(np.log10(T) - t))) for t in logT_targets]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True)
    logD = np.log10(D)
    for ax, iT in zip(axes.flat, iTs):
        # Shade the density range where the baryon P is <= 0 in the cell
        neg = negp[:, iT, iYe]
        if neg.any():
            ax.fill_between(logD, 0, 1, where=neg, transform=ax.get_xaxis_transform(),
                            color="#efeeea", linewidth=0, step="mid", zorder=0)
        for m, name in enumerate(METHODS):
            ax.plot(logD, gamma[:, iT, iYe, m], STYLES[name], color=COLORS[name],
                    label=name, zorder=3 - m * 0.5)
        ax.set_title(f"log$_{{10}}$ T = {np.log10(T[iT]):.2f} K")
        ax.set_ylabel(r"$\Gamma_1$")
        ax.set_yscale("log")
    for ax in axes[-1]:
        ax.set_xlabel(r"log$_{10}\,\rho$ [g cm$^{-3}$]")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    handles.append(plt.Rectangle((0, 0), 1, 1, color="#efeeea"))
    labels.append(r"cells with $P_\mathrm{bary} \leq 0$")
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4)
    fig.suptitle(rf"$\Gamma_1$ vs density, Y$_e$ = {Ye[iYe]:.3f}", y=1.07)

    path = os.path.join(out, f"gamma1_vs_rho_Ye{Ye[iYe]:.3f}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_error_maps(D, T, Ye, gamma, negp, iYe, out):
    """Relative Gamma1 error on the (log rho, log T) plane, shared color scale."""
    logD, logT = np.log10(D), np.log10(T)
    errs = {name: rel_err(gamma[:, :, iYe, m], gamma[:, :, iYe, 0])
            for m, name in enumerate(METHODS) if m > 0}
    norm = LogNorm(vmin=1e-4, vmax=1.0, clip=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True, constrained_layout=True)
    for ax, (name, err) in zip(axes, errs.items()):
        ax.grid(False)
        mesh = ax.pcolormesh(logD, logT, np.ma.masked_invalid(err).T, cmap=ERROR_CMAP,
                             norm=norm, shading="nearest", rasterized=True)
        ax.contour(logD, logT, negp[:, :, iYe].T.astype(float), levels=[0.5],
                   colors=INK, linewidths=1.0)
        ax.set_title(f"{name}:  median {np.nanmedian(err):.1e},  "
                     f"> 10 % at {100 * np.nanmean(err > 0.1):.1f} % of nodes")
        ax.set_xlabel(r"log$_{10}\,\rho$ [g cm$^{-3}$]")
    axes[0].set_ylabel(r"log$_{10}$ T [K]")
    cbar = fig.colorbar(mesh, ax=axes, extend="both", pad=0.01)
    cbar.set_label(r"$|\Gamma_1 - \Gamma_1^\mathrm{ref}| / \Gamma_1^\mathrm{ref}$")
    cbar.outline.set_visible(False)
    fig.suptitle(rf"Relative $\Gamma_1$ error, Y$_e$ = {Ye[iYe]:.3f}  "
                 r"(black contour: cells with $P_\mathrm{bary} \leq 0$; gray: no data)")

    path = os.path.join(out, f"gamma1_error_map_Ye{Ye[iYe]:.3f}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_histograms(gamma, cs, out):
    """Distribution of log10 relative error over all nodes, both methods."""
    bins = np.linspace(-8, 2, 101)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True, constrained_layout=True)
    for ax, (label, X) in zip(axes, [(r"$\Gamma_1$", gamma), (r"$c_s$", cs)]):
        for m, name in enumerate(METHODS):
            if m == 0:
                continue
            err = rel_err(X[..., m], X[..., 0]).ravel()
            err = err[np.isfinite(err)]
            le = np.log10(np.maximum(err, 1e-8))
            ax.hist(le, bins=bins, histtype="step", linewidth=2, color=COLORS[name],
                    label=f"{name} (median {np.median(err):.1e}, mean {np.mean(err):.1e})")
        ax.set_xlabel(f"log$_{{10}}$ relative error in {label}")
        ax.set_title(f"{label}, all reference nodes")
        ax.legend(loc="upper left")
    axes[0].set_ylabel("number of nodes")

    path = os.path.join(out, "error_histograms.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", nargs="?", default=os.path.join(here, "Executables", "SoundSpeed.bin"))
    parser.add_argument("--ye", type=float, default=0.3, help="Ye of the slices (nearest node)")
    parser.add_argument("--out", default=None, help="output directory (default: next to the .bin)")
    args = parser.parse_args()

    out = args.out or os.path.dirname(os.path.abspath(args.path))
    os.makedirs(out, exist_ok=True)

    D, T, Ye, gamma, cs, negp = read_bin(args.path)
    iYe = int(np.argmin(np.abs(Ye - args.ye)))

    for path in (plot_gamma_vs_rho(D, T, Ye, gamma, negp, iYe, out),
                 plot_error_maps(D, T, Ye, gamma, negp, iYe, out),
                 plot_histograms(gamma, cs, out)):
        print("wrote", path)


if __name__ == "__main__":
    main()
