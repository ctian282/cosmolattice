#!/usr/bin/env python3
"""Plot the simultaneous-nucleation scalar spectra in Fig. 6 coordinates.

The input directory must contain ``spectra_scalar_0.txt`` and
``average_spectra_times.txt`` from the ``vacuum_bubbles`` model.  The spectrum
file is a sequence of blank-line-separated measurements; its first column is
the wavenumber and its second column is the scalar power spectrum.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TARGET_TIMES = np.array([0.0, 0.47, 1.0, 4.0, 8.0])
CURVE_LABELS = ["$t/R_* = 0.0$", "$t/R_* = 0.47$", "$t/R_* = 1.0$", "$t/R_* = 4.0$", "$t/R_* = 8.0$"]
CURVE_COLORS = ["#1474b8", "#00a383", "#e57917", "#d37aa8", "#f0d329"]


def read_times(path: Path) -> np.ndarray:
    values = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.append(float(line.split()[0]))
    return np.asarray(values, dtype=float)


def read_spectrum_blocks(path: Path) -> list[np.ndarray]:
    blocks: list[np.ndarray] = []
    rows: list[tuple[float, float]] = []

    def finish() -> None:
        nonlocal rows
        if rows:
            blocks.append(np.asarray(rows, dtype=float))
            rows = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            finish()
            continue
        if line.startswith("#") or line[0].isalpha() or line.startswith("1:"):
            continue
        fields = line.split()
        if len(fields) >= 2:
            rows.append((float(fields[0]), float(fields[1])))
    finish()
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="vacuum_bubbles output directory")
    parser.add_argument("--output", type=Path, default=Path("fig6_left.png"), help="PNG path")
    parser.add_argument("--phi-b", type=float, default=2.45, help="bubble amplitude")
    parser.add_argument("--r-star", type=float, default=28.16, help="mean bubble separation R*")
    parser.add_argument("--x-max", type=float, default=3e3, help="upper x-axis limit (paper: 3e3)")
    args = parser.parse_args()

    times = read_times(args.output_dir / "average_spectra_times.txt")
    blocks = read_spectrum_blocks(args.output_dir / "spectra_scalar_0.txt")
    if len(times) != len(blocks):
        raise SystemExit(
            f"found {len(times)} spectrum times but {len(blocks)} spectrum blocks in {args.output_dir}"
        )
    if not blocks:
        raise SystemExit("the spectrum file contains no measurements")

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for target, label, color in zip(TARGET_TIMES, CURVE_LABELS, CURVE_COLORS):
        index = int(np.argmin(np.abs(times / args.r_star - target)))
        data = blocks[index]
        x = data[:, 0] * args.r_star
        y = data[:, 1] / args.phi_b**2
        positive = (x > 0) & (y > 0)
        ax.loglog(x[positive], y[positive], color=color, lw=1.45, label=label)

    ax.axvline(2 * np.pi, color="black", ls=":", lw=0.9)
    ax.set_xlim(1e-1, args.x_max)
    ax.set_ylim(1e-8, 1e0)
    ax.set_xlabel(r"$kR_*$")
    ax.set_ylabel(r"$P_{\phi}/\phi_b^2$")
    ax.legend(loc="lower left", frameon=True, framealpha=0.9, fontsize=8.5, handlelength=2.5)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.grid(False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.5, 0.015, "(a) Simultaneous", ha="center", va="bottom")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    print(f"wrote {args.output}")
    for target in TARGET_TIMES:
        index = int(np.argmin(np.abs(times / args.r_star - target)))
        print(f"t/R* target {target:g}: using t={times[index]:g} (t/R*={times[index] / args.r_star:g})")


if __name__ == "__main__":
    main()
