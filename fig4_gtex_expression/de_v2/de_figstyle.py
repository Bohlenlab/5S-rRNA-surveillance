# -----------------------------------------------------------------------------
# de_figstyle.py — Shared matplotlib style for the differential-expression figures.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Shared figure style: vector PDF, 4x4-inch panels, 300 DPI, 8 pt font,
1 pt lines, editable text (fonttype 42)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANEL = 4.0  # inches per panel


def setup():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,          # editable vector text, not outlines
        "svg.fonttype": "none",
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "figure.titlesize": 8,
        "axes.linewidth": 1.0, "lines.linewidth": 1.0,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "patch.linewidth": 1.0, "grid.linewidth": 1.0,
        "lines.markersize": 4.0,
        "savefig.dpi": 300, "figure.dpi": 300, "savefig.bbox": "tight",
        "axes.spines.top": False, "axes.spines.right": False,
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    })


def save(fig, path):
    fig.savefig(path, format="pdf")  # vector PDF
    plt.close(fig)
    print("wrote", path)
