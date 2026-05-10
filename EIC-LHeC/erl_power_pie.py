#!/usr/bin/env python3
"""
erl_power_pie.py  –  Single-pass ERL wall-plug power budget pie chart
Standalone figure for the IPAC WEV1301 paper.
Sized for a single JACoW column (~8.6 cm wide).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# ── Computed power values from erl_simulation.py ─────────────────────────
# Single-pass ERL, 10 GeV, 60 mA, 801.58 MHz RF
components = [
    ("Beam dump  (7 MeV \u00d7 60 mA)",         420.8),
    ("Injector RF  (\u03b7 = 35\u202f%)",        1202.4),
    ("SRF cryoplant  (714 cavities, 2\u202fK)",  4388.9),
    ("SR in return arc",                          88.9),
    ("Residual RF  (beam-loading)",              925.6),
]
labels = [c[0] for c in components]
values = np.array([c[1] for c in components])   # kW
total  = values.sum()                             # kW = 7026.6 kW ≈ 7.0 MW

# ── Colours (colour-blind-safe, matches simulation palette) ───────────────
colors = ["#1f77b4",   # blue   – dump
          "#17becf",   # cyan   – injector
          "#2ca02c",   # green  – cryo
          "#9467bd",   # purple – SR
          "#7f7f7f"]   # grey   – residual RF

# ── Figure layout ─────────────────────────────────────────────────────────
# JACoW single column: 8.6 cm = 3.386 in wide.
fig, ax = plt.subplots(figsize=(3.38, 2.60))   # inches: 8.59 × 6.6 cm
fig.subplots_adjust(left=0.0, right=0.56, top=0.97, bottom=0.02)

# ── Pie chart ─────────────────────────────────────────────────────────────
# Sort by value so legend reads large-to-small
order  = np.argsort(values)[::-1]
vals_s = values[order]
lbls_s = [labels[i] for i in order]
cols_s = [colors[i] for i in order]

# Explode the dominant slice (SRF cryoplant) slightly
explode = [0.06 if v == values.max() else 0.0 for v in vals_s]

wedges, texts, autotexts = ax.pie(
    vals_s,
    labels=None,
    colors=cols_s,
    explode=explode,
    autopct=lambda p: f"{p:.1f}%" if p > 4 else "",
    pctdistance=0.68,
    startangle=105,
    wedgeprops=dict(linewidth=0.8, edgecolor="white"),
    textprops=dict(fontsize=7.5),
)

for at in autotexts:
    at.set_fontsize(7.2)
    at.set_fontweight("bold")
    at.set_color("white")
    at.set_path_effects([pe.withStroke(linewidth=1.2, foreground="black")])

# ── Legend (right of pie) ─────────────────────────────────────────────────
legend_labels = [f"{lbl}\\n{v:.0f} kW" for lbl, v in
                 zip(lbls_s, vals_s)]
# Build cleaner single-line legend entries
legend_entries = [
    (w, f"{lbl}  —  {v:.0f} kW")
    for w, lbl, v in zip(wedges,
                         [l.replace("\n", " ") for l in lbls_s],
                         vals_s)
]
leg = ax.legend(
    [e[0] for e in legend_entries],
    [e[1] for e in legend_entries],
    loc="center left",
    bbox_to_anchor=(1.04, 0.50),
    fontsize=6.5,
    frameon=True,
    framealpha=0.95,
    edgecolor="#bbbbbb",
    handlelength=1.1,
    handleheight=0.95,
    labelspacing=0.45,
    borderpad=0.5,
)

# ── Centre annotation ─────────────────────────────────────────────────────
ax.text(0, 0.08, "Total", ha="center", va="center",
        fontsize=7.5, color="#333333")
ax.text(0, -0.18, f"{total/1e3:.2f} MW", ha="center", va="center",
        fontsize=9.5, fontweight="bold", color="#111111")

# ── Title ─────────────────────────────────────────────────────────────────
# Title removed — caption provided in LaTeX \caption{}

fig.savefig("ERL_power_pie.png",
            dpi=300, bbox_inches="tight",
            facecolor="white", edgecolor="none")
fig.savefig("ERL_power_pie.pdf",
            bbox_inches="tight",
            facecolor="white", edgecolor="none")

print(f"Saved ERL_power_pie.png  and  ERL_power_pie.pdf")
print(f"\nPower breakdown (single-pass ERL):")
for lbl, v in zip([l.replace('\n',' ') for l in lbls_s], vals_s):
    print(f"  {lbl:<38s}  {v:6.1f} kW  ({100*v/total:.1f}%)")
print(f"  {'TOTAL':<38s}  {total:6.1f} kW  =  {total/1e3:.3f} MW")
