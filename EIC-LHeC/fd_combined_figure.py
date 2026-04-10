# -*- coding: utf-8 -*-
"""
EIC-LHeC Final Doublet -- Combined Beam Envelope Figure
=========================================================
Produces a single two-panel figure:
  Top panel   : electron FD  (X and Y envelopes on the same axes)
  Bottom panel: proton  FD  (X and Y envelopes on the same axes)

No title header, no gradient footer -- those go into the paper caption/text.
Quadrupole strengths are taken directly from the fd_design.py solution.

Output: FD_combined.png

Usage:
    /opt/local/bin/python3.4 fd_combined_figure.py
"""

import numpy as np
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# FD geometry
# ─────────────────────────────────────────────────────────────────────────────
LSTAR = 2.0    # IP → entrance of Q1  [m]
LQ    = 2.0    # quadrupole length  [m]
LGAP  = 0.5    # gap between Q1 and Q2  [m]

# ─────────────────────────────────────────────────────────────────────────────
# Solved quadrupole strengths (from fd_design.py)
#   Electron: QD-QF  (k1 < 0, k2 > 0)
#   Proton:   QF-QD  (k1 > 0, k2 < 0)
# ─────────────────────────────────────────────────────────────────────────────
BEAMS = {
    'Electron': {
        'E_GeV':        10.0,
        'theta_x_urad': 30.6,
        'theta_y_urad': 149.9,
        'k1': -0.318996,   # m^-2
        'k2': +0.167134,
        'polarity': 'QD-QF',
    },
    'Proton': {
        'E_GeV':        275.0,
        'theta_x_urad': 306.3,
        'theta_y_urad': 130.6,
        'k1': +0.318996,   # m^-2
        'k2': -0.167134,
        'polarity': 'QF-QD',
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# Transfer matrices
# ─────────────────────────────────────────────────────────────────────────────

def drift(L):
    return np.array([[1., L], [0., 1.]])


def quad(k, L):
    """Thick-lens matrix; k > 0 focuses, k < 0 defocuses."""
    if abs(k) < 1e-14:
        return drift(L)
    if k > 0:
        sq = np.sqrt(k);   phi = sq * L
        return np.array([[ np.cos(phi),       np.sin(phi) / sq],
                         [-sq * np.sin(phi),  np.cos(phi)     ]])
    else:
        sq = np.sqrt(-k);  phi = sq * L
        return np.array([[ np.cosh(phi),      np.sinh(phi) / sq],
                         [ sq * np.sinh(phi), np.cosh(phi)     ]])


def trace(k1, k2, theta_x, theta_y, n=4000):
    """Trace 4 rays; returns s, xp, xm, yp, ym [m]."""
    L_tot = LSTAR + LQ + LGAP + LQ
    s     = np.linspace(0., L_tot, n)
    D0  = drift(LSTAR)
    Q1x = quad( k1, LQ);  Q1y = quad(-k1, LQ)
    Dg  = drift(LGAP)
    xp = np.empty(n); xm = np.empty(n)
    yp = np.empty(n); ym = np.empty(n)
    for i, si in enumerate(s):
        if si <= LSTAR:
            Mx = drift(si);  My = drift(si)
        elif si <= LSTAR + LQ:
            d  = si - LSTAR
            Mx = quad( k1, d).dot(D0)
            My = quad(-k1, d).dot(D0)
        elif si <= LSTAR + LQ + LGAP:
            d  = si - LSTAR - LQ
            Mx = drift(d).dot(Q1x.dot(D0))
            My = drift(d).dot(Q1y.dot(D0))
        else:
            d  = si - LSTAR - LQ - LGAP
            Mx = quad( k2, d).dot(Dg.dot(Q1x.dot(D0)))
            My = quad(-k2, d).dot(Dg.dot(Q1y.dot(D0)))
        xp[i] = Mx.dot([0.,  theta_x])[0]
        xm[i] = Mx.dot([0., -theta_x])[0]
        yp[i] = My.dot([0.,  theta_y])[0]
        ym[i] = My.dot([0., -theta_y])[0]
    return s, xp, xm, yp, ym


# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

def make_combined_figure(out_path):

    # ── font sizes ────────────────────────────────────────────────
    FS_AXIS   = 13    # axis labels
    FS_TICK   = 11    # tick labels
    FS_LEGEND = 11    # legend text
    FS_TITLE  = 12    # subplot title
    FS_QLABEL = 11    # Q1/Q2 labels inside plot
    FS_ANNOT  = 10.5  # max-size annotation boxes

    # Quadrupole landmark positions [m]
    q1s = LSTAR
    q1e = LSTAR + LQ
    q2s = LSTAR + LQ + LGAP
    q2e = LSTAR + LQ + LGAP + LQ

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    for ax, (beam_name, bp) in zip(axes, BEAMS.items()):

        k1 = bp['k1'];  k2 = bp['k2']
        tx = bp['theta_x_urad'] * 1e-6
        ty = bp['theta_y_urad'] * 1e-6

        s, xp, xm, yp, ym = trace(k1, k2, tx, ty)

        mx = max(np.max(np.abs(xp)), np.max(np.abs(xm)))
        my = max(np.max(np.abs(yp)), np.max(np.abs(ym)))

        # ── traces ──────────────────────────────────────────────
        lbl_xp = '+$\\theta_x$ = {0:.0f} $\\mu$rad'.format(
            bp['theta_x_urad'])
        lbl_yp = '+$\\theta_y$ = {0:.0f} $\\mu$rad'.format(
            bp['theta_y_urad'])

        ax.plot(s, xp * 1e3, color='steelblue', lw=1.8, label=lbl_xp)
        ax.plot(s, xm * 1e3, color='steelblue', lw=1.8, ls='--',
                label='$-\\theta_x$')
        ax.fill_between(s, xm * 1e3, xp * 1e3,
                        color='steelblue', alpha=0.13)

        ax.plot(s, yp * 1e3, color='seagreen', lw=1.8, label=lbl_yp)
        ax.plot(s, ym * 1e3, color='seagreen', lw=1.8, ls='--',
                label='$-\\theta_y$')
        ax.fill_between(s, ym * 1e3, yp * 1e3,
                        color='seagreen', alpha=0.13)

        # ── quadrupole shading ───────────────────────────────────
        q1_type = 'QF' if k1 > 0 else 'QD'
        q2_type = 'QF' if k2 > 0 else 'QD'
        c1 = '#AED6F1' if k1 > 0 else '#FAD7A0'
        c2 = '#AED6F1' if k2 > 0 else '#FAD7A0'

        ax.axvspan(q1s, q1e, color=c1, alpha=0.48, zorder=0,
                   label='Q1 ({0})'.format(q1_type))
        ax.axvspan(q2s, q2e, color=c2, alpha=0.48, zorder=0,
                   label='Q2 ({0})'.format(q2_type))
        for xs in (q1s, q1e, q2s, q2e):
            ax.axvline(xs, color='grey', lw=0.7, ls=':', zorder=1)
        ax.axvline(0., color='crimson', lw=1.2, ls='--', alpha=0.8,
                   label='IP', zorder=1)
        ax.axhline(0., color='black', lw=0.5, zorder=1)
        ax.grid(True, alpha=0.25, zorder=0)

        # ── axis labels & title ──────────────────────────────────
        ax.set_ylabel('x  or  y  [mm]', fontsize=FS_AXIS)
        ax.set_title(
            '{0} beam,  E = {1:.0f} GeV,  polarity {2}  '
            '|  max|x| = {3:.2f} mm,  max|y| = {4:.2f} mm'.format(
                beam_name, bp['E_GeV'], bp['polarity'],
                mx * 1e3, my * 1e3),
            fontsize=FS_TITLE)

        ax.tick_params(labelsize=FS_TICK)

        # ── max-size annotation boxes ────────────────────────────
        ax.text(0.97, 0.96,
                'max|x| = {0:.2f} mm'.format(mx * 1e3),
                transform=ax.transAxes,
                ha='right', va='top', fontsize=FS_ANNOT,
                color='steelblue', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.25',
                          facecolor='white', edgecolor='steelblue',
                          alpha=0.85))
        ax.text(0.97, 0.80,
                'max|y| = {0:.2f} mm'.format(my * 1e3),
                transform=ax.transAxes,
                ha='right', va='top', fontsize=FS_ANNOT,
                color='seagreen', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.25',
                          facecolor='white', edgecolor='seagreen',
                          alpha=0.85))

        # ── Q1 / Q2 labels ───────────────────────────────────────
        ylim = ax.get_ylim()
        yspan = ylim[1] - ylim[0]
        ytop = ylim[1] - 0.04 * yspan
        ax.text((q1s + q1e) / 2., ytop,
                'Q1\n({0})'.format(q1_type),
                ha='center', va='top', fontsize=FS_QLABEL,
                color='navy', fontweight='bold')
        ax.text((q2s + q2e) / 2., ytop,
                'Q2\n({0})'.format(q2_type),
                ha='center', va='top', fontsize=FS_QLABEL,
                color='navy', fontweight='bold')

        # ── legend ──────────────────────────────────────────────
        ax.legend(fontsize=FS_LEGEND, loc='upper left',
                  ncol=3, framealpha=0.85)

    axes[-1].set_xlabel('s  [m]   (s = 0 at IP)', fontsize=FS_AXIS)
    axes[-1].tick_params(labelsize=FS_TICK)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print('Saved --> {0}'.format(out_path))
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import os
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, 'FD_combined.png')
    make_combined_figure(out_path)
