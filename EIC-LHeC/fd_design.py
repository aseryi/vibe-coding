# -*- coding: utf-8 -*-
"""
EIC-LHeC Final Doublet (FD) Design
====================================
Traces beam from the IP (diverging point source) through a two-quadrupole
Final Doublet and adjusts quadrupole strengths so the beam exits as a
parallel beam in both planes simultaneously.

Beam parameters from Table 2 of WEV1301
(Hybrid EIC-LHeC, Initial Configuration, ~2e33 cm^-2 s^-1).

Physics
-------
* Transfer matrices: thick-lens (exact for uniform quadrupoles).
* "Parallel beam" condition: M[1,1] = 0 (0-indexed) in both x and y planes.
  A ray starting at position 0 with angle theta exits with angle 0.
* Four rays: initial angles +/-theta_x (horizontal) and +/-theta_y (vertical),
  all starting at position 0 at the IP.
* Both QF-QD and QD-QF polarities are tried; the one giving the smallest
  maximum beam size inside the FD is chosen.

Usage
-----
    /opt/local/bin/python3.4 fd_design.py
"""

import numpy as np
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# ─────────────────────────────────────────────────────────────────────────────
# Beam parameters  (Table 2, Hybrid EIC-LHeC, Initial Config)
# ─────────────────────────────────────────────────────────────────────────────
BEAMS = {
    'Electron': {
        'E_GeV':        10.0,
        'mass_GeV':     0.000511,
        'theta_x_urad': 30.6,    # horizontal divergence at IP [urad]
        'theta_y_urad': 149.9,   # vertical   divergence at IP [urad]
    },
    'Proton': {
        'E_GeV':        275.0,
        'mass_GeV':     0.938272,
        'theta_x_urad': 306.3,   # horizontal divergence at IP [urad]
        'theta_y_urad': 130.6,   # vertical   divergence at IP [urad]
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# FD geometry  (initial assumptions)
# ─────────────────────────────────────────────────────────────────────────────
LSTAR = 2.0   # IP to entrance of Q1  [m]
LQ    = 2.0   # length of each quadrupole  [m]
LGAP  = 0.5   # drift gap between Q1 and Q2  [m]


# ─────────────────────────────────────────────────────────────────────────────
# Thick-lens transfer matrices
# ─────────────────────────────────────────────────────────────────────────────

def drift(L):
    """Transfer matrix for a drift of length L."""
    return np.array([[1., L], [0., 1.]])


def quad(k, L):
    """
    Thick-lens transfer matrix for a quadrupole of length L.
    k > 0 -> focusing   (cos/sin functions)
    k < 0 -> defocusing (cosh/sinh functions)
    k = 0 -> pure drift
    Convention: k > 0 focuses in the plane of this matrix.
    """
    if abs(k) < 1e-14:
        return drift(L)
    if k > 0:
        sq  = np.sqrt(k)
        phi = sq * L
        return np.array([[ np.cos(phi),       np.sin(phi) / sq],
                         [-sq * np.sin(phi),  np.cos(phi)     ]])
    else:
        sq  = np.sqrt(-k)
        phi = sq * L
        return np.array([[ np.cosh(phi),      np.sinh(phi) / sq],
                         [ sq * np.sinh(phi), np.cosh(phi)     ]])


def sys_matrix(k1, k2):
    """
    Full transfer matrices (x-plane, y-plane) from IP to exit of Q2.
    Convention: k > 0 focuses in x (horizontal); Q1 is closer to IP.
    """
    D0  = drift(LSTAR)
    Q1x = quad( k1, LQ)
    Q1y = quad(-k1, LQ)   # opposite sign in y-plane
    Dg  = drift(LGAP)
    Q2x = quad( k2, LQ)
    Q2y = quad(-k2, LQ)
    Mx  = Q2x.dot(Dg.dot(Q1x.dot(D0)))
    My  = Q2y.dot(Dg.dot(Q1y.dot(D0)))
    return Mx, My


# ─────────────────────────────────────────────────────────────────────────────
# Parallel-beam condition:  M[1,1] = 0 in both planes
# ─────────────────────────────────────────────────────────────────────────────

def residuals(k_abs, sign1, sign2):
    """
    Residuals passed to fsolve.
    sign1, sign2 (+1 or -1) fix the polarity of Q1 and Q2.
    k_abs contains the magnitudes.
    """
    k1 = sign1 * abs(k_abs[0])
    k2 = sign2 * abs(k_abs[1])
    Mx, My = sys_matrix(k1, k2)
    return [Mx[1, 1], My[1, 1]]


# ─────────────────────────────────────────────────────────────────────────────
# Ray tracing: 4 trajectories through the full beamline
# ─────────────────────────────────────────────────────────────────────────────

def trace(k1, k2, theta_x, theta_y, n=4000):
    """
    Trace 4 rays from the IP through the FD.
    Initial conditions:
        x = 0,  x' = +/-theta_x  (horizontal)
        y = 0,  y' = +/-theta_y  (vertical)
    Returns
        s          position along beamline [m]
        xp, xm     horizontal trajectories for +/- theta_x [m]
        yp, ym     vertical   trajectories for +/- theta_y [m]
    """
    L_tot = LSTAR + LQ + LGAP + LQ
    s     = np.linspace(0., L_tot, n)

    # Precompute segment-end matrices
    D0  = drift(LSTAR)
    Q1x = quad( k1, LQ)
    Q1y = quad(-k1, LQ)
    Dg  = drift(LGAP)

    xp = np.empty(n)
    xm = np.empty(n)
    yp = np.empty(n)
    ym = np.empty(n)

    for i, si in enumerate(s):
        if si <= LSTAR:                           # drift L*
            Mx = drift(si)
            My = drift(si)
        elif si <= LSTAR + LQ:                    # inside Q1
            d  = si - LSTAR
            Mx = quad( k1, d).dot(D0)
            My = quad(-k1, d).dot(D0)
        elif si <= LSTAR + LQ + LGAP:             # gap drift
            d  = si - LSTAR - LQ
            Mx = drift(d).dot(Q1x.dot(D0))
            My = drift(d).dot(Q1y.dot(D0))
        else:                                     # inside Q2
            d  = si - LSTAR - LQ - LGAP
            Mx = quad( k2, d).dot(Dg.dot(Q1x.dot(D0)))
            My = quad(-k2, d).dot(Dg.dot(Q1y.dot(D0)))

        xp[i] = Mx.dot([0.,  theta_x])[0]
        xm[i] = Mx.dot([0., -theta_x])[0]
        yp[i] = My.dot([0.,  theta_y])[0]
        ym[i] = My.dot([0., -theta_y])[0]

    return s, xp, xm, yp, ym


# ─────────────────────────────────────────────────────────────────────────────
# Solver: try both polarities, pick smallest max beam size
# ─────────────────────────────────────────────────────────────────────────────

def solve_fd(theta_x, theta_y):
    """
    Find k1, k2 giving a parallel beam exit in both planes.
    Tries QF-QD and QD-QF polarities over a grid of initial guesses.
    Returns a dict with solution data, or None if nothing converges.
    """
    best = {'max_size': np.inf}

    init_vals = [0.01, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]

    for sign1, sign2, label in [(+1, -1, 'QF-QD'), (-1, +1, 'QD-QF')]:
        for ki in init_vals:
            for kj in init_vals:
                try:
                    sol, _, ier, _ = fsolve(
                        residuals, [ki, kj],
                        args=(sign1, sign2),
                        full_output=True,
                        xtol=1e-13,
                    )
                    if ier != 1:
                        continue
                    if any(x <= 1e-6 for x in sol):
                        continue
                    res = residuals(sol, sign1, sign2)
                    if max(abs(r) for r in res) > 1e-8:
                        continue

                    k1 = sign1 * sol[0]
                    k2 = sign2 * sol[1]

                    s, xp, xm, yp, ym = trace(k1, k2, theta_x, theta_y)
                    mx = max(np.max(np.abs(xp)), np.max(np.abs(xm)))
                    my = max(np.max(np.abs(yp)), np.max(np.abs(ym)))
                    overall = max(mx, my)

                    if overall < best['max_size']:
                        best = {
                            'max_size': overall,
                            'k1': k1, 'k2': k2,
                            's': s,
                            'xp': xp, 'xm': xm,
                            'yp': yp, 'ym': ym,
                            'mx': mx, 'my': my,
                            'label': label,
                        }
                except Exception:
                    pass

    return best if best['max_size'] < np.inf else None


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def make_plot(beam_name, bpar, sol, out_path):
    """Create and save the X/Y envelope plot for one beam."""

    k1, k2   = sol['k1'], sol['k2']
    s        = sol['s']
    xp, xm   = sol['xp'], sol['xm']
    yp, ym   = sol['yp'], sol['ym']
    mx, my   = sol['mx'], sol['my']
    polarity = sol['label']

    # Magnetic rigidity  Brho [T*m]
    Brho = np.sqrt(bpar['E_GeV']**2 - bpar['mass_GeV']**2) / 0.29979
    G1 = abs(k1) * Brho
    G2 = abs(k2) * Brho

    # Element positions
    q1s = LSTAR
    q1e = LSTAR + LQ
    q2s = LSTAR + LQ + LGAP
    q2e = LSTAR + LQ + LGAP + LQ

    q1_type = 'QF' if k1 > 0 else 'QD'
    q2_type = 'QF' if k2 > 0 else 'QD'

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    fig.suptitle(
        '{0} Beam -- Final Doublet Envelope\n'
        'E = {1:.0f} GeV,  L* = {2:.1f} m,  Lq = {3:.1f} m,  '
        'gap = {4:.1f} m,  polarity: {5}'.format(
            beam_name, bpar['E_GeV'], LSTAR, LQ, LGAP, polarity),
        fontsize=12, y=0.99,
    )

    def draw_elements(ax):
        """Draw quadrupole shading, IP marker, axis lines, grid."""
        c1 = '#AED6F1' if k1 > 0 else '#FAD7A0'   # blue=QF, orange=QD
        c2 = '#AED6F1' if k2 > 0 else '#FAD7A0'
        ax.axvspan(q1s, q1e, color=c1, alpha=0.50, zorder=0,
                   label='Q1 ({0})'.format(q1_type))
        ax.axvspan(q2s, q2e, color=c2, alpha=0.50, zorder=0,
                   label='Q2 ({0})'.format(q2_type))
        for xs in (q1s, q1e, q2s, q2e):
            ax.axvline(xs, color='grey', lw=0.7, ls=':', zorder=1)
        ax.axvline(0., color='crimson', lw=1.2, ls='--', alpha=0.8,
                   label='IP', zorder=1)
        ax.axhline(0., color='black', lw=0.5, zorder=1)
        ax.grid(True, alpha=0.25, zorder=0)

    def add_quad_labels(ax):
        ylim = ax.get_ylim()
        span = ylim[1] - ylim[0]
        y_label = ylim[1] - 0.06 * span
        ax.text((q1s + q1e) / 2., y_label,
                'Q1\n({0})'.format(q1_type),
                ha='center', va='top', fontsize=8,
                color='navy', fontweight='bold')
        ax.text((q2s + q2e) / 2., y_label,
                'Q2\n({0})'.format(q2_type),
                ha='center', va='top', fontsize=8,
                color='navy', fontweight='bold')

    def max_size_box(ax, max_mm, color):
        """Add a text box showing the maximum beam size."""
        ax.text(0.97, 0.95,
                'max = {0:.2f} mm'.format(max_mm),
                transform=ax.transAxes,
                ha='right', va='top', fontsize=10,
                color=color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='white', edgecolor=color,
                          alpha=0.85))

    # ── Horizontal plane ─────────────────────────────────────────────────
    ax1.plot(s, xp * 1e3, color='steelblue', lw=1.8,
             label='+$\\theta_x$ = {0:.1f} $\\mu$rad'.format(
                 bpar['theta_x_urad']))
    ax1.plot(s, xm * 1e3, color='steelblue', lw=1.8, ls='--',
             label='$-\\theta_x$')
    ax1.fill_between(s, xm * 1e3, xp * 1e3,
                     color='steelblue', alpha=0.12)
    draw_elements(ax1)
    ax1.set_ylabel('x  [mm]', fontsize=11)
    ax1.set_title(
        'Horizontal plane  --  max |x| in FD = {0:.2f} mm'.format(mx * 1e3),
        fontsize=10)
    ax1.legend(fontsize=8, loc='upper left', ncol=3)
    max_size_box(ax1, mx * 1e3, 'steelblue')
    add_quad_labels(ax1)

    # ── Vertical plane ───────────────────────────────────────────────────
    ax2.plot(s, yp * 1e3, color='seagreen', lw=1.8,
             label='+$\\theta_y$ = {0:.1f} $\\mu$rad'.format(
                 bpar['theta_y_urad']))
    ax2.plot(s, ym * 1e3, color='seagreen', lw=1.8, ls='--',
             label='$-\\theta_y$')
    ax2.fill_between(s, ym * 1e3, yp * 1e3,
                     color='seagreen', alpha=0.12)
    draw_elements(ax2)
    ax2.set_ylabel('y  [mm]', fontsize=11)
    ax2.set_xlabel('s  [m]   (s = 0 at IP)', fontsize=11)
    ax2.set_title(
        'Vertical plane  --  max |y| in FD = {0:.2f} mm'.format(my * 1e3),
        fontsize=10)
    ax2.legend(fontsize=8, loc='upper left', ncol=3)
    max_size_box(ax2, my * 1e3, 'seagreen')
    add_quad_labels(ax2)

    # Footer: magnet parameters
    foot = (
        'k1 = {0:+.5f} m**-2  (|G1| = {1:.1f} T/m),   '
        'k2 = {2:+.5f} m**-2  (|G2| = {3:.1f} T/m),   '
        'Brho = {4:.1f} T*m    '
        '[overall max beam size = {5:.2f} mm]'.format(
            k1, G1, k2, G2, Brho, sol['max_size'] * 1e3)
    )
    fig.text(
        0.5, 0.003, foot,
        ha='center', va='bottom', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.3',
                  facecolor='#FFFDE7', alpha=0.85),
    )

    plt.tight_layout(rect=[0, 0.045, 1, 0.97])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print('  Saved --> {0}'.format(out_path))


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    import os
    out_dir = os.path.dirname(os.path.abspath(__file__))

    print('\nEIC-LHeC Final Doublet Design')
    print('Geometry:  L* = {0} m,  Lq = {1} m,  gap = {2} m\n'.format(
        LSTAR, LQ, LGAP))

    for beam_name, bpar in BEAMS.items():
        tx = bpar['theta_x_urad'] * 1e-6   # [rad]
        ty = bpar['theta_y_urad'] * 1e-6

        print('{0}'.format('-' * 60))
        print('{0} beam:  E = {1} GeV'.format(beam_name, bpar['E_GeV']))
        print('  theta_x = {0:.1f} urad,  theta_y = {1:.1f} urad'.format(
            bpar['theta_x_urad'], bpar['theta_y_urad']))

        sol = solve_fd(tx, ty)

        if sol is None:
            print('  *** No solution found! ***')
            continue

        Brho = np.sqrt(bpar['E_GeV']**2 - bpar['mass_GeV']**2) / 0.29979
        G1 = abs(sol['k1']) * Brho
        G2 = abs(sol['k2']) * Brho

        print('  Best polarity  : {0}'.format(sol['label']))
        print('  k1 = {0:+.6f} m**-2   |G1| = {1:.2f} T/m'.format(
            sol['k1'], G1))
        print('  k2 = {0:+.6f} m**-2   |G2| = {1:.2f} T/m'.format(
            sol['k2'], G2))
        print('  Brho = {0:.2f} T*m'.format(Brho))
        print('  Max |x| in FD  = {0:.3f} mm'.format(sol['mx'] * 1e3))
        print('  Max |y| in FD  = {0:.3f} mm'.format(sol['my'] * 1e3))
        print('  Overall max    = {0:.3f} mm'.format(sol['max_size'] * 1e3))

        fname = os.path.join(out_dir, 'FD_{0}_beam.png'.format(
            beam_name.lower()))
        make_plot(beam_name, bpar, sol, fname)

    print('\n{0}'.format('-' * 60))
    print('Done.')
    plt.show()    # show all figures at once after saving both


if __name__ == '__main__':
    main()
