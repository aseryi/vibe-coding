# -*- coding: utf-8 -*-
"""
EIC-LHeC Final Doublet -- Magnet Technology Analysis
======================================================
For each quadrupole in the electron and proton Final Doublets, computes:
  R_max = B_bore_max / G      (maximum aperture radius allowed by technology)
  N_X   = R_max / max_x_FD   (aperture / max horizontal beam size in FD)
  N_Y   = R_max / max_y_FD   (aperture / max vertical   beam size in FD)

Beam parameters and FD gradients come from fd_design.py results.
N >> 1 means comfortable margin; N ~ 1 means barely enough room.

Usage:
    /opt/local/bin/python3.4 fd_magnet_analysis.py
"""

# ── FD results from fd_design.py ─────────────────────────────────────────────
# Gradients and max beam sizes (1-sigma envelope) in each FD

FD = {
    'Proton (275 GeV)': {
        'Q1': {'G_Tm': 292.62, 'type': 'QF'},
        'Q2': {'G_Tm': 153.31, 'type': 'QD'},
        'max_x_mm': 0.818,
        'max_y_mm': 1.434,
    },
    'Electron (10 GeV)': {
        'Q1': {'G_Tm': 10.64, 'type': 'QD'},
        'Q2': {'G_Tm':  5.58, 'type': 'QF'},
        'max_x_mm': 0.336,
        'max_y_mm': 0.400,
    },
}

# ── Superconducting magnet technologies ──────────────────────────────────────
# B_bore_max: maximum achievable bore field (aperture field) [Tesla]
# For a quadrupole: B_bore = G * R  =>  R_max = B_bore_max / G

TECHNOLOGIES = [
    ('NbTi,  1.9 K', 8.5),   # superfluid helium
    ('NbTi,  4.2 K', 7.0),   # liquid helium
    ('Nb3Sn, 4.2 K', 12.0),  # liquid helium, advanced conductor
]

# ── Compute and print ─────────────────────────────────────────────────────────

def print_table(fd_name, fd_data):
    Q1 = fd_data['Q1'];  Q2 = fd_data['Q2']
    mx = fd_data['max_x_mm'];  my = fd_data['max_y_mm']

    print()
    print('=' * 88)
    print('  {0}'.format(fd_name))
    print('  G(Q1/{0}) = {1:.2f} T/m,   G(Q2/{2}) = {3:.2f} T/m'.format(
        Q1['type'], Q1['G_Tm'], Q2['type'], Q2['G_Tm']))
    print('  max|x| in FD = {0:.3f} mm,   max|y| in FD = {1:.3f} mm'.format(mx, my))
    print('=' * 88)

    hdr = '{:<18s}  {:>9s} {:>7s} {:>7s}   {:>9s} {:>7s} {:>7s}'.format(
        'Technology',
        'Rmax_Q1', 'NX_Q1', 'NY_Q1',
        'Rmax_Q2', 'NX_Q2', 'NY_Q2')
    print(hdr)
    print('{:<18s}  {:>9s} {:>7s} {:>7s}   {:>9s} {:>7s} {:>7s}'.format(
        '[B_bore_max]', '[mm]', '', '', '[mm]', '', ''))
    print('-' * 88)

    for tech, B in TECHNOLOGIES:
        R1 = B / Q1['G_Tm'] * 1e3   # mm
        R2 = B / Q2['G_Tm'] * 1e3
        NX1 = R1 / mx;  NY1 = R1 / my
        NX2 = R2 / mx;  NY2 = R2 / my
        print('{:<18s}  {:>9.1f} {:>7.1f} {:>7.1f}   {:>9.1f} {:>7.1f} {:>7.1f}'.format(
            tech, R1, NX1, NY1, R2, NX2, NY2))

    print()
    print('  Rmax = B_bore_max / G  [mm]')
    print('  NX   = Rmax / max|x|_FD,   NY = Rmax / max|y|_FD')


for fd_name, fd_data in FD.items():
    print_table(fd_name, fd_data)

print()
print('=' * 88)
print('  Note: NX, NY >> 1  =>  large margin (comfortable design)')
print('        NX, NY ~ 1   =>  aperture barely fits the beam core')
print('        For collider IR quads, typically NX, NY > 10 is required')
print('        (to accommodate orbit offsets, halo, and tolerances).')
print('=' * 88)
