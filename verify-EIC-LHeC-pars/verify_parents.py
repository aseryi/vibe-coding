"""
Independent verification of LHeC and EIC parameters quoted in Table 1 of
"ACCELERATOR DESIGN EDUCATIONAL PRIMER -- CONCEPTUALIZING AND OPTIMIZING
THE HYBRID LHEC-LIKE ELECTRON-ION COLLIDER DESIGN" (IPAC'26 draft, WEV1301).

References used by the paper:
  [6] EIC Parameter List Summary, S. Peggs, T. Satogata, April 4, 2026
      https://eic.jlab.org/Documents/EIC-General/EIC_ParameterList.pdf
  [7] Phase-One LHeC submission to ESPP 2026, Piotrzkowski et al.
      https://indico.cern.ch/event/1439855/contributions/6461558/   (cf. arXiv:2503.20475)

For each parent collider we recompute, from the parameters quoted in the
paper's Table 1, the derived quantities:
   - geometric (un-normalised) RMS emittances
   - RMS spot sizes sigma_x*, sigma_y* at the IP  (= sqrt(beta* * eps_geo))
   - electron and proton beam currents from N, f_collision
   - peak head-on luminosity for unequal Gaussian beams:
        L = f_coll * N_e * N_p / (2 pi sqrt(sx_e^2+sx_p^2) sqrt(sy_e^2+sy_p^2))
     (this reduces to f_coll N_e N_p / (4 pi sigma_x sigma_y) when the two
      beams are matched in size at the IP, which is the equal-beam case)
   - proton beam-beam tune shifts xi_x, xi_y  using
        xi_x = N_e r_p beta_p_x / [2 pi gamma_p sx_e (sx_e+sy_e)]
   - centre-of-mass energy sqrt(s) = 2 sqrt(E_e E_p)

The values found are then compared to the values quoted in Table 1.
Inconsistencies in Table 1 are flagged.

Note on luminosity units:  formula uses SI metres -> result in m^-2 s^-1.
                          Convert with 1 m^-2 = 1e-4 cm^-2  (since 1 m^2 = 1e4 cm^2).
"""

import math

# Physical constants (PDG)
M_E   = 0.5109989461e-3   # GeV/c^2
M_P   = 0.93827208816     # GeV/c^2
R_E   = 2.8179403262e-15  # m  (classical electron radius)
R_P   = 1.5346982671e-18  # m  (classical proton radius = R_e * m_e/m_p)
ELEM  = 1.602176634e-19   # C

CM2_PER_M2 = 1e-4         # multiply m^-2 luminosity by this to get cm^-2

def gamma_from_E(E_GeV, m_GeV):
    return E_GeV / m_GeV

def geo_emit(eps_norm_um, gamma):
    """Geometric RMS emittance [m] from normalised eps [µm·rad]."""
    return eps_norm_um * 1e-6 / gamma

def sigma_at_IP(eps_norm_um, beta_m, gamma):
    return math.sqrt(geo_emit(eps_norm_um, gamma) * beta_m)

def luminosity_unequal(N1, N2, f_coll, sx1, sy1, sx2, sy2):
    """
    Head-on luminosity for unequal Gaussian beams:
        L = f_coll N1 N2 / [ 2 pi sqrt(sx1^2+sx2^2) sqrt(sy1^2+sy2^2) ]
    Result has units of  N^2 m^-2 s^-1  -> convert at print time.
    """
    Sx = math.sqrt(sx1*sx1 + sx2*sx2)
    Sy = math.sqrt(sy1*sy1 + sy2*sy2)
    return f_coll * N1 * N2 / (2.0 * math.pi * Sx * Sy)

def beam_current_mA(N, f_coll):
    return ELEM * N * f_coll * 1e3

def xi_proton(N_e, beta_p_x, beta_p_y, gamma_p, sx_e, sy_e):
    pref = N_e * R_P / (2.0 * math.pi * gamma_p)
    return pref * beta_p_x / (sx_e * (sx_e + sy_e)),\
           pref * beta_p_y / (sy_e * (sx_e + sy_e))

def xi_electron(N_p, beta_e_x, beta_e_y, gamma_e, sx_p, sy_p):
    pref = N_p * R_E / (2.0 * math.pi * gamma_e)
    return pref * beta_e_x / (sx_p * (sx_p + sy_p)),\
           pref * beta_e_y / (sy_p * (sx_p + sy_p))

def sqrt_s(E_e, E_p):
    return 2.0 * math.sqrt(E_e * E_p)

def line(label, computed, quoted, fmt="{:>9.3f}", units=""):
    """Print computed vs quoted in aligned columns."""
    if isinstance(computed, tuple):
        c1, c2 = computed
        q1, q2 = quoted
        print(f"  {label:<32}  {fmt.format(c1)} / {fmt.format(c2)} {units:<10}"
              f"  Table 1: {q1} / {q2}")
    else:
        print(f"  {label:<32}  {fmt.format(computed):>20} {units:<10}"
              f"  Table 1: {quoted}")


# ==================================================================
# 1) LHeC Phase-1  (20 GeV e on 7 TeV p)  --  reference [7]
# ==================================================================
print("=" * 80)
print(" LHeC Phase-1 (20 GeV e × 7 TeV p) -- as quoted in Table 1 of IPAC26 paper")
print("=" * 80)

E_e   = 20.0;     E_p   = 7000.0
n_b   = 2808
f_coll = 4.0e7                  # Hz  (Table 1; original LHeC CDR uses 32 MHz)
N_e   = 1.19e10
N_p   = 2.2e11
eps_e_n_x = 22.0; eps_e_n_y = 22.0     # µm
eps_p_n_x =  2.5; eps_p_n_y =  2.5
beta_e_x  = 0.20; beta_e_y  = 0.20
beta_p_x  = 0.35; beta_p_y  = 0.35

gam_e = gamma_from_E(E_e, M_E);  gam_p = gamma_from_E(E_p, M_P)

sx_e = sigma_at_IP(eps_e_n_x, beta_e_x, gam_e)
sy_e = sigma_at_IP(eps_e_n_y, beta_e_y, gam_e)
sx_p = sigma_at_IP(eps_p_n_x, beta_p_x, gam_p)
sy_p = sigma_at_IP(eps_p_n_y, beta_p_y, gam_p)

I_e = beam_current_mA(N_e, f_coll)
I_p = beam_current_mA(N_p, f_coll)
sqs  = sqrt_s(E_e, E_p)

L_m  = luminosity_unequal(N_e, N_p, f_coll, sx_e, sy_e, sx_p, sy_p)
L_cm = L_m * CM2_PER_M2

xi_p = xi_proton(N_e, beta_p_x, beta_p_y, gam_p, sx_e, sy_e)
xi_e = xi_electron(N_p, beta_e_x, beta_e_y, gam_e, sx_p, sy_p)

print(f"  Inputs:  E_e={E_e} GeV, E_p={E_p} GeV, n_b={n_b}, f_coll={f_coll:.2e} Hz")
print(f"           N_e={N_e:.2e},  N_p={N_p:.2e}")
print(f"           eps_n^e={eps_e_n_x}/{eps_e_n_y} µm,  eps_n^p={eps_p_n_x}/{eps_p_n_y} µm")
print(f"           beta*_e={beta_e_x}/{beta_e_y} m,    beta*_p={beta_p_x}/{beta_p_y} m\n")

line("sigma*_e (h/v)",     (sx_e*1e6, sy_e*1e6), ("10.6", "10.6"), "{:>5.2f}", "µm")
line("sigma*_p (h/v)",     (sx_p*1e6, sy_p*1e6), ("10.7", "10.7"), "{:>5.2f}", "µm")
line("I_e",                I_e, "60.1",  "{:>7.2f}", "mA")
line("I_p",                I_p, "1110",  "{:>7.2f}", "mA")
line("sqrt(s)",            sqs, "749",   "{:>7.2f}", "GeV")
line("L (cm^-2 s^-1)",     L_cm,"5.78e33","{:>9.3e}", "")
line("xi_p (h/v)",         xi_p, ("0.001","0.001"),   "{:>7.4f}", "")
line("xi_e (h/v)",         xi_e, ("--","--"),         "{:>7.3f}", "  (not quoted)")

# Cross-check on f_coll consistency
print()
print("  Consistency check on f_coll vs quoted currents:")
f_from_Ie = 60.1e-3 / (ELEM * N_e)
f_from_Ip = 1110e-3 / (ELEM * N_p)
print(f"     f_coll implied by I_e=60.1 mA & N_e=1.19e10 :  {f_from_Ie:.3e} Hz "
      f" -> {f_from_Ie/1e6:.1f} MHz")
print(f"     f_coll implied by I_p=1110 mA & N_p=2.2e11  :  {f_from_Ip:.3e} Hz "
      f" -> {f_from_Ip/1e6:.1f} MHz")
print( "     ==> Table 1 quotes 40 MHz, but the I_e and I_p values are "
       "consistent with 32 MHz")
print( "     ==> note matches the comment in the .tex source: 'CHECK: LHeC rep "
       "rate 40 MHz and why initial Grok code gave it 32 MHz value?'\n")


# ==================================================================
# 2) Baseline EIC  (10 GeV e on 275 GeV p)  --  reference [6]
# ==================================================================
print("=" * 80)
print(" Baseline EIC (10 GeV e × 275 GeV p, 1e34 design point)")
print("=" * 80)

E_e   = 10.0;    E_p = 275.0
n_b   = 1160
f_coll = 9.1e7
N_e   = 1.72e11
N_p   = 6.9e10
eps_e_n_x = 391.0; eps_e_n_y = 26.0
eps_p_n_x =   3.3; eps_p_n_y =  0.3
beta_e_x  = 0.45; beta_e_y  = 0.06
beta_p_x  = 0.80; beta_p_y  = 0.07

gam_e = gamma_from_E(E_e, M_E);   gam_p = gamma_from_E(E_p, M_P)

sx_e = sigma_at_IP(eps_e_n_x, beta_e_x, gam_e)
sy_e = sigma_at_IP(eps_e_n_y, beta_e_y, gam_e)
sx_p = sigma_at_IP(eps_p_n_x, beta_p_x, gam_p)
sy_p = sigma_at_IP(eps_p_n_y, beta_p_y, gam_p)

I_e = beam_current_mA(N_e, f_coll);  I_p = beam_current_mA(N_p, f_coll)
sqs = sqrt_s(E_e, E_p)
L_m = luminosity_unequal(N_e, N_p, f_coll, sx_e, sy_e, sx_p, sy_p)
L_cm = L_m * CM2_PER_M2
xi_p = xi_proton(N_e, beta_p_x, beta_p_y, gam_p, sx_e, sy_e)
xi_e = xi_electron(N_p, beta_e_x, beta_e_y, gam_e, sx_p, sy_p)

print(f"  Inputs:  E_e={E_e} GeV, E_p={E_p} GeV, n_b={n_b}, f_coll={f_coll:.2e} Hz")
print(f"           N_e={N_e:.2e},  N_p={N_p:.2e}")
print(f"           eps_n^e={eps_e_n_x}/{eps_e_n_y} µm,  eps_n^p={eps_p_n_x}/{eps_p_n_y} µm")
print(f"           beta*_e={beta_e_x}/{beta_e_y} m,    beta*_p={beta_p_x}/{beta_p_y} m\n")

line("sigma*_e (h/v)",     (sx_e*1e6, sy_e*1e6), ("94.8", "8.6"),  "{:>5.2f}", "µm")
line("sigma*_p (h/v)",     (sx_p*1e6, sy_p*1e6), ("94.9", "8.6"),  "{:>5.2f}", "µm")
line("I_e",                I_e, "1555", "{:>7.2f}", "mA")
line("I_p",                I_p, "1000", "{:>7.2f}", "mA")
line("sqrt(s)",            sqs, "105",  "{:>7.2f}", "GeV")
line("L (cm^-2 s^-1)",     L_cm,"1.05e34","{:>9.3e}", "")
line("xi_p (h/v)",         xi_p, ("0.012","0.012"), "{:>7.4f}", "")
line("xi_e (h/v)",         xi_e, ("--","--"),       "{:>7.3f}", "  (not quoted)")

print()
print("  Consistency check on quoted I_e=1555 mA:")
N_e_implied = 1555e-3 / (ELEM * f_coll)
print(f"     N_e implied by I_e=1555 mA & f_coll=91 MHz :  {N_e_implied:.3e}")
print(f"     But Table 1 quotes N_e=1.72e11.  These differ by factor "
      f"{N_e/N_e_implied:.2f}.")
print( "     The EIC CDR design point is I_e=2.5 A (and I_p=1 A) for 10 GeV;"
       " 1555 mA in Table 1 looks like a typo")
print(f"     (with N_e=1.72e11 and f_coll=91 MHz, I_e = "
       f"{ELEM*1.72e11*9.1e7*1e3:.0f} mA = 2.5 A).")
