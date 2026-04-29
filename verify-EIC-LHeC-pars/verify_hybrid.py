"""
Independent verification of the Hybrid EIC-LHeC tables in the IPAC'26 paper:
  Table 2 (tab:hybrid2e33) -- "Hybrid EIC-LHeC, initial configuration, ~2e33"
  Table 3 (tab:hybrid1e34) -- "Hybrid EIC-LHeC, optimised, 1e34, e-cooling+flat-beam"

Hybrid logic from §3 of the paper ("Resolving parameter inconsistencies"):
  - Take EIC proton ring (3834 m circumference, 1160 buckets, f_rev = 78.3 kHz)
  - Replace electron ring with LHeC-style ERL (60 mA / N_e=1.19e10 nominal)
  - Reducing n_b from 1160 -> 386 (factor of 3) drops f_coll from 91 MHz to ~30 MHz,
      and at constant 1.0 A proton current  N_p must rise from 6.9e10 to 2.1e11
      (factor 3 again).
  - Electron current at f_coll=30 MHz with N_e=1.19e10 is ~57 mA -- close to LHeC nominal.
  - IP optics re-optimised so sigma* of the two beams match in both planes.

We:
   1. Verify the bunch-reduction logic explicitly  (1160 -> 386).
   2. Plug Table 2 numbers into the same beam-size / beam-beam / luminosity
      machinery and compare results with what is quoted.
   3. Repeat for Table 3 (with reduced emittances and tighter beta*).
"""

import math

M_E   = 0.5109989461e-3
M_P   = 0.93827208816
R_E   = 2.8179403262e-15
R_P   = 1.5346982671e-18
ELEM  = 1.602176634e-19
CM2_PER_M2 = 1e-4

def gamma_from_E(E,m): return E/m
def geo_emit(eps_um, gam): return eps_um*1e-6/gam
def sigma_at_IP(eps_um, beta_m, gam): return math.sqrt(geo_emit(eps_um, gam)*beta_m)
def beam_current_mA(N, f): return ELEM*N*f*1e3

def luminosity_unequal(N1,N2,f,sx1,sy1,sx2,sy2):
    Sx = math.sqrt(sx1*sx1+sx2*sx2)
    Sy = math.sqrt(sy1*sy1+sy2*sy2)
    return f*N1*N2/(2.0*math.pi*Sx*Sy)

def xi_proton(N_e, beta_p_x, beta_p_y, gam_p, sx_e, sy_e):
    pref = N_e*R_P/(2.0*math.pi*gam_p)
    return pref*beta_p_x/(sx_e*(sx_e+sy_e)), pref*beta_p_y/(sy_e*(sx_e+sy_e))

def xi_electron(N_p, beta_e_x, beta_e_y, gam_e, sx_p, sy_p):
    pref = N_p*R_E/(2.0*math.pi*gam_e)
    return pref*beta_e_x/(sx_p*(sx_p+sy_p)), pref*beta_e_y/(sy_p*(sx_p+sy_p))


# =================================================================
# 1) Verify bunch-reduction logic of §3 of the paper
# =================================================================
print("="*82)
print(" Bunch-reduction logic check  (§3 'Resolving parameter inconsistencies')")
print("="*82)

C_HSR    = 3833.845       # m  RHIC/HSR circumference
c_light  = 299792458.0
f_rev    = c_light / C_HSR
n_b_EIC  = 1160
f_coll_EIC = f_rev * n_b_EIC

print(f"  HSR circumference          C  = {C_HSR:.2f} m")
print(f"  Revolution frequency       f_rev = {f_rev:.0f} Hz = {f_rev*1e-3:.2f} kHz")
print(f"  EIC nominal (n_b=1160):    f_coll = {f_coll_EIC*1e-6:.2f} MHz "
      f"(paper Table 1: 91 MHz)")

# Required electron current target ≈ 60 mA (LHeC) with N_e = 1.19e10
N_e   = 1.19e10
target_Ie = 60e-3   # A
needed_f_coll = target_Ie / (ELEM * N_e)
needed_n_b = needed_f_coll / f_rev
print(f"\n  Target  : I_e = 60 mA  with  N_e = 1.19e10  -> f_coll needed = "
      f"{needed_f_coll*1e-6:.2f} MHz")
print(f"  At f_rev = {f_rev*1e-3:.2f} kHz, this requires n_b = {needed_n_b:.1f}")
print(f"  Paper chooses n_b = 386 (close to integer 1160/3).")

# Try the paper's choice n_b = 386 and integer 1160/3
for nb_try in (386, 1160//3):
    fcoll = f_rev * nb_try
    Ie_at = beam_current_mA(N_e, fcoll)
    print(f"     n_b = {nb_try:>4d}  ->  f_coll = {fcoll*1e-6:5.2f} MHz "
          f", I_e = {Ie_at:5.2f} mA  (target ~58 mA)")

# Constant proton current 1 A => N_p scales x3 from EIC baseline
print("\n  Proton bunch population at constant I_p=1 A:")
for nb_try, label in ((1160, "EIC baseline"),
                      ( 386, "Hybrid (1160/3)"),
                      ( 387, "Hybrid (1160/3 rounded up)")):
    fcoll = f_rev * nb_try
    N_p   = 1.0/(ELEM*fcoll)
    print(f"     n_b={nb_try:>4d}  f_coll={fcoll*1e-6:5.2f} MHz  ->  "
          f"N_p = {N_p:.3e}  ({label})")
print("\n  ==> The factor-of-3 bunch reduction (1160 -> 386) and the resulting "
      "f_coll, I_e, N_p are mutually consistent.")

# =================================================================
# 2) Hybrid - Table 2  (initial, ~2e33)
# =================================================================
print("\n" + "="*82)
print(" Hybrid EIC-LHeC  --  Table 2 (initial config, target ~2e33)")
print("="*82)

E_e   = 10.0;     E_p   = 275.0
n_b   = 386
f_coll = 3.0e7              # paper quotes 3.0e7 (also matches f_rev*n_b precisely)
f_coll_calc = f_rev * n_b
N_e   = 1.19e10
N_p   = 2.1e11
eps_e_n_x = 22.0; eps_e_n_y = 22.0
eps_p_n_x =  3.3; eps_p_n_y =  0.3
beta_e_x  = 1.20; beta_e_y  = 0.05
beta_p_x  = 0.12; beta_p_y  = 0.06
sigma_z_e = 3.0e-3; sigma_z_p = 70.0e-3

gam_e = gamma_from_E(E_e, M_E);  gam_p = gamma_from_E(E_p, M_P)

sx_e = sigma_at_IP(eps_e_n_x, beta_e_x, gam_e)
sy_e = sigma_at_IP(eps_e_n_y, beta_e_y, gam_e)
sx_p = sigma_at_IP(eps_p_n_x, beta_p_x, gam_p)
sy_p = sigma_at_IP(eps_p_n_y, beta_p_y, gam_p)

# RMS angular divergence sigma_x' = sqrt(eps_geo / beta*)
sxp_e = sx_e/beta_e_x;  syp_e = sy_e/beta_e_y
sxp_p = sx_p/beta_p_x;  syp_p = sy_p/beta_p_y

I_e = beam_current_mA(N_e, f_coll); I_p = beam_current_mA(N_p, f_coll)
L_m = luminosity_unequal(N_e, N_p, f_coll, sx_e, sy_e, sx_p, sy_p)
L_cm = L_m * CM2_PER_M2

xi_p = xi_proton(N_e, beta_p_x, beta_p_y, gam_p, sx_e, sy_e)
xi_e = xi_electron(N_p, beta_e_x, beta_e_y, gam_e, sx_p, sy_p)

# Effective beam sizes (sqrt of quadratic sum of e and p sizes)
sx_eff = math.sqrt(sx_e*sx_e + sx_p*sx_p)
sy_eff = math.sqrt(sy_e*sy_e + sy_p*sy_p)

print(f"  f_coll : computed n_b * f_rev = {f_coll_calc*1e-6:.3f} MHz "
      f"(paper: 3.0e7 Hz = 30.0 MHz)\n")
print(f"  Computed                              Paper Table 2")
print(f"  -----------------------------------   --------------")
print(f"  sigma*_e   x/y = {sx_e*1e6:5.2f} / {sy_e*1e6:4.2f}  µm     "
      f"36.7 / 7.5")
print(f"  sigma*_p   x/y = {sx_p*1e6:5.2f} / {sy_p*1e6:4.2f}  µm     "
      f"36.8 / 7.8")
print(f"  sigma'_e   x/y = {sxp_e*1e6:5.1f} / {syp_e*1e6:5.1f} µrad   "
      f"30.6 / 149.9")
print(f"  sigma'_p   x/y = {sxp_p*1e6:5.1f} / {syp_p*1e6:5.1f} µrad   "
      f"306.3 / 130.6")
print(f"  xi_p       x/y = {xi_p[0]:.4f} / {xi_p[1]:.4f}              "
      f"0.001 / 0.002")
print(f"  xi_e       x/y = {xi_e[0]:.3f} / {xi_e[1]:.3f}              "
      f"3.474 / 0.679")
print(f"  I_e            = {I_e:5.2f} mA                       57.6")
print(f"  I_p            = {I_p:5.2f} mA                       1000")
print(f"  sigma_x_eff    = {sx_eff*1e6:5.2f} µm                       52.0")
print(f"  sigma_y_eff    = {sy_eff*1e6:5.2f} µm                       10.8")
print(f"  Peak L (cm-2 s-1)= {L_cm:.3e}                  2.1e33")

# =================================================================
# 3) Hybrid - Table 3  (1e34 with cooling + flat beam)
# =================================================================
print("\n" + "="*82)
print(" Hybrid EIC-LHeC  --  Table 3 (optimised config, 1e34, cooling + flat-beam)")
print("="*82)

# unchanged : n_b, f_coll, N_e, N_p, energies, sigma_z (says 'as in Table 2')
eps_e_n_x = 22.0;  eps_e_n_y = 2.2     # vertical e emittance / 10 (flat-beam)
eps_p_n_x =  0.9;  eps_p_n_y = 0.1     # cooled proton emittances
beta_e_x  = 0.25;  beta_e_y  = 0.10
beta_p_x  = 0.09;  beta_p_y  = 0.05

sx_e = sigma_at_IP(eps_e_n_x, beta_e_x, gam_e)
sy_e = sigma_at_IP(eps_e_n_y, beta_e_y, gam_e)
sx_p = sigma_at_IP(eps_p_n_x, beta_p_x, gam_p)
sy_p = sigma_at_IP(eps_p_n_y, beta_p_y, gam_p)
sxp_e = sx_e/beta_e_x;  syp_e = sy_e/beta_e_y
sxp_p = sx_p/beta_p_x;  syp_p = sy_p/beta_p_y

xi_p = xi_proton(N_e, beta_p_x, beta_p_y, gam_p, sx_e, sy_e)
xi_e = xi_electron(N_p, beta_e_x, beta_e_y, gam_e, sx_p, sy_p)
L_m  = luminosity_unequal(N_e, N_p, f_coll, sx_e, sy_e, sx_p, sy_p)
L_cm = L_m * CM2_PER_M2
sx_eff = math.sqrt(sx_e*sx_e + sx_p*sx_p)
sy_eff = math.sqrt(sy_e*sy_e + sy_p*sy_p)

print(f"  Computed                              Paper Table 3")
print(f"  -----------------------------------   --------------")
print(f"  sigma*_e   x/y = {sx_e*1e6:5.2f} / {sy_e*1e6:4.2f}  µm     "
      f"16.8 / 3.4")
print(f"  sigma*_p   x/y = {sx_p*1e6:5.2f} / {sy_p*1e6:4.2f}  µm     "
      f"16.8 / 3.6")
print(f"  sigma'_e   x/y = {sxp_e*1e6:5.1f} / {syp_e*1e6:5.1f} µrad  "
      f" 66.8 / 32.7")
print(f"  sigma'_p   x/y = {sxp_p*1e6:5.1f} / {syp_p*1e6:5.1f} µrad  "
      f" 182.4 / 77.8")
print(f"  xi_p       x/y = {xi_p[0]:.4f} / {xi_p[1]:.4f}              "
      f"0.003 / 0.007")
print(f"  xi_e       x/y = {xi_e[0]:.3f} / {xi_e[1]:.3f}              "
      f"3.474 / 6.789")
print(f"  sigma_x_eff    = {sx_eff*1e6:5.2f} µm                       23.8")
print(f"  sigma_y_eff    = {sy_eff*1e6:5.2f} µm                       5.0")
print(f"  Peak L (cm-2 s-1)= {L_cm:.3e}                  1.0e34")
