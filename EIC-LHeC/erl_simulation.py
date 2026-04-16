#!/usr/bin/env python3
"""
erl_simulation.py  –  Energy Recovery Linac Comprehensive Simulation
=====================================================================
Models single-pass and two-pass Energy Recovery Linacs (ERLs) for
electron-proton collisions at EIC-like parameters.

Beam parameters
---------------
  Bunch charge       : 1.5 nC
  Bunch length (rms) : 1.5 mm
  Injection energy   : 7 MeV
  Single-pass gain   : 10 GeV  (one linac traversal)
  Two-pass gain      : 5 GeV per pass  (same 10 GeV total)
  RF frequency       : 801.58 MHz
  Bunch rep rate     : 40.08 MHz  →  I_arc = 60 mA,  I_linac = 120 mA

Physics covered
---------------
  * Longitudinal phase-space evolution (macro-particle tracking)
  * 180-degree phase reversal in the return arc for deceleration
  * Electron-proton beam-beam interaction at the collision point
  * Energy recovery efficiency vs. deceleration-phase error
  * SRF cavity power budget (dynamic + cryogenic losses)
  * Synchrotron radiation in the return arcs
  * Wall-plug power comparison: ERL vs. conventional storage ring
  * Two-pass topology: linac operates at 5 GeV/pass (half-gradient)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Arc
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 – Physical constants
# ─────────────────────────────────────────────────────────────────────────────
c_light  = 2.99792458e8      # m/s
e_charge = 1.60217663e-19    # C
m_e_MeV  = 0.51099895        # MeV/c²  (electron rest mass)
m_p_MeV  = 938.27208816      # MeV/c²  (proton rest mass)
r_e      = 2.8179403227e-15  # m       (classical electron radius)
r_p      = 1.5346982e-18     # m       (classical proton radius)
C_gamma  = 8.85e-5           # m/GeV³  (synchrotron-radiation constant)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 – ERL / RF parameters
# ─────────────────────────────────────────────────────────────────────────────
Q_bunch   = 1.5e-9          # C        bunch charge
sigma_z   = 1.5e-3          # m        rms bunch length
E_inj     = 7.0             # MeV      injection energy
N_e       = Q_bunch / e_charge          # electrons per bunch ≈ 9.375e9

f_RF      = 801.58e6        # Hz       RF frequency
f_rep     = 40.08e6         # Hz       bunch repetition rate
omega_RF  = 2.0 * np.pi * f_RF
lambda_RF = c_light / f_RF  # m ≈ 374 mm
T_rep     = 1.0 / f_rep     # s        bunch spacing ≈ 24.95 ns
harmonic  = round(f_RF / f_rep)         # ≈ 20

I_arc    = Q_bunch * f_rep  # A        average arc current  = 60 mA
I_linac  = 2.0 * I_arc      # A        linac current (acc+dec) = 120 mA

# Single-pass ERL
E_gain_sp  = 10_000.0       # MeV      linac gain (single-pass)
E_peak_sp  = E_inj + E_gain_sp          # MeV ≈ 10007

# Two-pass ERL  (same total energy, but the linac voltage is halved)
E_gain_pp  = 5_000.0        # MeV      gain per pass
E_mid_tp   = E_inj + E_gain_pp          # MeV after pass 1
E_peak_tp  = E_mid_tp + E_gain_pp       # MeV after pass 2 ≈ 10007

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 – SRF linac parameters  (802 MHz, 5-cell cavity)
# ─────────────────────────────────────────────────────────────────────────────
G_acc        = 15.0     # MV/m     accelerating gradient
L_cell       = 0.187    # m        λ/2 cell length at 801.58 MHz
N_cells_cav  = 5        # cells per cavity
L_cav        = N_cells_cav * L_cell   # m ≈ 0.935 m  active length per cavity
V_cav        = G_acc * L_cav          # MV/cavity ≈ 14 MV

Q_0          = 2e10     # intrinsic quality factor
R_over_Q     = 400.0    # Ω   (linac convention, 5-cell at 802 MHz)
eta_cryo     = 250.0    # W_wall-plug / W_heat @ 2 K  (cryoplant COP⁻¹)
eta_RFPA     = 0.65     # RF power-amplifier efficiency
eta_inj      = 0.35     # injector (photogun + bunching) wall-plug efficiency

# Number of cavities
N_cav_sp  = int(np.ceil(E_gain_sp / V_cav))    # ≈ 715 for 10 GeV
N_cav_pp  = int(np.ceil(E_gain_pp / V_cav))    # ≈ 358 for 5 GeV

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 – Proton beam (EIC-like) and IP optics
# ─────────────────────────────────────────────────────────────────────────────
E_p        = 275.0e3        # MeV   proton energy = 275 GeV
gamma_p    = E_p / m_p_MeV  # ≈ 293
N_p        = 6.88e10        # protons per bunch
eps_n_p    = 1.0e-6         # m·rad normalised proton emittance
beta_star  = 0.05           # m     5 cm β* at IP (both planes, round)

# Electron normalised emittance at collision (cathode-limited, not SR equilibrium)
eps_n_e_sp = 20e-6          # m·rad  single-pass ERL (larger, single-pass wash-out)
eps_n_e_tp = 10e-6          # m·rad  two-pass ERL  (slightly better control)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 – Helper / physics functions
# ─────────────────────────────────────────────────────────────────────────────

def sigma_IP(eps_n, beta_s, gamma):
    """RMS beam size at the IP (round beam): σ = √(ε_n β* / γ)."""
    return np.sqrt(eps_n * beta_s / gamma)


def beam_beam_xi(N_opp, r_class, gamma_self, sigma):
    """
    Beam-beam tune-shift parameter (round, strong-weak):
        ξ = r · N / (4π γ σ²)
    """
    return r_class * N_opp / (4.0 * np.pi * gamma_self * sigma**2)


def luminosity(N1, N2, f, sigma1, sigma2, H=1.0):
    """
    Geometric luminosity  L = H N1 N2 f / (4π σ1 σ2).
    H includes hourglass and crossing-angle reduction.
    """
    return H * N1 * N2 * f / (4.0 * np.pi * sigma1 * sigma2)


def hourglass(sigma_z_m, beta_s):
    """Approximate hourglass factor  F ≈ 1 / √(1 + (σ_z/β*)²)."""
    return 1.0 / np.sqrt(1.0 + (sigma_z_m / beta_s) ** 2)


def disruption(N_opp, r_class, sz, gamma_self, sigma_T):
    """
    Disruption parameter  D = 2 r N σ_z / (γ σ_T²).
    Measures beam-focusing during the collision.
    """
    return 2.0 * r_class * N_opp * sz / (gamma_self * sigma_T**2)


def SR_loss_per_turn(E_GeV, rho_m):
    """
    Synchrotron-radiation energy loss per revolution [MeV]:
        U₀ = C_γ E⁴ / ρ    (C_γ = 8.85×10⁻⁵ m/GeV³)
    """
    return C_gamma * E_GeV**4 / rho_m * 1.0e3   # keV → MeV


def SRF_dynamic_2K(V_cav_MV, Rq, Q0, N_cav):
    """
    Total SRF dynamic heat load at 2 K [W]:
        P = V² / (R/Q · Q₀)  per cavity
    """
    return N_cav * (V_cav_MV * 1e6)**2 / (Rq * Q0)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 – Macro-particle longitudinal tracking
# ─────────────────────────────────────────────────────────────────────────────

def track_bunch(n_part=60_000, n_passes=1, E_gain=10_000.0,
                phi_acc_deg=0.0, delta_in_rms=4e-3, seed=42):
    """
    Track a macro-particle ensemble through the full ERL cycle:
        injection → n_passes acceleration → collision → 180° phase flip
        → n_passes deceleration → dump

    Coordinates: φ [rad] (RF phase relative to crest), δ = ΔE/E_ref
    Uses the thin-lens linac kick:  ΔE = eV₀[cos(φ_s + φ) − cos(φ_s)]

    Returns a dict with (φ, δ) snapshots at each stage.
    """
    rng  = np.random.default_rng(seed)
    snaps = {}

    phi_acc = np.deg2rad(phi_acc_deg)
    V_pass  = E_gain / n_passes              # MeV per acceleration pass

    # Initial distribution at injection
    sigma_phi = omega_RF * sigma_z / c_light  # rad, rms RF phase extent
    phi  = rng.normal(0.0, sigma_phi,  n_part)
    dE   = rng.normal(0.0, delta_in_rms * E_inj, n_part)  # MeV
    E_ref = E_inj                              # MeV, reference energy
    snaps["injection"] = (phi.copy(), dE.copy(), E_ref)

    # ── Acceleration passes ────────────────────────────────────────────────
    for k in range(n_passes):
        # Energy kick of each macro-particle
        kick = V_pass * (np.cos(phi_acc + phi) - np.cos(phi_acc))
        dE  += kick
        E_ref += V_pass * np.cos(phi_acc)
        snaps[f"acc_pass_{k+1}"] = (phi.copy(), dE.copy(), E_ref)

    snaps["post_accel"] = snaps[f"acc_pass_{n_passes}"]

    # ── Electron-proton collision (beam-beam kick) ─────────────────────────
    # A head-on e-p collision imparts a transverse kick that couples into
    # the longitudinal plane via nonlinear chromaticity and Oide effect.
    # Here we use a simplified linearised model: the longitudinal-phase
    # smearing from beamstrahlung / disruption.
    gamma_e_coll = E_ref / m_e_MeV
    sigma_e_coll = sigma_IP(eps_n_e_sp, beta_star, gamma_e_coll)
    D_e = disruption(N_p, r_e, sigma_z, gamma_e_coll, sigma_e_coll)
    xi_e = beam_beam_xi(N_p, r_e, gamma_e_coll, sigma_e_coll)
    # Energy spread added by beamstrahlung (Yokoya-Chen parametrisation, simplified)
    # δ_bs ≈ ξ_e × (σ_z / β*)  (very rough linear model)
    delta_bb = xi_e * (sigma_z / beta_star) * E_ref
    dE += rng.normal(0.0, delta_bb * 0.5, n_part)   # symmetric smearing
    snaps["post_collision"] = (phi.copy(), dE.copy(), E_ref)

    # ── 180-degree phase reversal in the return arc ────────────────────────
    # The arc transport is chosen so that the RF phase seen by the bunch
    # shifts by π, turning the accelerating crest into the decelerating trough.
    # For a single-pass ERL this is a single return arc with R_56 tuned to
    # produce exactly a π phase advance at the linac entrance.
    phi_dec = phi_acc + np.pi          # synchronous decelerating phase

    # ── Deceleration passes ────────────────────────────────────────────────
    for k in range(n_passes):
        kick = V_pass * (np.cos(phi_dec + phi) - np.cos(phi_dec))
        dE  += kick
        E_ref += V_pass * np.cos(phi_dec)   # reference loses energy (cos π = −1)
        snaps[f"dec_pass_{k+1}"] = (phi.copy(), dE.copy(), E_ref)

    snaps["dump"] = snaps[f"dec_pass_{n_passes}"]
    snaps["E_dump_ref"] = E_ref          # should be ≈ E_inj for ideal recovery
    return snaps


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 – Energy-recovery efficiency vs. phase error
# ─────────────────────────────────────────────────────────────────────────────

def recovery_efficiency(phi_err_deg, n_passes=1, E_gain=10_000.0):
    """
    Analytical energy recovery efficiency when the deceleration phase
    deviates by phi_err_deg from the ideal 180-degree reversal.

        E_recovered = E_gain × |cos(π + φ_err)|
        η_rec = E_recovered / E_gain
        E_dump = E_peak − E_recovered  (ideal crest acceleration)

    For n_passes each of E_gain/n_passes:
        Apply the same phase error on every deceleration pass.
    """
    phi_err = np.deg2rad(phi_err_deg)
    E_gain_per = E_gain / n_passes
    E = E_inj + E_gain          # peak energy (acc at crest)

    for _ in range(n_passes):
        # deceleration phase: π + error
        E -= E_gain_per * np.cos(np.pi + phi_err)   # cos(π+ε) = −cos ε
    E_dump = max(E, 0.001)

    E_recovered = (E_inj + E_gain) - E_dump
    eta = E_recovered / E_gain   # fraction of linac energy recovered
    return eta, E_dump


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8 – Power budgets
# ─────────────────────────────────────────────────────────────────────────────

def power_budget_ERL(n_passes=1, E_gain=10_000.0, rho_arc=300.0):
    """
    Wall-plug power budget for an ERL.

    Components
    ----------
    1. Beam-dump power  : I × E_dump  (residual kinetic energy at dump)
    2. Injector RF      : I × E_inj / η_inj
    3. SRF cryo         : P_2K × η_cryo  (dynamic cavity losses)
    4. SR in return arc : I × U₀(E_peak, ρ)
    5. Residual RF      : compensates imperfect beam-loading cancellation
    """
    bgt = {}
    E_peak  = E_inj + E_gain              # MeV
    E_peak_GeV = E_peak / 1e3

    # 1. Beam-dump power (ideal recovery → dump at E_inj)
    P_dump = I_arc * E_inj * 1e6 * 1e-6  # MW
    bgt["P_dump_MW"]    = P_dump

    # 2. Injector RF wall-plug
    P_inj_plug = I_arc * E_inj * 1e6 / eta_inj * 1e-6   # MW
    bgt["P_injector_MW"] = P_inj_plug

    # 3. SRF dynamic losses + cryoplant
    N_cav = int(np.ceil(E_gain / n_passes / V_cav))  # cavities (one linac segment)
    P_2K  = SRF_dynamic_2K(V_cav, R_over_Q, Q_0, N_cav)  # W at 2 K
    P_cryo_plug = P_2K * eta_cryo * 1e-6               # MW
    bgt["P_srf_2K_W"]   = P_2K
    bgt["P_cryo_MW"]    = P_cryo_plug
    bgt["N_cav"]        = N_cav

    # 4. Synchrotron radiation in the high-energy return arc
    #    (radiation in injection / dump lines at 7 MeV is negligible)
    U0_arc = SR_loss_per_turn(E_peak_GeV, rho_arc)      # MeV per turn
    #    A single-turn arc (not a ring), so SR loss ≈ U0 × arc_fraction
    #    We approximate as a full quarter-turn (θ ≈ π/2):
    arc_fraction = 0.5                                   # half of a full ring
    P_SR_arc = U0_arc * arc_fraction * 1e6 * I_arc * 1e-6  # MW
    bgt["P_SR_arc_MW"]  = P_SR_arc

    # 5. Residual RF (beam-loading cancellation is not perfect)
    #    Residual ≈ 0.1% of full beam power, amortised by η_RFPA
    P_beam_full = I_arc * E_peak * 1e6 * 1e-6           # MW
    P_rf_resid  = 1e-3 * P_beam_full / eta_RFPA
    bgt["P_RF_residual_MW"] = P_rf_resid

    # Totals
    P_total = (P_dump + P_inj_plug + P_cryo_plug
               + P_SR_arc + P_rf_resid)
    bgt["P_total_MW"]  = P_total
    bgt["eta_recovery"] = E_gain / (E_inj + E_gain)     # ~99.93%

    # For context: power needed WITHOUT energy recovery
    P_no_recovery = I_arc * E_peak * 1e6 * 1e-6         # MW
    bgt["P_no_recovery_MW"] = P_no_recovery
    bgt["savings_factor"]   = P_no_recovery / P_total

    return bgt


def power_budget_ring(E_GeV=10.007, I_A=0.060, rho_m=300.0, eta_RF=0.50):
    """
    Wall-plug power budget for a conventional electron storage ring
    at the same energy and current.

    Components
    ----------
    1. RF to compensate SR losses (dominant)
    2. Normal-conducting arc magnets
    3. Injector (replacement rate ~1%)
    """
    bgt = {}

    U0  = SR_loss_per_turn(E_GeV, rho_m)          # MeV/turn
    P_SR = U0 * 1e6 * I_A * 1e-6                  # MW (beam power)
    P_RF = P_SR / eta_RF                           # wall-plug RF
    bgt["U0_MeV"]           = U0
    bgt["P_SR_beam_MW"]     = P_SR
    bgt["P_RF_wallplug_MW"] = P_RF

    # Arc magnet power: rough ~1 kW/m of dipole for a 10-GeV ring
    C_ring = 2.0 * np.pi * rho_m * 4              # circumference (4-bend)
    P_mag  = 1e-3 * C_ring * 1e-3                 # MW (1 kW/m × C metres)
    bgt["P_magnets_MW"] = P_mag

    # Injector (top-up): beam lifetime in a 10-GeV ring is ~2 hours
    # Average replacement current = I / (lifetime_s × f_rep) × ...
    # Simplified: P_inj = I × E / (tau_beam × f_rep) / eta
    # More directly: P_inj [MW] = I_A × E_GeV×1e9 / (tau_s × eta_inj) × 1e-6
    tau_s = 2.0 * 3600.0          # 2-hour beam lifetime (realistic for 10 GeV ring)
    P_inj = I_A * E_GeV * 1e9 / (tau_s * eta_inj) * 1e-6  # MW
    bgt["P_injector_MW"] = P_inj

    bgt["P_total_MW"] = P_RF + P_mag + P_inj
    return bgt


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9 – Beam-beam analysis
# ─────────────────────────────────────────────────────────────────────────────

def bb_analysis(E_e_MeV, eps_n_e, N_e_b, N_p_b,
                eps_n_p_val, beta_s, f_r, label=""):
    """Full beam-beam characterisation for one configuration."""
    bb = {"label": label}
    gamma_e = E_e_MeV / m_e_MeV
    gamma_p_val = E_p / m_p_MeV

    sig_e = sigma_IP(eps_n_e, beta_s, gamma_e)
    sig_p = sigma_IP(eps_n_p_val, beta_s, gamma_p_val)
    sig_c = np.sqrt(sig_e**2 + sig_p**2)        # convolved size

    bb["sigma_e_um"]  = sig_e * 1e6
    bb["sigma_p_um"]  = sig_p * 1e6
    bb["sigma_c_um"]  = sig_c * 1e6

    xi_e = beam_beam_xi(N_p_b, r_e, gamma_e, sig_c)
    xi_p = beam_beam_xi(N_e_b, r_p, gamma_p_val, sig_c)
    # NOTE: In a storage ring ξ_e must be < 0.05-0.1 to avoid multi-turn
    # beam-beam resonances.  In a single-pass ERL each electron traverses
    # the IP only once, so resonant build-up is absent and ξ >> 0.1 is
    # achievable (limited instead by disruption / emittance growth, not resonance).
    bb["xi_e"] = xi_e
    bb["xi_p"] = xi_p

    D_e = disruption(N_p_b, r_e, sigma_z, gamma_e, sig_e)
    D_p = disruption(N_e_b, r_p, sigma_z, gamma_p_val, sig_p)
    bb["D_e"] = D_e
    bb["D_p"] = D_p

    H = hourglass(sigma_z, beta_s)
    bb["H_hg"] = H

    # Pinch enhancement (empirical: H_pinch ≈ 1 + 0.5 D^0.45 for D < 10)
    H_pinch = 1.0 + 0.5 * D_e**0.45 if D_e > 1e-3 else 1.0
    bb["H_pinch"] = H_pinch

    L_geo = luminosity(N_e_b, N_p_b, f_r, sig_e, sig_p, H)
    L_pin = L_geo * H_pinch
    bb["L_geo_cm2s"]   = L_geo  * 1e-4
    bb["L_pinch_cm2s"] = L_pin  * 1e-4

    E_cm = np.sqrt(4 * E_e_MeV * E_p) / 1e3    # GeV
    bb["E_cm_GeV"] = E_cm
    return bb


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10 – Run all analyses
# ─────────────────────────────────────────────────────────────────────────────

print("=" * 72)
print("  ENERGY RECOVERY LINAC  –  COMPREHENSIVE SIMULATION")
print("=" * 72)

print(f"\nRF / timing")
print(f"  RF frequency          {f_RF/1e6:.2f} MHz")
print(f"  Bunch rep rate        {f_rep/1e6:.2f} MHz")
print(f"  RF wavelength         {lambda_RF*1e3:.2f} mm")
print(f"  Harmonic number       {harmonic}")
print(f"  Bunch spacing         {T_rep*1e9:.2f} ns")
print(f"  Bunch charge          {Q_bunch*1e9:.1f} nC  ({N_e:.3e} e⁻/bunch)")
print(f"  I_arc (arcs)          {I_arc*1e3:.2f} mA")
print(f"  I_linac (acc+dec)     {I_linac*1e3:.2f} mA")

# ── Macro-particle tracking ────────────────────────────────────────────────
print("\nTracking bunches (50 000 macro-particles each)…")
snap_sp = track_bunch(n_part=50_000, n_passes=1, E_gain=E_gain_sp)
snap_tp = track_bunch(n_part=50_000, n_passes=2, E_gain=E_gain_sp)

# ── Power budgets ──────────────────────────────────────────────────────────
bgt_sp  = power_budget_ERL(n_passes=1, E_gain=E_gain_sp)
bgt_tp  = power_budget_ERL(n_passes=2, E_gain=E_gain_sp)
# Three storage-ring scenarios for comparison
# (a) Large ring 300 m – low SR, same current: favourable to ring
# (b) Compact ring 50 m – realistic 10 GeV design, same current
# (c) Compact ring 50 m at 2.5 A – high luminosity target (EIC design point)
bgt_ring_lrg  = power_budget_ring(E_GeV=E_peak_sp/1e3, I_A=I_arc,    rho_m=300.0)
bgt_ring_cmp  = power_budget_ring(E_GeV=E_peak_sp/1e3, I_A=I_arc,    rho_m=50.0)
bgt_ring_hiL  = power_budget_ring(E_GeV=E_peak_sp/1e3, I_A=2.5,      rho_m=50.0)
bgt_ring      = bgt_ring_cmp   # default for backward-compat references

# ── Beam-beam analysis ─────────────────────────────────────────────────────
bb_sp = bb_analysis(E_peak_sp, eps_n_e_sp, N_e, N_p, eps_n_p,
                    beta_star, f_rep, "1-pass ERL")
bb_tp = bb_analysis(E_peak_tp, eps_n_e_tp, N_e, N_p, eps_n_p,
                    beta_star, f_rep, "2-pass ERL")

# ── Phase-error efficiency sweep ───────────────────────────────────────────
phi_err_arr = np.linspace(-15, 15, 1001)          # degrees
eta_sp_arr   = np.array([recovery_efficiency(p, 1, E_gain_sp)[0] for p in phi_err_arr])
eta_tp_arr   = np.array([recovery_efficiency(p, 2, E_gain_sp)[0] for p in phi_err_arr])
Edump_sp_arr = np.array([recovery_efficiency(p, 1, E_gain_sp)[1] for p in phi_err_arr])
Edump_tp_arr = np.array([recovery_efficiency(p, 2, E_gain_sp)[1] for p in phi_err_arr])

# ─────────────────────────────────────────────────────────────────────────────
# Console summary
# ─────────────────────────────────────────────────────────────────────────────
for label, bgt, bb, snap, npasses in [
        ("SINGLE-PASS ERL (10 GeV gain)",  bgt_sp, bb_sp, snap_sp, 1),
        ("TWO-PASS ERL    (5 GeV / pass)", bgt_tp, bb_tp, snap_tp, 2)]:
    print(f"\n{'─'*65}")
    print(f"  {label}")
    print(f"{'─'*65}")
    print(f"  Injection energy         {E_inj:.1f} MeV")
    if npasses == 2:
        print(f"  Energy after pass 1      {E_mid_tp/1e3:.3f} GeV")
    print(f"  Peak (collision) energy  {E_peak_sp/1e3:.3f} GeV")
    print(f"  Dump energy (ideal)      {E_inj:.1f} MeV")
    print(f"  Energy recovery η        {bgt['eta_recovery']*100:.4f} %")
    print(f"  Savings vs. no recovery  {bgt['savings_factor']:.0f}×")
    print(f"\n  ─ Power budget ─")
    print(f"  Beam dump               {bgt['P_dump_MW']*1e3:7.1f} kW")
    print(f"  Injector (RF plug)      {bgt['P_injector_MW']*1e3:7.1f} kW")
    print(f"  SRF dyn. heat @ 2 K     {bgt['P_srf_2K_W']:7.1f} W  ({bgt['N_cav']} cavities)")
    print(f"  Cryoplant (wall-plug)   {bgt['P_cryo_MW']*1e3:7.1f} kW")
    print(f"  SR in return arc        {bgt['P_SR_arc_MW']*1e3:7.1f} kW")
    print(f"  Residual RF             {bgt['P_RF_residual_MW']*1e3:7.2f} kW")
    print(f"  TOTAL wall-plug         {bgt['P_total_MW']:7.3f} MW")
    print(f"  [Without ERL: would need {bgt['P_no_recovery_MW']:.0f} MW]")
    print(f"\n  ─ Beam-beam at IP ─")
    print(f"  Centre-of-mass energy   {bb['E_cm_GeV']:.1f} GeV")
    print(f"  σ_e at IP               {bb['sigma_e_um']:.2f} μm")
    print(f"  σ_p at IP               {bb['sigma_p_um']:.2f} μm")
    print(f"  Beam-beam param ξ_e     {bb['xi_e']:.4f}")
    print(f"  Beam-beam param ξ_p     {bb['xi_p']:.6f}")
    print(f"  Disruption D_e          {bb['D_e']:.3f}")
    print(f"  Hourglass H             {bb['H_hg']:.4f}")
    print(f"  Pinch enhancement       {bb['H_pinch']:.3f}")
    print(f"  L (geometric)           {bb['L_geo_cm2s']:.2e} cm⁻²s⁻¹")
    print(f"  L (with pinch)          {bb['L_pinch_cm2s']:.2e} cm⁻²s⁻¹")

print(f"\n{'─'*65}")
print(f"  STORAGE RING COMPARISON")
print(f"{'─'*65}")
for ring_lbl, bgt_r, I_r, rho_r in [
        ("Large ring  (ρ=300 m, 60 mA)",    bgt_ring_lrg, I_arc, 300),
        ("Compact ring (ρ= 50 m, 60 mA)",   bgt_ring_cmp, I_arc,  50),
        ("Compact ring (ρ= 50 m, 2.5 A ← high-lumi)", bgt_ring_hiL, 2.5, 50)]:
    print(f"\n  {ring_lbl}")
    print(f"  SR loss / turn (U₀)    {bgt_r['U0_MeV']*1e3:.1f} keV")
    print(f"  SR beam power          {bgt_r['P_SR_beam_MW']*1e3:.1f} kW")
    print(f"  RF wall-plug           {bgt_r['P_RF_wallplug_MW']:.3f} MW")
    print(f"  Arc magnets            {bgt_r['P_magnets_MW']:.3f} MW")
    print(f"  Injector               {bgt_r['P_injector_MW']:.3f} MW")
    print(f"  TOTAL wall-plug        {bgt_r['P_total_MW']:.3f} MW")
print(f"""
  ─ Comparison at equal current (60 mA) ─
  Note: at 60 mA the ERL injector dominates its budget (1.2 MW, η_inj=35%).
  The ERL's primary advantage is luminosity per unit power, not raw power.

  ─ At high luminosity (2.5 A storage ring scenario) ─
  Ring (2.5 A, ρ=50 m) total: {bgt_ring_hiL['P_total_MW']:.1f} MW
  ERL (1-pass,  60 mA)  total: {bgt_sp['P_total_MW']:.3f} MW
  ERL power advantage vs hi-lum ring: {bgt_ring_hiL['P_total_MW']/bgt_sp['P_total_MW']:.0f}×

  ─ Savings vs. operating WITHOUT energy recovery ─
  No-ERL power (1-pass design):  {bgt_sp['P_no_recovery_MW']:.0f} MW
  With ERL (1-pass):             {bgt_sp['P_total_MW']:.3f} MW
  Savings:                       {bgt_sp['savings_factor']:.0f}×
""")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11 – Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

BLUE   = "#1f77b4"
RED    = "#d62728"
GREEN  = "#2ca02c"
ORANGE = "#ff7f0e"
PURPLE = "#9467bd"
CYAN   = "#17becf"
GREY   = "#7f7f7f"

def set_panel(ax, title="", xlabel="", ylabel="", grid=True, fs=11):
    if title:  ax.set_title(title, fontsize=fs, fontweight="bold", pad=4)
    if xlabel: ax.set_xlabel(xlabel, fontsize=fs-1)
    if ylabel: ax.set_ylabel(ylabel, fontsize=fs-1)
    if grid:   ax.grid(True, ls="--", alpha=0.35, lw=0.7)
    ax.tick_params(labelsize=9)


def phase_to_mm(phi_arr):
    """Convert RF phase [rad] to longitudinal position [mm]."""
    return phi_arr / omega_RF * c_light * 1e3   # mm


def delta_to_keV(dE_arr, E_ref_MeV):
    return dE_arr * 1e3   # dE is already in MeV → keV


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 – ERL Schematic and Energy Profile (single-pass)
# ─────────────────────────────────────────────────────────────────────────────

fig1, axes1 = plt.subplots(2, 1, figsize=(13, 9),
                            gridspec_kw={"height_ratios": [1.6, 1.0]})
fig1.suptitle("Single-Pass Energy Recovery Linac – Lattice & Energy Profile",
              fontsize=13, fontweight="bold", y=0.98)

# ── Top panel: schematic ─────────────────────────────────────────────────
ax_sch = axes1[0]
ax_sch.set_xlim(0, 14)
ax_sch.set_ylim(-2.8, 2.8)
ax_sch.axis("off")

def draw_box(ax, x, y, w, h, fc, ec, lw=1.5, label="", fs=9, va="center", zorder=3):
    rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle="round,pad=0.06", fc=fc, ec=ec, lw=lw, zorder=zorder)
    ax.add_patch(rect)
    if label:
        ax.text(x, y, label, ha="center", va=va, fontsize=fs,
                fontweight="bold", zorder=zorder+1)

def arrow(ax, x0, y0, x1, y1, color=GREY, lw=1.8, style="->"):
    ax_sch.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle=style, color=color, lw=lw))

# Components
draw_box(ax_sch,  1.2, 0.0, 1.6, 0.9, "#d0e8ff", BLUE,  label="INJECTOR\n7 MeV")
draw_box(ax_sch,  4.5, 1.6, 4.2, 0.85,"#d0ffe8", GREEN, label="SRF LINAC   10 GeV gain")
draw_box(ax_sch,  4.5,-1.6, 4.2, 0.85,"#ffe8d0", RED,   label="SRF LINAC   (deceleration)")
draw_box(ax_sch,  9.5, 0.0, 1.0, 0.9, "#fff0d0", ORANGE,label="IP\ncollision", fs=8)
draw_box(ax_sch, 12.5, 0.0, 1.4, 0.9, "#ffd0d0", RED,   label="DUMP\n7 MeV")

# Return arcs (top = acceleration arc → IP)
arc_top = Arc((7.0, 0.0), 6.0, 3.5, angle=0, theta1=0, theta2=180,
              color=GREEN, lw=2.2)
ax_sch.add_patch(arc_top)

# Return arc from IP back through bottom linac to dump
arc_bot = Arc((6.0, 0.0), 7.0, 3.5, angle=0, theta1=180, theta2=360,
              color=RED, lw=2.2, linestyle="--")
ax_sch.add_patch(arc_bot)

# Arrows showing beam direction
arrow(ax_sch, 2.0, 0.5,  2.4, 1.55, GREEN, 2.0)
arrow(ax_sch, 6.6, 2.05, 8.8, 0.45, GREEN, 2.0)
arrow(ax_sch,10.0, 0.0, 11.8, 0.0, RED,   2.0)
arrow(ax_sch,10.0,-0.45, 8.8,-1.55, RED,   2.0, "->")
arrow(ax_sch, 2.4,-1.55, 2.0,-0.45, RED,   2.0)

# Energy labels along path
for txt, xp, yp, col in [
        ("7 MeV", 2.3, 2.25, GREEN),
        ("10.007 GeV", 7.0, 2.42, GREEN),
        ("10.007 GeV", 10.5, 0.55, ORANGE),
        ("10.007 GeV", 10.5,-0.55, RED),
        ("7 MeV",  2.3,-2.25, RED),
        ("7 MeV", 12.5, 0.60, RED)]:
    ax_sch.text(xp, yp, txt, ha="center", va="center",
                fontsize=8.5, color=col, fontweight="bold")

# Phase-flip annotation
ax_sch.annotate("180° phase\nreversal in arc",
                xy=(2.4, -1.0), xytext=(0.3, -2.3),
                fontsize=8, color="darkred",
                arrowprops=dict(arrowstyle="->", color="darkred", lw=1.2))

# Current labels
ax_sch.text(4.5, 2.35, f"I = {I_linac*1e3:.0f} mA  (acc + dec bunches)",
            ha="center", fontsize=8.5, color=GREY)
ax_sch.text(9.5, 1.45, f"I = {I_arc*1e3:.0f} mA",
            ha="center", fontsize=8.5, color=GREEN)
ax_sch.text(9.5,-1.45, f"I = {I_arc*1e3:.0f} mA",
            ha="center", fontsize=8.5, color=RED)
ax_sch.set_title("ERL Layout", fontsize=11, fontweight="bold", pad=2)

# ── Bottom panel: energy profile along ERL ─────────────────────────────
ax_e = axes1[1]

# Parametric position along the ERL path
#  0→1  injection line
#  1→2  linac (acceleration)
#  2→3  top arc to IP
#  3→4  IP region
#  4→5  bottom arc from IP
#  5→6  linac (deceleration)
#  6→7  dump line

pos_nodes = np.array([0, 1, 2, 3,     3.3,   4, 5, 6, 7])
E_nodes   = np.array([E_inj, E_inj, E_peak_sp, E_peak_sp,
                       E_peak_sp*0.9998,   # small loss at IP (beamstrahlung ~0.02%)
                       E_peak_sp, E_peak_sp, E_inj, E_inj])

# Smooth interpolation
from scipy.interpolate import PchipInterpolator
pos_fine = np.linspace(0, 7, 2000)
E_fine   = PchipInterpolator(pos_nodes, E_nodes)(pos_fine)

ax_e.plot(pos_fine, E_fine / 1e3, color=BLUE, lw=2.5)
ax_e.fill_between(pos_fine, E_inj/1e3, E_fine/1e3,
                  where=(E_fine > E_inj), alpha=0.15, color=BLUE, label="Energy stored in beam")

# Marks
markers = {
    "Injection\n7 MeV":  (0.0,  E_inj),
    "Linac exit\n10.007 GeV": (2.0, E_peak_sp),
    "IP\ncollision":     (3.15, E_peak_sp),
    "Dump\n7 MeV":       (7.0,  E_inj),
}
for lbl, (xp, ep) in markers.items():
    ax_e.scatter(xp, ep/1e3, zorder=5, s=60, color=RED)
    ax_e.annotate(lbl, (xp, ep/1e3), textcoords="offset points",
                  xytext=(0, 12), ha="center", fontsize=8, color="darkred")

# Shade the linac sections
for xs, xe, fc in [(1, 2, "#d0ffe8"), (5, 6, "#ffe8d0")]:
    ax_e.axvspan(xs, xe, alpha=0.25, color=fc)

ax_e.text(1.5, 5.5, "SRF Linac\n(accel.)", ha="center", fontsize=8, color=GREEN)
ax_e.text(5.5, 5.5, "SRF Linac\n(decel.)", ha="center", fontsize=8, color=RED)

# Annotate recovered energy
ax_e.annotate("", xy=(6.0, E_inj/1e3 + 0.4),
              xytext=(6.0, E_peak_sp/1e3 - 0.4),
              arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=2))
ax_e.text(6.35, E_peak_sp/2e3, f"Recovered\n{E_gain_sp/1e3:.0f} GeV",
          ha="left", fontsize=8.5, color=ORANGE)

xtick_labels = ["Cathode", "Linac\nentrance", "Linac\nexit",
                "IP", "", "Return\narc", "Linac\nentrance", "Linac\nexit", "Dump"]
ax_e.set_xticks(range(8))
ax_e.set_xticklabels(xtick_labels[:8], fontsize=8)
set_panel(ax_e, "Energy Profile Along the ERL Beam Path",
          "Position along lattice →", "Beam Energy (GeV)")
ax_e.set_ylim(-0.3, E_peak_sp/1e3 + 1.5)
ax_e.legend(fontsize=9, loc="center right")

fig1.tight_layout(rect=[0, 0, 1, 0.97])
fig1.savefig("ERL_fig1_layout_energy.png", dpi=150, bbox_inches="tight")
print("\nSaved: ERL_fig1_layout_energy.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 – Longitudinal Phase Space Evolution (single-pass)
# ─────────────────────────────────────────────────────────────────────────────

fig2, axes2 = plt.subplots(2, 4, figsize=(16, 7))
fig2.suptitle("Single-Pass ERL – Longitudinal Phase Space Evolution",
              fontsize=13, fontweight="bold")

stage_keys = ["injection", "acc_pass_1", "post_collision", "dump"]
stage_lbl  = ["(a) Injection  7 MeV",
               "(b) After linac  10.007 GeV",
               "(c) After e-p collision",
               "(d) After deceleration  7 MeV"]
colors_ps  = [BLUE, GREEN, ORANGE, RED]

phi_sp, dE_sp, Eref_sp = snap_sp["injection"]

for col, (key, lbl, clr) in enumerate(zip(stage_keys, stage_lbl, colors_ps)):
    phi_k, dE_k, Eref_k = snap_sp[key]
    z_mm   = phase_to_mm(phi_k)       # longitudinal position [mm]
    dE_MeV = dE_k                      # energy deviation [MeV]

    # Top row: z vs ΔE scatter
    ax_top = axes2[0, col]
    ax_top.scatter(z_mm, dE_MeV * 1e3, s=0.15, c=clr, alpha=0.25, rasterized=True)
    set_panel(ax_top, lbl, "z (mm)", "ΔE (keV)", fs=9)
    ax_top.set_xlim(-6, 6)

    # Bottom row: histogram of ΔE
    ax_bot = axes2[1, col]
    ax_bot.hist(dE_MeV * 1e3, bins=120, color=clr, alpha=0.7, density=True)
    rms = np.std(dE_MeV * 1e3)
    ax_bot.axvline(0, color="k", lw=1, ls="--")
    ax_bot.text(0.97, 0.93, f"σ = {rms:.1f} keV", transform=ax_bot.transAxes,
                ha="right", fontsize=8.5, color=clr,
                bbox=dict(fc="white", alpha=0.7, ec="none"))
    set_panel(ax_bot, "", "ΔE (keV)", "Density", fs=9)

fig2.tight_layout()
fig2.savefig("ERL_fig2_phase_space.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig2_phase_space.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 – RF Phase Relationship
# ─────────────────────────────────────────────────────────────────────────────

fig3, (ax_rf1, ax_rf2, ax_rf3) = plt.subplots(1, 3, figsize=(15, 5))
fig3.suptitle("RF Phase Relationship in the ERL", fontsize=13, fontweight="bold")

phi_rf = np.linspace(-np.pi, np.pi, 1000)
V_rf   = np.cos(phi_rf)            # normalised RF voltage

# ── Panel (a): Accelerating bucket ────────────────────────────────────────
ax_rf1.plot(np.degrees(phi_rf), V_rf, "b-", lw=2.2, label="RF voltage")
ax_rf1.axvline(0, color=GREEN, lw=2, ls="-", label="Sync. phase (crest, φ=0)")
ax_rf1.fill_between(np.degrees(phi_rf), V_rf,
                    where=(np.abs(phi_rf) < np.pi/2), alpha=0.15, color=GREEN)

# Show bunch spread
phi_bunch = np.linspace(-0.12, 0.12, 200)
ax_rf1.fill_between(np.degrees(phi_bunch), np.cos(phi_bunch),
                    color=BLUE, alpha=0.55, label="Bunch distribution")
ax_rf1.set_xlim(-180, 180)
ax_rf1.set_ylim(-1.25, 1.3)
set_panel(ax_rf1, "(a) Acceleration (crest, φ_s = 0°)",
          "RF Phase φ (deg)", "Normalised RF voltage")
ax_rf1.legend(fontsize=8.5)

# ── Panel (b): 180-degree shift → decelerating trough ─────────────────────
ax_rf2.plot(np.degrees(phi_rf), V_rf, "r--", lw=2.2, label="RF voltage (same cavity)")
ax_rf2.axvline(180, color=RED, lw=2, ls="-", label="Sync. phase (trough, φ=180°)")
ax_rf2.fill_between(np.degrees(phi_rf), V_rf,
                    where=(np.cos(phi_rf) < 0), alpha=0.15, color=RED)

phi_bunch_d = np.linspace(np.pi - 0.12, np.pi + 0.12, 200)
ax_rf2.fill_between(np.degrees(phi_bunch_d), np.cos(phi_bunch_d),
                    color=ORANGE, alpha=0.55, label="Returned bunch")
ax_rf2.set_xlim(0, 360)
ax_rf2.set_ylim(-1.25, 1.3)
set_panel(ax_rf2, "(b) Deceleration (trough, φ_s = 180°)\nafter return arc phase flip",
          "RF Phase φ (deg)", "Normalised RF voltage")
ax_rf2.legend(fontsize=8.5)

# ── Panel (c): Phasor diagram showing cancellation of beam loading ─────────
theta  = np.linspace(0, 2*np.pi, 400)
ax_rf3.plot(np.cos(theta), np.sin(theta), "k-", lw=0.8, alpha=0.3)

# Accelerating bunch phasor
ax_rf3.annotate("", xy=(1.0, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=2.5))
ax_rf3.text(1.05, 0.12, "Acc. beam\nloading\n(I_acc)", color=GREEN, fontsize=9)

# Decelerating bunch phasor (anti-phase → cancels)
ax_rf3.annotate("", xy=(-1.0, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color=RED, lw=2.5))
ax_rf3.text(-2.05, 0.12, "Dec. beam\nloading\n(I_dec)", color=RED, fontsize=9)

# Net phasor ≈ 0
ax_rf3.scatter(0, 0, s=120, color=ORANGE, zorder=5, label="Net beam loading ≈ 0")
ax_rf3.set_xlim(-2.3, 2.3)
ax_rf3.set_ylim(-1.6, 1.6)
ax_rf3.set_aspect("equal")
set_panel(ax_rf3, "(c) Beam-loading phasor cancellation\nin ERL linac",
          "Re", "Im")
ax_rf3.legend(fontsize=8.5, loc="lower center")
ax_rf3.text(0, -1.5,
            f"I_linac = {I_linac*1e3:.0f} mA  (I_acc + I_dec)\n"
            f"Net generator power ≈ 0  (ideal ERL)",
            ha="center", fontsize=9, color=GREY)

fig3.tight_layout()
fig3.savefig("ERL_fig3_RF_phase.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig3_RF_phase.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 – Energy Recovery Efficiency Analysis
# ─────────────────────────────────────────────────────────────────────────────

fig4, axes4 = plt.subplots(2, 2, figsize=(13, 9))
fig4.suptitle("Energy Recovery Efficiency Analysis", fontsize=13, fontweight="bold")

# ── (a) Recovery efficiency vs. phase error ────────────────────────────────
ax = axes4[0, 0]
ax.plot(phi_err_arr, eta_sp_arr * 100, color=BLUE,  lw=2.2, label="Single-pass ERL")
ax.plot(phi_err_arr, eta_tp_arr * 100, color=RED,   lw=2.2, ls="--", label="Two-pass ERL")
ax.axhline(99.93, color=GREY, lw=1, ls=":", label="Ideal (100%)")
ax.axvline(0, color="k", lw=0.8)
set_panel(ax, "(a) Energy Recovery Efficiency vs. Phase Error",
          "Deceleration phase error Δφ (deg)", "Recovery efficiency η (%)")
ax.set_ylim(60, 101)
ax.legend(fontsize=9)

# ── (b) Dump energy vs. phase error ───────────────────────────────────────
ax = axes4[0, 1]
ax.semilogy(phi_err_arr, Edump_sp_arr, color=BLUE,  lw=2.2, label="Single-pass")
ax.semilogy(phi_err_arr, Edump_tp_arr, color=RED,   lw=2.2, ls="--", label="Two-pass")
ax.axhline(E_inj, color=GREEN, lw=1.5, ls=":", label=f"E_inj = {E_inj} MeV (ideal)")
ax.axvline(0, color="k", lw=0.8)
set_panel(ax, "(b) Dump Energy vs. Phase Error",
          "Phase error Δφ (deg)", "Dump energy E_dump (MeV)")
ax.legend(fontsize=9)

# ── (c) Wall-plug power vs. phase error ────────────────────────────────────
ax = axes4[1, 0]
P_dump_sp = I_arc * Edump_sp_arr * 1e6 * 1e-3    # kW
P_dump_tp = I_arc * Edump_tp_arr * 1e6 * 1e-3    # kW
ax.plot(phi_err_arr, P_dump_sp, color=BLUE,  lw=2.2, label="Single-pass dump power")
ax.plot(phi_err_arr, P_dump_tp, color=RED,   lw=2.2, ls="--", label="Two-pass dump power")
ax.axhline(I_arc * E_inj * 1e6 * 1e-3, color=GREEN, lw=1.5, ls=":",
           label=f"Ideal: {I_arc*E_inj*1e6*1e-3:.0f} kW")
set_panel(ax, "(c) Beam Dump Power vs. Phase Error",
          "Phase error Δφ (deg)", "Dump power (kW)")
ax.legend(fontsize=9)

# ── (d) Energy-recovery gain relative to ideal ────────────────────────────
ax = axes4[1, 1]
# Compare wall-plug of ERL (with phase error) vs storage ring baseline
P_ring_tot = bgt_ring["P_total_MW"] * 1e3   # kW
P_erl_base = bgt_sp["P_total_MW"]   * 1e3   # kW (ideal)
P_err_sp   = P_dump_sp + (bgt_sp["P_total_MW"]*1e3 - bgt_sp["P_dump_MW"]*1e3)
P_err_tp   = P_dump_tp + (bgt_tp["P_total_MW"]*1e3 - bgt_tp["P_dump_MW"]*1e3)

ax.plot(phi_err_arr, P_err_sp, color=BLUE,  lw=2.2, label="Single-pass ERL")
ax.plot(phi_err_arr, P_err_tp, color=RED,   lw=2.2, ls="--", label="Two-pass ERL")
ax.axhline(P_ring_tot, color=ORANGE, lw=2.0, ls="-.",
           label=f"Storage ring: {P_ring_tot:.0f} kW")
ax.axhline(P_erl_base, color=GREEN, lw=1.5, ls=":",
           label=f"ERL ideal: {P_erl_base:.0f} kW")
set_panel(ax, "(d) Total ERL Wall-Plug Power vs. Phase Error",
          "Phase error Δφ (deg)", "Total wall-plug power (kW)")
ax.legend(fontsize=8.5)

fig4.tight_layout()
fig4.savefig("ERL_fig4_efficiency.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig4_efficiency.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5 – Power Budget Comparison
# ─────────────────────────────────────────────────────────────────────────────

fig5, axes5 = plt.subplots(1, 3, figsize=(15, 6))
fig5.suptitle("Wall-Plug Power Budget: ERL vs Storage Ring", fontsize=13, fontweight="bold")

def draw_stacked_bar(ax, configs, width=0.55):
    """configs: list of (label, ordered_dict of {component: value_kW})"""
    keys_all = set()
    for _, d in configs:
        keys_all |= set(d.keys())
    keys_all = list(keys_all)

    bar_colors = {"Beam dump":     BLUE,
                  "Injector RF":   CYAN,
                  "SRF cryoplant": GREEN,
                  "SR in arcs":    PURPLE,
                  "Residual RF":   GREY,
                  "Ring RF":       ORANGE,
                  "Arc magnets":   RED}
    x_positions = np.arange(len(configs))
    bottoms = np.zeros(len(configs))

    for comp in ["Beam dump", "Injector RF", "SRF cryoplant",
                 "SR in arcs", "Residual RF", "Ring RF", "Arc magnets"]:
        vals = []
        for _, d in configs:
            vals.append(d.get(comp, 0.0))
        bars = ax.bar(x_positions, vals, bottom=bottoms, width=width,
                      color=bar_colors.get(comp, GREY), label=comp,
                      edgecolor="white", lw=0.7)
        # Label if significant
        for i, (bar, v) in enumerate(zip(bars, vals)):
            if v > 30:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bottoms[i] + v/2, f"{v:.0f}", ha="center",
                        va="center", fontsize=8.5, color="white", fontweight="bold")
        bottoms += np.array(vals)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([c[0] for c in configs], fontsize=10, fontweight="bold")
    ax.set_ylabel("Power (kW)", fontsize=10)
    ax.legend(fontsize=8, loc="upper right", ncol=1)

# Build data dicts (kW)
data_sp = {
    "Beam dump":     bgt_sp["P_dump_MW"] * 1e3,
    "Injector RF":   bgt_sp["P_injector_MW"] * 1e3,
    "SRF cryoplant": bgt_sp["P_cryo_MW"] * 1e3,
    "SR in arcs":    bgt_sp["P_SR_arc_MW"] * 1e3,
    "Residual RF":   bgt_sp["P_RF_residual_MW"] * 1e3,
}
data_tp = {
    "Beam dump":     bgt_tp["P_dump_MW"] * 1e3,
    "Injector RF":   bgt_tp["P_injector_MW"] * 1e3,
    "SRF cryoplant": bgt_tp["P_cryo_MW"] * 1e3,
    "SR in arcs":    bgt_tp["P_SR_arc_MW"] * 1e3,
    "Residual RF":   bgt_tp["P_RF_residual_MW"] * 1e3,
}
def ring_data_dict(bgt_r):
    return {"Ring RF":    bgt_r["P_RF_wallplug_MW"] * 1e3,
            "Arc magnets":bgt_r["P_magnets_MW"]     * 1e3,
            "Injector RF":bgt_r["P_injector_MW"]    * 1e3}

data_ring_lrg  = ring_data_dict(bgt_ring_lrg)
data_ring_cmp  = ring_data_dict(bgt_ring_cmp)
data_ring_hiL  = ring_data_dict(bgt_ring_hiL)
data_ring      = data_ring_cmp  # backward compat

ax_bar = axes5[0]
draw_stacked_bar(ax_bar,
    [("1-pass\nERL", data_sp),
     ("2-pass\nERL", data_tp),
     ("Ring\n300 m\n60 mA", data_ring_lrg),
     ("Ring\n50 m\n60 mA",  data_ring_cmp),
     ("Ring\n50 m\n2.5 A",  data_ring_hiL)])
set_panel(ax_bar, "Wall-Plug Power Breakdown", grid=False)
ax_bar.set_title("(a) Component breakdown\n(3 ring scenarios shown)",
                  fontsize=10, fontweight="bold")

# ── (b) Pie chart – single-pass ERL ───────────────────────────────────────
ax_pie1 = axes5[1]
labels_pie = list(data_sp.keys())
vals_pie   = [data_sp[k] for k in labels_pie]
clrs_pie   = [BLUE, CYAN, GREEN, PURPLE, GREY]
wedges, texts, autotexts = ax_pie1.pie(
    vals_pie, labels=None, autopct="%1.1f%%",
    colors=clrs_pie, startangle=90,
    textprops={"fontsize": 8}, pctdistance=0.75,
    wedgeprops={"edgecolor": "white", "lw": 1.2})
ax_pie1.legend(wedges, labels_pie, fontsize=8, loc="lower center",
               bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_pie1.set_title(f"(b) Single-pass ERL\nTotal = {bgt_sp['P_total_MW']*1e3:.0f} kW",
                  fontsize=11, fontweight="bold")

# ── (c) Savings summary bar ────────────────────────────────────────────────
ax_sav = axes5[2]
categories  = ["No ERL\n(602 MW!)",
               "Ring 300 m\n60 mA",
               "Ring 50 m\n60 mA",
               "Ring 50 m\n2.5 A (hi-lum)",
               "2-pass ERL\n60 mA",
               "1-pass ERL\n60 mA"]
P_values_MW = [
    bgt_sp["P_no_recovery_MW"],
    bgt_ring_lrg["P_total_MW"],
    bgt_ring_cmp["P_total_MW"],
    bgt_ring_hiL["P_total_MW"],
    bgt_tp["P_total_MW"],
    bgt_sp["P_total_MW"],
]
bar_cols = [RED, GREY, GREY, ORANGE, PURPLE, BLUE]
yb = ax_sav.barh(categories, P_values_MW, color=bar_cols, edgecolor="white", lw=1.2,
                 height=0.6)
for bar, v in zip(yb, P_values_MW):
    ax_sav.text(v * 1.1, bar.get_y() + bar.get_height()/2,
                f"{v:.2f} MW", va="center", fontsize=8.5)
set_panel(ax_sav, "(c) Total Wall-Plug Power  [log scale]", "Power (MW)", grid=True)
ax_sav.set_xscale("log")
ax_sav.set_xlim(0.01, 2000)
# Annotate the savings on ERL bars vs no-ERL
for bar, v_this in [(yb[-2], P_values_MW[-2]), (yb[-1], P_values_MW[-1])]:
    saving = bgt_sp["P_no_recovery_MW"] / v_this
    ax_sav.text(0.02, bar.get_y() + bar.get_height()/2,
                f"  ×{saving:.0f} vs no ERL",
                va="center", fontsize=8, color="white", fontweight="bold")
# Annotate ERL advantage vs high-lum ring
for bar, v_this in [(yb[-2], P_values_MW[-2]), (yb[-1], P_values_MW[-1])]:
    adv = bgt_ring_hiL["P_total_MW"] / v_this
    ax_sav.text(bgt_ring_hiL["P_total_MW"] * 1.1,
                bar.get_y() + bar.get_height()/2,
                f"  ←{adv:.0f}× vs hi-lum ring",
                va="center", fontsize=7.5, color="darkred")

fig5.tight_layout()
fig5.savefig("ERL_fig5_power_budget.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig5_power_budget.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 6 – Beam-Beam Analysis
# ─────────────────────────────────────────────────────────────────────────────

fig6, axes6 = plt.subplots(2, 3, figsize=(16, 10))
fig6.suptitle("Electron-Proton Beam-Beam Analysis at the IP",
              fontsize=13, fontweight="bold")

# ── (a) Beam sizes vs electron energy ─────────────────────────────────────
ax = axes6[0, 0]
E_e_sweep = np.linspace(500, 12000, 400)   # MeV
gamma_e_sw = E_e_sweep / m_e_MeV
sig_e_sw1 = sigma_IP(eps_n_e_sp, beta_star, gamma_e_sw) * 1e6  # μm
sig_e_sw2 = sigma_IP(eps_n_e_tp, beta_star, gamma_e_sw) * 1e6
gamma_p_sw = gamma_p  # fixed proton energy
sig_p_sw   = sigma_IP(eps_n_p, beta_star, gamma_p_sw) * 1e6

ax.semilogy(E_e_sweep/1e3, sig_e_sw1, color=BLUE,  lw=2.2, label=f"e⁻ (ε_n={eps_n_e_sp*1e6:.0f} mm·mrad, 1-pass)")
ax.semilogy(E_e_sweep/1e3, sig_e_sw2, color=CYAN,  lw=2.2, ls="--", label=f"e⁻ (ε_n={eps_n_e_tp*1e6:.0f} mm·mrad, 2-pass)")
ax.axhline(sig_p_sw, color=RED, lw=1.8, ls=":", label=f"p (ε_n={eps_n_p*1e6:.0f} mm·mrad, fixed)")
ax.axvline(E_peak_sp/1e3, color=GREY, lw=1.2, ls="--", alpha=0.6)
ax.set_xlabel("Electron energy (GeV)", fontsize=10)
ax.set_ylabel("σ at IP (μm)", fontsize=10)
set_panel(ax, "(a) Beam Sizes at IP vs Electron Energy")
ax.legend(fontsize=8)

# ── (b) Beam-beam parameter ξ_e vs beta* ──────────────────────────────────
ax = axes6[0, 1]
bstar_sw = np.linspace(0.01, 0.30, 300)   # m
gamma_e_peak_sp = E_peak_sp / m_e_MeV
gamma_e_peak_tp = E_peak_tp / m_e_MeV

xi_sp_bstar = np.array([beam_beam_xi(N_p, r_e, gamma_e_peak_sp,
               sigma_IP(eps_n_e_sp, bs, gamma_e_peak_sp)) for bs in bstar_sw])
xi_tp_bstar = np.array([beam_beam_xi(N_p, r_e, gamma_e_peak_tp,
               sigma_IP(eps_n_e_tp, bs, gamma_e_peak_tp)) for bs in bstar_sw])

ax.loglog(bstar_sw * 100, xi_sp_bstar, color=BLUE, lw=2.2, label="1-pass ERL")
ax.loglog(bstar_sw * 100, xi_tp_bstar, color=RED,  lw=2.2, ls="--", label="2-pass ERL")
ax.axhline(0.05, color=ORANGE, lw=1.5, ls=":", label="Storage-ring limit ξ ≈ 0.05")
ax.axhline(0.10, color=ORANGE, lw=1.0, ls="--", alpha=0.5, label="ERL tolerance ξ ~ 0.1")
ax.axvline(beta_star * 100, color=GREY, lw=1.2, ls="--", alpha=0.7, label=f"β* = {beta_star*100:.0f} cm")
ax.set_xlabel("β* (cm)", fontsize=10)
ax.set_ylabel("Beam-beam parameter ξ_e", fontsize=10)
set_panel(ax, "(b) Electron Beam-Beam Parameter vs β*")
ax.legend(fontsize=8)

# ── (c) Luminosity vs beta* ────────────────────────────────────────────────
ax = axes6[0, 2]
L_sp_bstar = []
L_tp_bstar = []
for bs in bstar_sw:
    H = hourglass(sigma_z, bs)
    s_e1 = sigma_IP(eps_n_e_sp, bs, gamma_e_peak_sp)
    s_e2 = sigma_IP(eps_n_e_tp, bs, gamma_e_peak_tp)
    s_p  = sigma_IP(eps_n_p, bs, gamma_p)
    L_sp_bstar.append(luminosity(N_e, N_p, f_rep, s_e1, s_p, H) * 1e-4)
    L_tp_bstar.append(luminosity(N_e, N_p, f_rep, s_e2, s_p, H) * 1e-4)

ax.loglog(bstar_sw * 100, L_sp_bstar, color=BLUE, lw=2.2, label="1-pass ERL")
ax.loglog(bstar_sw * 100, L_tp_bstar, color=RED,  lw=2.2, ls="--", label="2-pass ERL")
ax.axvline(beta_star * 100, color=GREY, lw=1.2, ls="--", alpha=0.7,
           label=f"Design β* = {beta_star*100:.0f} cm")
ax.set_xlabel("β* (cm)", fontsize=10)
ax.set_ylabel("Luminosity (cm⁻²s⁻¹)", fontsize=10)
set_panel(ax, "(c) Luminosity vs β*")
ax.legend(fontsize=8)

# ── (d) Disruption parameter vs electron energy ────────────────────────────
ax = axes6[1, 0]
D_e_sw1 = np.array([disruption(N_p, r_e, sigma_z,
              Ee/m_e_MeV, sigma_IP(eps_n_e_sp, beta_star, Ee/m_e_MeV))
              for Ee in E_e_sweep])
D_e_sw2 = np.array([disruption(N_p, r_e, sigma_z,
              Ee/m_e_MeV, sigma_IP(eps_n_e_tp, beta_star, Ee/m_e_MeV))
              for Ee in E_e_sweep])

ax.semilogy(E_e_sweep/1e3, D_e_sw1, color=BLUE,  lw=2.2, label="1-pass ERL")
ax.semilogy(E_e_sweep/1e3, D_e_sw2, color=RED,   lw=2.2, ls="--", label="2-pass ERL")
ax.axhline(1.0, color=ORANGE, lw=1.5, ls=":", label="D = 1 (pinch regime)")
ax.axvline(E_peak_sp/1e3, color=GREY, lw=1.2, ls="--", alpha=0.6)
ax.set_xlabel("Electron energy (GeV)", fontsize=10)
ax.set_ylabel("Disruption parameter D_e", fontsize=10)
set_panel(ax, "(d) Electron Disruption Parameter")
ax.legend(fontsize=8)

# ── (e) ERL vs storage ring: luminosity at same current ───────────────────
ax = axes6[1, 1]
I_sweep = np.logspace(-3, 0, 300)   # A  (1 mA to 1 A)

# ERL luminosity: L ∝ N_e × f_rep = I/e
# ERL emittance is gun-limited (doesn't grow with current like a ring)
# Storage ring emittance: SR equilibrium (independent of current for Gaussian)
# But ring beam-beam limit restricts ξ < 0.05

eps_n_ring = 1.0e-6   # m·rad (SR equilibrium emittance – much smaller!)
# At 10 GeV, E^4 dependence gives tiny SR emittance ~ 0.1 mm·mrad (normalised)
# but beam-beam limits # of bunches / current, so let's show current scaling

# For ERL: N_e scales with I; emittance is fixed by gun
N_e_sw = I_sweep / e_charge / f_rep
L_erl_sw = np.array([
    luminosity(Ne, N_p, f_rep,
               sigma_IP(eps_n_e_sp, beta_star, gamma_e_peak_sp),
               sigma_IP(eps_n_p, beta_star, gamma_p),
               hourglass(sigma_z, beta_star)) * 1e-4
    for Ne in N_e_sw])

# For ring: same, but beam-beam limits ξ → once ξ > 0.05, lum saturates
gamma_ring = E_peak_sp / m_e_MeV
sig_e_ring = sigma_IP(eps_n_ring, beta_star, gamma_ring)
# BB limit: N_e < 4π γ σ² ξ_max / r_e
N_e_bb_lim = 4*np.pi * gamma_ring * sig_e_ring**2 * 0.05 / r_e
I_bb_lim   = N_e_bb_lim * e_charge * f_rep

L_ring_sw = np.where(N_e_sw <= N_e_bb_lim,
    luminosity(N_e_sw, N_p, f_rep, sig_e_ring,
               sigma_IP(eps_n_p, beta_star, gamma_p),
               hourglass(sigma_z, beta_star)) * 1e-4,
    luminosity(N_e_bb_lim, N_p, f_rep, sig_e_ring,
               sigma_IP(eps_n_p, beta_star, gamma_p),
               hourglass(sigma_z, beta_star)) * 1e-4)

ax.loglog(I_sweep * 1e3, L_erl_sw, color=BLUE,  lw=2.2, label="ERL (no BB limit)")
ax.loglog(I_sweep * 1e3, L_ring_sw, color=ORANGE, lw=2.2, ls="--",
          label=f"Storage ring (ξ_max=0.05, ε_n={eps_n_ring*1e6:.1f} mm·mrad)")
ax.axvline(I_arc * 1e3, color=GREY, lw=1.2, ls="--", alpha=0.7,
           label=f"ERL design point ({I_arc*1e3:.0f} mA)")
ax.scatter([I_arc*1e3], [bb_sp["L_pinch_cm2s"]], s=80, color=RED, zorder=5,
           label=f"ERL design L = {bb_sp['L_pinch_cm2s']:.1e}")
ax.set_xlabel("Average current (mA)", fontsize=10)
ax.set_ylabel("Luminosity (cm⁻²s⁻¹)", fontsize=10)
set_panel(ax, "(e) Luminosity vs Current: ERL vs Ring")
ax.legend(fontsize=8)

# ── (f) Summary table ─────────────────────────────────────────────────────
ax = axes6[1, 2]
ax.axis("off")
table_data = [
    ["Parameter", "1-pass ERL", "2-pass ERL", "Unit"],
    ["e⁻ energy",    f"{E_peak_sp/1e3:.3f}", f"{E_peak_tp/1e3:.3f}", "GeV"],
    ["p energy",     f"{E_p/1e6:.0f}",       f"{E_p/1e6:.0f}",       "GeV"],
    ["√s",           f"{bb_sp['E_cm_GeV']:.1f}", f"{bb_tp['E_cm_GeV']:.1f}", "GeV"],
    ["σ_e at IP",    f"{bb_sp['sigma_e_um']:.2f}", f"{bb_tp['sigma_e_um']:.2f}", "μm"],
    ["σ_p at IP",    f"{bb_sp['sigma_p_um']:.2f}", f"{bb_tp['sigma_p_um']:.2f}", "μm"],
    ["ξ_e",          f"{bb_sp['xi_e']:.4f}",  f"{bb_tp['xi_e']:.4f}",  "—"],
    ["ξ_p",          f"{bb_sp['xi_p']:.2e}",  f"{bb_tp['xi_p']:.2e}",  "—"],
    ["D_e",          f"{bb_sp['D_e']:.3f}",   f"{bb_tp['D_e']:.3f}",   "—"],
    ["H_hg",         f"{bb_sp['H_hg']:.4f}",  f"{bb_tp['H_hg']:.4f}",  "—"],
    ["H_pinch",      f"{bb_sp['H_pinch']:.3f}",f"{bb_tp['H_pinch']:.3f}","—"],
    ["L (geo)",      f"{bb_sp['L_geo_cm2s']:.2e}", f"{bb_tp['L_geo_cm2s']:.2e}", "cm⁻²s⁻¹"],
    ["L (pinch)",    f"{bb_sp['L_pinch_cm2s']:.2e}",f"{bb_tp['L_pinch_cm2s']:.2e}","cm⁻²s⁻¹"],
]
tab = ax.table(cellText=table_data[1:], colLabels=table_data[0],
               loc="center", cellLoc="center")
tab.auto_set_font_size(False)
tab.set_fontsize(8.5)
tab.scale(1.25, 1.45)
# Header row styling
for j in range(4):
    tab[0, j].set_facecolor("#1f4e79")
    tab[0, j].set_text_props(color="white", fontweight="bold")
for i in range(1, len(table_data)):
    bg = "#dce6f1" if i % 2 == 0 else "white"
    for j in range(4):
        tab[i, j].set_facecolor(bg)
ax.set_title("(f) Beam-Beam Summary", fontsize=11, fontweight="bold")

fig6.tight_layout()
fig6.savefig("ERL_fig6_beam_beam.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig6_beam_beam.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 7 – Two-Pass ERL Analysis
# ─────────────────────────────────────────────────────────────────────────────

fig7, axes7 = plt.subplots(2, 3, figsize=(16, 10))
fig7.suptitle("Two-Pass ERL Analysis (5 GeV per pass)", fontsize=13, fontweight="bold")

# ── (a) Energy profile along two-pass ERL ─────────────────────────────────
ax = axes7[0, 0]
# Stations: inj → linac1_out → arc1 → linac2_out → IP → arc_dec1 → linac_dec1
#           → arc_dec2 → linac_dec2 → dump
pos_tp_n = [0,    1,      2,    3,      3.4,   4,      5,    6,      7,    8]
E_tp_n   = [E_inj, E_inj, E_mid_tp, E_mid_tp,
             E_peak_tp, E_peak_tp,
             E_peak_tp, E_mid_tp, E_mid_tp, E_inj]

E_tp_fine = PchipInterpolator(pos_tp_n, E_tp_n)(np.linspace(0, 8, 2000))
pos_fine8 = np.linspace(0, 8, 2000)

ax.plot(pos_fine8, E_tp_fine / 1e3, color=RED, lw=2.5, label="2-pass ERL")
ax.plot(pos_fine,  E_fine   / 1e3,  color=BLUE, lw=1.8, ls="--", alpha=0.6,
        label="1-pass ERL")
ax.fill_between(pos_fine8, E_inj/1e3, E_tp_fine/1e3,
                where=(E_tp_fine > E_inj), alpha=0.12, color=RED)

# Mark key points
for txt, xp, ep, col in [
    ("7 MeV",       0.0, E_inj,       "black"),
    (f"{E_mid_tp/1e3:.1f} GeV", 2.0, E_mid_tp, ORANGE),
    (f"{E_peak_tp/1e3:.3f} GeV", 3.0, E_peak_tp, RED),
    ("IP",          3.4, E_peak_tp,    RED),
    (f"{E_mid_tp/1e3:.1f} GeV", 6.0, E_mid_tp, ORANGE),
    ("7 MeV",       8.0, E_inj,       "black"),
]:
    ax.scatter(xp, ep/1e3, s=60, color=col, zorder=5)
    ax.annotate(txt, (xp, ep/1e3), xytext=(0, 10), textcoords="offset points",
                ha="center", fontsize=7.5, color=col)

# Shade the two linac segments
for xs, xe, fc in [(1, 2, "#ffe8d0"), (2, 3, "#ffe0b0"),
                   (5, 6, "#ffd0d0"), (6, 7, "#ffc0c0")]:
    ax.axvspan(xs, xe, alpha=0.20, color=fc)

ax.text(1.5, 3.0, "Pass 1\n+5 GeV", ha="center", fontsize=8, color=ORANGE)
ax.text(2.5, 7.5, "Pass 2\n+5 GeV", ha="center", fontsize=8, color=RED)
ax.text(5.5, 7.5, "Dec 1\n−5 GeV", ha="center", fontsize=8, color=RED)
ax.text(6.5, 3.0, "Dec 2\n−5 GeV", ha="center", fontsize=8, color=ORANGE)

set_panel(ax, "(a) 2-Pass ERL Energy Profile",
          "Position along lattice →", "Beam Energy (GeV)")
ax.legend(fontsize=8.5)

# ── (b) Phase space: 4 key stages for 2-pass ──────────────────────────────
tp_stages = [("injection", "Injection\n7 MeV"),
             ("acc_pass_2", "After pass 2\n10.007 GeV"),
             ("post_collision", "Post collision"),
             ("dec_pass_2", "After dec. pass 2\n7 MeV")]
for idx, (key, lbl) in enumerate(tp_stages):
    if idx >= 3:
        break

ax_ps = axes7[0, 1]
for (key, lbl), col in zip(tp_stages[:3], [BLUE, GREEN, RED]):
    phi_k, dE_k, Eref_k = snap_tp[key]
    z_mm = phase_to_mm(phi_k)
    ax_ps.scatter(z_mm, dE_k * 1e3, s=0.2, c=col, alpha=0.3, rasterized=True,
                  label=lbl)
set_panel(ax_ps, "(b) Phase Space at 3 Stages (2-pass ERL)",
          "z (mm)", "ΔE (keV)")
ax_ps.set_xlim(-6, 6)
ax_ps.legend(fontsize=8.5, markerscale=15)

# ── (c) Energy-spread histograms comparing 1-pass vs 2-pass ───────────────
ax_hist = axes7[0, 2]
for (key_sp, key_tp, stage, col1, col2) in [
        ("dump", "dec_pass_2", "At dump", BLUE, RED)]:
    phi_sp_k, dE_sp_k, E_sp_k = snap_sp[key_sp]
    phi_tp_k, dE_tp_k, E_tp_k = snap_tp[key_tp]
    ax_hist.hist(dE_sp_k * 1e3, bins=120, density=True,
                 color=col1, alpha=0.6, label=f"1-pass: σ={np.std(dE_sp_k)*1e3:.1f} keV")
    ax_hist.hist(dE_tp_k * 1e3, bins=120, density=True,
                 color=col2, alpha=0.6, label=f"2-pass: σ={np.std(dE_tp_k)*1e3:.1f} keV")
set_panel(ax_hist, f"(c) Energy Spread at Dump: 1-pass vs 2-pass",
          "ΔE at dump (keV)", "Density")
ax_hist.legend(fontsize=9)

# ── (d) Cavity count and gradient comparison ───────────────────────────────
ax_cav = axes7[1, 0]
configs_cav = ["1-pass ERL\n(10 GeV)", "2-pass ERL\n(5 GeV/pass)"]
N_cav_vals  = [bgt_sp["N_cav"], bgt_tp["N_cav"]]
V_cav_vals  = [V_cav, V_cav]  # same design cavities
G_vals      = [G_acc, G_acc]

x_pos = [0, 1]
bars  = ax_cav.bar(x_pos, N_cav_vals, width=0.5, color=[BLUE, RED],
                   edgecolor="white", lw=1.5)
for bar, n in zip(bars, N_cav_vals):
    ax_cav.text(bar.get_x() + bar.get_width()/2, n + 10, f"{n}", ha="center",
                va="bottom", fontsize=11, fontweight="bold")
ax_cav.set_xticks(x_pos)
ax_cav.set_xticklabels(configs_cav, fontsize=10)
ax_cav.set_ylabel("Number of SRF cavities", fontsize=10)
set_panel(ax_cav, "(d) SRF Cavity Count")
ax_cav2 = ax_cav.twinx()
ax_cav2.plot(x_pos, V_cav_vals, "ko--", ms=8, label="Cavity voltage (MV)")
ax_cav2.set_ylabel("V_cavity (MV)", fontsize=10)
ax_cav2.set_ylim(0, 25)

# ── (e) Power budget side-by-side ─────────────────────────────────────────
ax_pwr = axes7[1, 1]
components = ["Beam\ndump", "Injector\nRF", "SRF\ncryoplant",
              "SR in\narcs", "Residual\nRF"]
vals_1p = [bgt_sp["P_dump_MW"]*1e3, bgt_sp["P_injector_MW"]*1e3,
           bgt_sp["P_cryo_MW"]*1e3, bgt_sp["P_SR_arc_MW"]*1e3,
           bgt_sp["P_RF_residual_MW"]*1e3]
vals_2p = [bgt_tp["P_dump_MW"]*1e3, bgt_tp["P_injector_MW"]*1e3,
           bgt_tp["P_cryo_MW"]*1e3, bgt_tp["P_SR_arc_MW"]*1e3,
           bgt_tp["P_RF_residual_MW"]*1e3]
x       = np.arange(len(components))
w       = 0.35
ax_pwr.bar(x - w/2, vals_1p, w, color=BLUE,  label=f"1-pass  ({bgt_sp['P_total_MW']*1e3:.0f} kW total)",
           edgecolor="white", lw=1)
ax_pwr.bar(x + w/2, vals_2p, w, color=RED,   label=f"2-pass  ({bgt_tp['P_total_MW']*1e3:.0f} kW total)",
           edgecolor="white", lw=1)
ax_pwr.set_xticks(x)
ax_pwr.set_xticklabels(components, fontsize=9)
set_panel(ax_pwr, "(e) Power Budget: 1-pass vs 2-pass",
          "", "Wall-plug power (kW)")
ax_pwr.legend(fontsize=9)

# ── (f) ERL topology diagram (2-pass) ─────────────────────────────────────
ax_topo = axes7[1, 2]
ax_topo.set_xlim(0, 10)
ax_topo.set_ylim(-3.5, 3.5)
ax_topo.axis("off")

def draw_arc(ax, cx, cy, rx, ry, t1, t2, col, lw=2.0, ls="-"):
    theta = np.linspace(np.radians(t1), np.radians(t2), 200)
    x = cx + rx * np.cos(theta)
    y = cy + ry * np.sin(theta)
    ax.plot(x, y, color=col, lw=lw, ls=ls)

# SRF Linac box (horizontal)
draw_box(ax_topo, 5.0, 0.0, 6.5, 0.75, "#e8f0ff", BLUE, label="SRF LINAC  (5 GeV/pass)", fs=8.5)

# Arcs for acceleration (pass 1 and 2)
draw_arc(ax_topo, 5.0, 0.0, 3.2, 1.5,  0, 180, GREEN, lw=2.0)  # pass 1
draw_arc(ax_topo, 5.0, 0.0, 3.2, 2.8,  0, 180, CYAN,  lw=2.0)  # pass 2

# Arcs for deceleration
draw_arc(ax_topo, 5.0, 0.0, 3.2, 1.5, 180, 360, ORANGE, lw=2.0, ls="--")  # dec 1
draw_arc(ax_topo, 5.0, 0.0, 3.2, 2.8, 180, 360, RED,    lw=2.0, ls="--")  # dec 2

# IP marker
draw_box(ax_topo, 8.2, 0.0, 0.9, 0.6, "#fff0d0", ORANGE, label="IP", fs=9)

# Injector and dump
draw_box(ax_topo, 1.0, 0.0, 1.0, 0.6, "#d0e8ff", BLUE,  label="INJ\n7 MeV", fs=7.5)
draw_box(ax_topo, 9.0, 2.9, 1.0, 0.6, "#ffd0d0", RED,   label="DUMP\n7 MeV", fs=7.5)

# Labels on arcs
ax_topo.text(1.65, 1.65, "Arc 1\n(+ΔΦ₁)", color=GREEN, fontsize=7.5, ha="center")
ax_topo.text(1.65, 2.95, "Arc 2\n(+ΔΦ₂)", color=CYAN,  fontsize=7.5, ha="center")
ax_topo.text(1.65,-1.65, "Dec arc 1\n(+π)", color=ORANGE, fontsize=7.5, ha="center")
ax_topo.text(1.65,-2.95, "Dec arc 2\n(+π)", color=RED,    fontsize=7.5, ha="center")
ax_topo.text(5.0,  3.25, "5 GeV / pass × 2 = 10 GeV total",
             color=BLUE, fontsize=8.5, ha="center", fontweight="bold")
ax_topo.set_title("(f) Two-Pass ERL Topology", fontsize=11, fontweight="bold")

fig7.tight_layout()
fig7.savefig("ERL_fig7_twopass.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig7_twopass.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8 – Comprehensive Summary Dashboard
# ─────────────────────────────────────────────────────────────────────────────

fig8 = plt.figure(figsize=(16, 11))
gs = gridspec.GridSpec(3, 4, figure=fig8, hspace=0.48, wspace=0.42)
fig8.suptitle("ERL Comprehensive Summary Dashboard", fontsize=14, fontweight="bold")

# ── 8.1  RF timing diagram ─────────────────────────────────────────────────
ax81 = fig8.add_subplot(gs[0, :2])
t_rf  = np.linspace(0, 5 / f_rep, 5000)
V_RF_wf = np.cos(omega_RF * t_rf)
ax81.plot(t_rf * 1e9, V_RF_wf, color=BLUE, lw=1.5, alpha=0.7, label="RF voltage")

# Mark bunch positions
t_bunch = np.arange(0, 5) / f_rep
for i, tb in enumerate(t_bunch):
    clr_b = GREEN if i % 2 == 0 else RED
    ax81.axvline(tb * 1e9, color=clr_b, lw=1.5, ls="--", alpha=0.75)
    lbl = "Acc." if i % 2 == 0 else "Dec."
    ax81.text(tb * 1e9, 0.95, lbl, ha="center", fontsize=7.5, color=clr_b,
              rotation=90, va="top")

ax81.set_xlabel("Time (ns)", fontsize=10)
ax81.set_ylabel("Norm. RF voltage", fontsize=10)
set_panel(ax81, f"RF Timing: f_RF = {f_RF/1e6:.2f} MHz,  "
                f"f_rep = {f_rep/1e6:.2f} MHz  (harmonic = {harmonic})",
          grid=False)
custom_lines = [Line2D([0], [0], color=GREEN, lw=1.5, ls="--"),
                Line2D([0], [0], color=RED,   lw=1.5, ls="--"),
                Line2D([0], [0], color=BLUE,  lw=1.5)]
ax81.legend(custom_lines, ["Accelerating bunch", "Decelerating bunch", "RF voltage"],
            fontsize=8.5, loc="lower right")

# ── 8.2  Bunch charge and current ─────────────────────────────────────────
ax82 = fig8.add_subplot(gs[0, 2])
ax82.axis("off")
summary_text = (
    f"{'Beam parameters':^30}\n"
    f"{'─'*30}\n"
    f"{'Bunch charge:':<22} {Q_bunch*1e9:.1f} nC\n"
    f"{'Electrons / bunch:':<22} {N_e:.3e}\n"
    f"{'Bunch length (σ_z):':<22} {sigma_z*1e3:.1f} mm\n"
    f"{'Inj. energy:':<22} {E_inj:.1f} MeV\n"
    f"{'RF frequency:':<22} {f_RF/1e6:.2f} MHz\n"
    f"{'Rep rate:':<22} {f_rep/1e6:.2f} MHz\n"
    f"{'I_arc:':<22} {I_arc*1e3:.1f} mA\n"
    f"{'I_linac:':<22} {I_linac*1e3:.1f} mA\n"
    f"{'─'*30}\n"
    f"{'1-pass peak energy:':<22} {E_peak_sp/1e3:.3f} GeV\n"
    f"{'2-pass peak energy:':<22} {E_peak_tp/1e3:.3f} GeV\n"
    f"{'η_recovery:':<22} {bgt_sp['eta_recovery']*100:.4f}%\n"
)
ax82.text(0.02, 0.98, summary_text, transform=ax82.transAxes,
          va="top", ha="left", fontsize=8.5,
          fontfamily="monospace",
          bbox=dict(fc="#f0f4ff", ec=BLUE, lw=1.5, boxstyle="round,pad=0.5"))

# ── 8.3  Power scorecard ──────────────────────────────────────────────────
ax83 = fig8.add_subplot(gs[0, 3])
ax83.axis("off")
power_text = (
    f"{'Power budget (kW)':^32}\n"
    f"{'─'*32}\n"
    f"{'Component':<18}{'1-pass':>6} {'2-pass':>7}\n"
    f"{'─'*32}\n"
    f"{'Beam dump':<18}{bgt_sp['P_dump_MW']*1e3:6.0f} {bgt_tp['P_dump_MW']*1e3:7.0f}\n"
    f"{'Injector RF':<18}{bgt_sp['P_injector_MW']*1e3:6.0f} {bgt_tp['P_injector_MW']*1e3:7.0f}\n"
    f"{'SRF cryoplant':<18}{bgt_sp['P_cryo_MW']*1e3:6.0f} {bgt_tp['P_cryo_MW']*1e3:7.0f}\n"
    f"{'SR in arcs':<18}{bgt_sp['P_SR_arc_MW']*1e3:6.0f} {bgt_tp['P_SR_arc_MW']*1e3:7.0f}\n"
    f"{'Residual RF':<18}{bgt_sp['P_RF_residual_MW']*1e3:6.1f} {bgt_tp['P_RF_residual_MW']*1e3:7.1f}\n"
    f"{'─'*32}\n"
    f"{'TOTAL ERL':<18}{bgt_sp['P_total_MW']*1e3:6.0f} {bgt_tp['P_total_MW']*1e3:7.0f}\n"
    f"{'─'*32}\n"
    f"{'Ring ρ=300m,60mA':<18}{bgt_ring_lrg['P_total_MW']*1e3:6.0f}\n"
    f"{'Ring ρ= 50m,60mA':<18}{bgt_ring_cmp['P_total_MW']*1e3:6.0f}\n"
    f"{'Ring ρ= 50m,2.5A':<18}{bgt_ring_hiL['P_total_MW']*1e3:6.0f}\n"
    f"{'w/o ERL (602 MW!)':<18}{bgt_sp['P_no_recovery_MW']*1e3:6.0f}\n"
    f"{'─'*32}\n"
    f"{'Saving vs no-ERL':<18}{bgt_sp['savings_factor']:6.0f}×\n"
    f"{'Saving vs hi-lum':<18}{bgt_ring_hiL['P_total_MW']/bgt_sp['P_total_MW']:6.0f}×\n"
)
ax83.text(0.02, 0.98, power_text, transform=ax83.transAxes,
          va="top", ha="left", fontsize=8.5, fontfamily="monospace",
          bbox=dict(fc="#fff4f0", ec=RED, lw=1.5, boxstyle="round,pad=0.5"))

# ── 8.4  Energy profile: 1-pass and 2-pass overlaid ──────────────────────
ax84 = fig8.add_subplot(gs[1, :2])
ax84.plot(pos_fine,  E_fine   / 1e3, color=BLUE, lw=2.2, label="1-pass ERL")
ax84.plot(pos_fine8, E_tp_fine/ 1e3, color=RED,  lw=2.2, ls="--", label="2-pass ERL")
ax84.axhline(E_inj/1e3, color=GREY, lw=1, ls=":", alpha=0.7, label=f"E_inj = {E_inj} MeV")
set_panel(ax84, "Energy Profile Comparison: 1-pass vs 2-pass",
          "Position along lattice", "Energy (GeV)")
ax84.legend(fontsize=9)

# ── 8.5  Phase-error sensitivity ─────────────────────────────────────────
ax85 = fig8.add_subplot(gs[1, 2:])
ax85.plot(phi_err_arr, eta_sp_arr * 100, color=BLUE, lw=2.2, label="1-pass")
ax85.plot(phi_err_arr, eta_tp_arr * 100, color=RED,  lw=2.2, ls="--", label="2-pass")
ax85.axhline(100, color=GREY, lw=0.8, ls=":", alpha=0.5)

# Shade the "acceptable" region (η > 99%)
mask_ok_sp = eta_sp_arr >= 0.99
if np.any(mask_ok_sp):
    phi_lo = phi_err_arr[mask_ok_sp][0]
    phi_hi = phi_err_arr[mask_ok_sp][-1]
    ax85.axvspan(phi_lo, phi_hi, alpha=0.10, color=BLUE)

set_panel(ax85, "Energy Recovery Efficiency vs. Phase Error",
          "Deceleration phase error Δφ (deg)", "Recovery efficiency (%)")
ax85.set_ylim(85, 101.5)
ax85.legend(fontsize=9)
ax85.text(0.02, 0.15,
          "Shaded region: η > 99%\n(phase control tolerance)",
          transform=ax85.transAxes, fontsize=8.5, color=BLUE,
          bbox=dict(fc="white", alpha=0.7, ec=BLUE, lw=1))

# ── 8.6  Luminosity potential ─────────────────────────────────────────────
ax86 = fig8.add_subplot(gs[2, :2])
bstar_fine = np.linspace(0.01, 0.25, 300)
L_sp_arr   = []
L_tp_arr   = []
for bs in bstar_fine:
    Hf   = hourglass(sigma_z, bs)
    s_e1 = sigma_IP(eps_n_e_sp, bs, E_peak_sp / m_e_MeV)
    s_e2 = sigma_IP(eps_n_e_tp, bs, E_peak_tp / m_e_MeV)
    s_p  = sigma_IP(eps_n_p,    bs, gamma_p)
    Hp1  = 1.0 + 0.5 * disruption(N_p, r_e, sigma_z,
                                    E_peak_sp/m_e_MeV, s_e1)**0.45
    Hp2  = 1.0 + 0.5 * disruption(N_p, r_e, sigma_z,
                                    E_peak_tp/m_e_MeV, s_e2)**0.45
    L_sp_arr.append(luminosity(N_e, N_p, f_rep, s_e1, s_p, Hf) * Hp1 * 1e-4)
    L_tp_arr.append(luminosity(N_e, N_p, f_rep, s_e2, s_p, Hf) * Hp2 * 1e-4)

ax86.semilogy(bstar_fine*100, L_sp_arr, color=BLUE, lw=2.2, label="1-pass (pinch incl.)")
ax86.semilogy(bstar_fine*100, L_tp_arr, color=RED,  lw=2.2, ls="--", label="2-pass (pinch incl.)")
ax86.axvline(beta_star*100, color=GREY, lw=1.5, ls="--", alpha=0.7,
             label=f"Design β* = {beta_star*100:.0f} cm")
ax86.scatter([beta_star*100], [bb_sp["L_pinch_cm2s"]], s=80, color=BLUE, zorder=5)
ax86.scatter([beta_star*100], [bb_tp["L_pinch_cm2s"]], s=80, color=RED,  zorder=5)
ax86.axhline(1e33, color=ORANGE, lw=1.2, ls=":", alpha=0.7, label="10³³ cm⁻²s⁻¹ reference")
set_panel(ax86, "Peak Luminosity vs. β* (with pinch enhancement)",
          "β* (cm)", "Luminosity (cm⁻²s⁻¹)")
ax86.legend(fontsize=9)

# ── 8.7  Cryo power breakdown ─────────────────────────────────────────────
ax87 = fig8.add_subplot(gs[2, 2:])
categories7 = ["No ERL (no recovery)\n602 MW!",
               "Ring ρ=300m, 60 mA",
               "Ring ρ= 50m, 60 mA",
               "Ring ρ= 50m, 2.5 A\n(high-lum scenario)",
               "2-pass ERL, 60 mA",
               "1-pass ERL, 60 mA"]
P_total7    = [bgt_sp["P_no_recovery_MW"] * 1e3,
               bgt_ring_lrg["P_total_MW"] * 1e3,
               bgt_ring_cmp["P_total_MW"] * 1e3,
               bgt_ring_hiL["P_total_MW"] * 1e3,
               bgt_tp["P_total_MW"] * 1e3,
               bgt_sp["P_total_MW"] * 1e3]
colors7     = [RED, GREY, "#888888", ORANGE, PURPLE, BLUE]
bars7 = ax87.barh(categories7, P_total7, color=colors7, edgecolor="white", lw=1.2,
                  height=0.65)
for bar, v in zip(bars7, P_total7):
    ax87.text(v * 1.08, bar.get_y() + bar.get_height()/2,
              f"{v:.0f} kW" if v < 1000 else f"{v/1000:.2f} MW",
              va="center", fontsize=8.5)
ax87.set_xscale("log")
ax87.set_xlim(1, 3e6)
set_panel(ax87, "Total Wall-Plug Power Comparison (log scale)",
          "Wall-plug power (kW)", grid=True)

# Annotate savings
for bar, v_this in [(bars7[-2], P_total7[-2]), (bars7[-1], P_total7[-1])]:
    s1 = bgt_sp["P_no_recovery_MW"] * 1e3 / v_this
    s2 = bgt_ring_hiL["P_total_MW"]  * 1e3 / v_this
    ax87.text(1.5, bar.get_y() + bar.get_height()/2,
              f"  ×{s1:.0f} vs no-ERL | ×{s2:.0f} vs hi-lum ring",
              va="center", fontsize=7.5, color="white", fontweight="bold")

fig8.savefig("ERL_fig8_dashboard.png", dpi=150, bbox_inches="tight")
print("Saved: ERL_fig8_dashboard.png")

print("\n" + "=" * 72)
print("  SIMULATION COMPLETE")
print("=" * 72)
print("\nOutput figures:")
for i, desc in enumerate([
    "fig1  – ERL schematic and energy profile",
    "fig2  – Longitudinal phase-space evolution (single-pass)",
    "fig3  – RF phase relationship and beam-loading cancellation",
    "fig4  – Energy recovery efficiency vs. phase error",
    "fig5  – Wall-plug power budget breakdown",
    "fig6  – Beam-beam analysis (e-p collisions)",
    "fig7  – Two-pass ERL analysis",
    "fig8  – Comprehensive summary dashboard",
], start=1):
    print(f"  ERL_{desc}")

print(f"""
Key results:
  Single-pass ERL  total wall-plug : {bgt_sp['P_total_MW']*1e3:.0f} kW
  Two-pass ERL     total wall-plug : {bgt_tp['P_total_MW']*1e3:.0f} kW
  Without energy recovery          : {bgt_sp['P_no_recovery_MW']*1e3:.0f} kW
  Savings (1-pass vs. no-recovery) : {bgt_sp['savings_factor']:.0f}×

  Storage-ring scenarios (10 GeV, β*=5 cm):
    Large ring  ρ=300 m, 60 mA    : {bgt_ring_lrg['P_total_MW']*1e3:.0f} kW
    Compact ring ρ= 50 m, 60 mA   : {bgt_ring_cmp['P_total_MW']*1e3:.0f} kW
    Compact ring ρ= 50 m, 2.5 A   : {bgt_ring_hiL['P_total_MW']:.1f} MW  ← high-lum design
  ERL advantage vs high-lum ring   : {bgt_ring_hiL['P_total_MW']/bgt_sp['P_total_MW']:.0f}×

  Luminosity (1-pass, with pinch)  : {bb_sp['L_pinch_cm2s']:.2e} cm⁻²s⁻¹
  Luminosity (2-pass, with pinch)  : {bb_tp['L_pinch_cm2s']:.2e} cm⁻²s⁻¹
  Centre-of-mass energy √s         : {bb_sp['E_cm_GeV']:.1f} GeV  (e × p)
""")
