"""
LWFA-fed Compton Source Parameter Calculator
=============================================
Based on: "LWFA-fed Compton Source Based on Dual-Energy Accumulator Ring
           with Cryocooled Nb3Sn Twin-Axis Cavities"

Calculates Compton photon energy, flux estimates, and identifies
key missing parameters needed for a full design study.
"""

import math

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
c      = 2.998e8          # speed of light [m/s]
h      = 6.626e-34        # Planck constant [J·s]
hbar   = h / (2 * math.pi)
eV     = 1.602e-19        # electron-volt [J]
me_eV  = 0.511e6          # electron rest mass [eV]
me_J   = me_eV * eV       # electron rest mass [J]
e_C    = 1.602e-19        # elementary charge [C]
sigma_T = 6.652e-29       # Thomson cross-section [m²]
nm     = 1e-9             # nanometre
fs     = 1e-15            # femtosecond
mJ     = 1e-3             # millijoule
mu_m   = 1e-6             # micrometre


# ─────────────────────────────────────────────────────────────────────────────
# SEPARATOR HELPERS
# ─────────────────────────────────────────────────────────────────────────────
W = 72  # table width

def section(title):
    print(f"\n{'═'*W}")
    print(f"  {title}")
    print(f"{'═'*W}")

def subsection(title):
    print(f"\n  {'─'*68}")
    print(f"  {title}")
    print(f"  {'─'*68}")

def row(label, value, unit, comment=""):
    label_col = f"  {label:<34}"
    value_col = f"{value:<18}"
    unit_col  = f"{unit:<12}"
    print(f"{label_col}{value_col}{unit_col}{comment}")

def comment_line(text):
    print(f"  ► {text}")

def missing_row(i, param, why_needed, estimate=""):
    print(f"  [{i:02d}]  {param:<35} {why_needed:<28} {estimate}")


# ═════════════════════════════════════════════════════════════════════════════
# 1.  KNOWN INPUT PARAMETERS (from the concept document)
# ═════════════════════════════════════════════════════════════════════════════
section("1.  KNOWN INPUT PARAMETERS  (from concept document)")

E_electron_MeV   = 3.0          # electron kinetic energy [MeV]
E_inject_MeV     = 0.2          # LWFA injection energy   [MeV]
lambda_laser_nm  = 600.0        # Compton laser wavelength [nm]
lambda_target_nm = 3.0          # target Compton photon wavelength [nm]
E_laser_pulse_mJ = 35.0         # LWFA laser pulse energy [mJ]
tau_laser_fs     = 35.0         # LWFA laser pulse duration [fs]
I_ring_A         = 1.0          # hypothesised ring current [A]

print()
row("Electron kinetic energy",  f"{E_electron_MeV:.1f}",      "MeV",   "Output of twin-axis Nb3Sn cavities")
row("LWFA injection energy",    f"{E_inject_MeV:.2f}",        "MeV",   "After laser-plasma acceleration")
row("Compton laser wavelength", f"{lambda_laser_nm:.0f}",     "nm",    "Stored in high-finesse optical cavity")
row("Target photon wavelength", f"{lambda_target_nm:.0f}",    "nm",    "Water-window X-ray goal")
row("LWFA laser pulse energy",  f"{E_laser_pulse_mJ:.0f}",    "mJ",    "Drive laser for plasma acceleration")
row("LWFA laser pulse duration",f"{tau_laser_fs:.0f}",        "fs",    "Ultra-short pulse")
row("Ring current (hypothesis)",f"{I_ring_A:.1f}",            "A",     "Key hypothesis — see Section 4")


# ═════════════════════════════════════════════════════════════════════════════
# 2.  RELATIVISTIC ELECTRON PARAMETERS
# ═════════════════════════════════════════════════════════════════════════════
section("2.  RELATIVISTIC ELECTRON PARAMETERS")

E_total_eV  = (E_electron_MeV + me_eV/1e6) * 1e6   # total energy [eV]
gamma       = E_total_eV / me_eV                     # Lorentz factor
beta        = math.sqrt(1 - 1/gamma**2)              # relativistic beta
p_MeVc      = gamma * beta * me_eV / 1e6            # momentum [MeV/c]
Brho_Tm     = p_MeVc * 1e6 * eV / (c * e_C)        # magnetic rigidity [T·m]
rho_min_m   = Brho_Tm / 1.5                         # min bend radius at B=1.5 T [m]

print()
row("Total electron energy",     f"{E_total_eV/1e6:.3f}",   "MeV",   "Kinetic + rest mass")
row("Lorentz factor γ",          f"{gamma:.3f}",             "—",     "γ = E_total / m_e c²")
row("Relativistic β",            f"{beta:.4f}",              "—",     "v/c")
row("Momentum p",                f"{p_MeVc:.3f}",            "MeV/c", "")
row("Magnetic rigidity Bρ",      f"{Brho_Tm:.4f}",           "T·m",   "Sets minimum ring bend radius")
row("Min. bend radius (B=1.5T)", f"{rho_min_m*100:.2f}",    "cm",    "Very compact ring feasible")


# ═════════════════════════════════════════════════════════════════════════════
# 3.  COMPTON PHOTON ENERGY
# ═════════════════════════════════════════════════════════════════════════════
section("3.  COMPTON PHOTON ENERGY  (head-on collision)")

subsection("3a. Forward calculation: electron energy → photon energy")

E_laser_photon_eV = (h * c) / (lambda_laser_nm * nm) / eV   # laser photon energy [eV]
# Head-on inverse Compton: E_gamma ≈ 4*gamma^2 * E_laser  (non-recoil limit)
E_gamma_eV_approx = 4 * gamma**2 * E_laser_photon_eV
lambda_gamma_nm   = (h * c) / (E_gamma_eV_approx * eV) / nm

# Full kinematics (exact, head-on, electron moves in +z, laser in -z):
# E_γ = E_laser * (1 + beta)^2 * gamma^2 / (1 + E_laser*(1+beta)*gamma / (m_e c^2))
numerator   = E_laser_photon_eV * (1 + beta)**2 * gamma**2
denominator = 1 + (E_laser_photon_eV * (1 + beta) * gamma) / me_eV
E_gamma_eV_exact = numerator / denominator
lambda_gamma_nm_exact = (h * c) / (E_gamma_eV_exact * eV) / nm

print()
row("Laser photon energy",          f"{E_laser_photon_eV:.3f}",        "eV",    "At 600 nm")
row("Approx. Compton photon energy",f"{E_gamma_eV_approx:.1f}",        "eV",    "4γ²·E_laser (no recoil)")
row("Approx. photon wavelength",    f"{lambda_gamma_nm:.2f}",          "nm",    "")
row("Exact Compton photon energy",  f"{E_gamma_eV_exact:.1f}",         "eV",    "Full kinematics")
row("Exact photon wavelength",      f"{lambda_gamma_nm_exact:.2f}",    "nm",    "Recoil correction small at 3 MeV")

print()
comment_line(f"Target wavelength from document: {lambda_target_nm:.0f} nm  →  "
             f"calculated: {lambda_gamma_nm_exact:.2f} nm")
comment_line(f"To hit exactly {lambda_target_nm} nm, need E_electron ≈ "
             f"{me_eV/1e6 * math.sqrt((h*c/(lambda_target_nm*nm)/eV) / E_laser_photon_eV) / 1e0:.2f} γ units.")

# Back-calculate required gamma for exactly 3 nm
E_gamma_target_eV = (h * c) / (lambda_target_nm * nm) / eV
gamma_needed = math.sqrt(E_gamma_target_eV / (4 * E_laser_photon_eV))
E_electron_needed_MeV = gamma_needed * me_eV / 1e6 - me_eV / 1e6

subsection("3b. Reverse: required electron energy for exactly 3 nm photons")
print()
row("Target photon energy",         f"{E_gamma_target_eV:.1f}",        "eV",    "At λ = 3 nm")
row("Required γ",                   f"{gamma_needed:.2f}",             "—",     "")
row("Required electron energy",     f"{E_electron_needed_MeV:.2f}",    "MeV",   "Kinetic energy needed")
row("Document specifies",           f"{E_electron_MeV:.1f}",           "MeV",   "Close — small gap, tunable")

print()
comment_line("CONCLUSION: ~3 nm photons are achievable; exact wavelength can be tuned")
comment_line("by adjusting electron energy within the cavity structure.")
comment_line(f"Water window definition: 2.3 – 4.4 nm (280 – 540 eV). ✓ Satisfied.")


# ═════════════════════════════════════════════════════════════════════════════
# 4.  RING PARAMETERS AND ELECTRON POPULATION
# ═════════════════════════════════════════════════════════════════════════════
section("4.  ACCUMULATOR RING PARAMETERS")

# Assume compact ring circumference (not stated — parametric)
C_ring_m  = 2.0          # assumed circumference [m]  ← NOT in document
f_rev     = beta * c / C_ring_m
N_e_ring  = I_ring_A / (e_C * f_rev)

print()
comment_line("Ring circumference NOT stated. Using assumed value C = 2 m.")
comment_line("This is the largest source of uncertainty in flux calculation.")
print()
row("Assumed ring circumference",  f"{C_ring_m:.1f}",              "m",     "⚠ ASSUMED — critical parameter")
row("Revolution frequency",        f"{f_rev/1e6:.1f}",             "MHz",   "f_rev = βc / C")
row("Electrons in ring at 1 A",    f"{N_e_ring:.2e}",              "e⁻",    "N = I / (e · f_rev)")
row("Ring current (document)",     f"{I_ring_A:.1f}",              "A",     "⚠ Hypothesis, not demonstrated")
row("Energy acceptance needed",    "~1%",                          "—",     "Typical for isochronous rings")

print()
comment_line("Ampere-scale current at 3 MeV is highly ambitious.")
comment_line("Comparable machines (e.g., MLS, BESSY VSR) achieve mA–100 mA.")
comment_line("LWFA injection emittance and energy spread will limit accumulation.")


# ═════════════════════════════════════════════════════════════════════════════
# 5.  LASER CAVITY — STORED PHOTON NUMBER
# ═════════════════════════════════════════════════════════════════════════════
section("5.  COMPTON LASER CAVITY (stored photon number)")

# Parametric: cavity finesse and input coupling
E_stored_mJ  = 1.0        # stored laser energy [mJ] — assumed
E_laser_photon_J = E_laser_photon_eV * eV

N_photons_stored = (E_stored_mJ * mJ) / E_laser_photon_J

# Laser cavity length ~ ring section (assume matched)
L_cavity_m   = 0.5        # assumed interaction region length [m]
f_laser_round = c / (2 * L_cavity_m)   # round-trip frequency

print()
comment_line("Laser cavity finesse (F) and input coupler not specified.")
comment_line(f"Assuming stored energy ~ {E_stored_mJ} mJ (moderate-finesse cavity).")
print()
row("Stored laser energy (assumed)",  f"{E_stored_mJ:.1f}",          "mJ",    "⚠ ASSUMED — depends on finesse F")
row("Stored photon number N_ph",      f"{N_photons_stored:.2e}",      "photons","Per cavity round-trip")
row("Laser photon energy",            f"{E_laser_photon_eV:.3f}",     "eV",    "At 600 nm")
row("Cavity round-trip freq. (est.)", f"{f_laser_round/1e6:.0f}",     "MHz",   "For L_cav ~ 0.5 m")

print()
comment_line("High-finesse enhancement cavities (e.g., F~10000) can store 100s of mJ.")
comment_line("This dramatically improves flux — key design parameter.")


# ═════════════════════════════════════════════════════════════════════════════
# 6.  COMPTON PHOTON FLUX ESTIMATE
# ═════════════════════════════════════════════════════════════════════════════
section("6.  COMPTON PHOTON FLUX ESTIMATE")

# Interaction point beam parameters — assumed
sigma_e_m   = 10e-6      # electron beam spot size (rms) [m] — assumed
sigma_ph_m  = 10e-6      # laser waist at IP [m] — assumed
A_IP        = math.pi * (sigma_e_m)**2    # effective area [m²]

# Luminosity per pass [m⁻²]
# For Gaussian beams: L = N_e * N_ph / (4*pi*sigma_x*sigma_y)
sigma_eff   = math.sqrt(sigma_e_m**2 + sigma_ph_m**2)
N_gamma_per_pass = N_e_ring * N_photons_stored * sigma_T / (4 * math.pi * sigma_eff**2)

# How often do they collide?
# Assume electrons meet laser photons once per revolution
f_collision = f_rev    # collisions per second

flux_total = N_gamma_per_pass * f_collision

# Bandwidth estimate (natural Compton bandwidth ≈ 1/gamma)
BW_percent  = (1 / gamma) * 100

print()
comment_line("Beam overlap geometry is critical and NOT specified in document.")
comment_line(f"Assuming σ_e = σ_ph = {sigma_e_m*1e6:.0f} μm at interaction point.")
print()
row("Electron beam spot σ_e",       f"{sigma_e_m*1e6:.0f}",          "μm",    "⚠ ASSUMED")
row("Laser beam waist σ_ph",        f"{sigma_ph_m*1e6:.0f}",         "μm",    "⚠ ASSUMED")
row("Effective overlap area",       f"{A_IP:.2e}",                   "m²",    "")
row("Photons per collision pass",   f"{N_gamma_per_pass:.2e}",        "γ/pass","N_e·N_ph·σ_T / A_eff")
row("Collision frequency",          f"{f_collision/1e6:.1f}",         "MHz",   "= f_rev (1 IP assumed)")
row("Total photon flux",            f"{flux_total:.2e}",              "ph/s",  "")
row("Natural bandwidth (1/γ)",      f"{BW_percent:.1f}",              "%",     "Monochromaticity ~17%")

print()
comment_line(f"Estimated flux:  ~{flux_total:.1e} photons/s")
comment_line("Comparable to: Thomson sources at synchrotrons (~10^12–10^14 ph/s)")
comment_line("State-of-art compact Compton sources: ~10^12 ph/s (e.g., CBETA, STAR)")
comment_line("This estimate strongly assumes 1 A ring current is achieved.")


# ═════════════════════════════════════════════════════════════════════════════
# 7.  SENSITIVITY ANALYSIS
# ═════════════════════════════════════════════════════════════════════════════
section("7.  FLUX SENSITIVITY TO KEY UNKNOWNS")

print()
print(f"  {'Parameter':<30} {'Variation':<20} {'Flux scaling':<20} {'Effect'}")
print(f"  {'─'*66}")

sensitivities = [
    ("Ring current I",          "×10 (10 A vs 1 A)",   "×10 (linear)",      "Dominant factor"),
    ("Ring current I",          "÷10 (100 mA)",        "÷10",               "More realistic baseline"),
    ("Stored laser energy",     "×10 (10 mJ)",         "×10 (linear)",      "High-finesse cavity"),
    ("Beam spot size σ",        "÷2 (5 μm)",           "×4 (quadratic)",    "Tight focusing critical"),
    ("Beam spot size σ",        "×2 (20 μm)",          "÷4",                "Poor overlap kills flux"),
    ("Ring circumference C",    "×2 (4 m)",            "÷2",                "Fewer revolutions/s"),
    ("Number of IPs",           "×4 (4 IPs)",          "×4 (linear)",       "Multiple cavities"),
]

for p, var, scale, effect in sensitivities:
    print(f"  {p:<30} {var:<20} {scale:<20} {effect}")


# ═════════════════════════════════════════════════════════════════════════════
# 8.  KEY MISSING PARAMETERS — NEXT STEPS
# ═════════════════════════════════════════════════════════════════════════════
section("8.  KEY MISSING PARAMETERS AND NEXT-STEP ESTIMATIONS REQUIRED")

print()
print(f"  {'#':<6} {'Parameter':<36} {'Why Needed':<29} {'First Estimate / Range'}")
print(f"  {'─'*70}")

missing = [
    ( 1, "Ring circumference C [m]",
          "Sets f_rev, N_e at given I",
          "Lattice design: ~1–5 m"),
    ( 2, "Ring lattice (bends, quads, sexts)",
          "Stability, emittance, tune",
          "FODO or DBA lattice study"),
    ( 3, "Ring momentum compaction α_c",
          "Isochronism, bunch length",
          "Need α_c < 10⁻³ for short bunches"),
    ( 4, "Beam emittance ε [nm·rad]",
          "IP spot size, brightness",
          "LWFA: ~0.1–1 mm·mrad (normalised)"),
    ( 5, "Beam energy spread δE/E [%]",
          "Chromatic effects, bandwidth",
          "LWFA: ~1–5%; cavity acceptance?"),
    ( 6, "Injection efficiency from LWFA",
          "How much charge accumulates",
          "~pC/pulse × rep-rate needed"),
    ( 7, "LWFA rep-rate [Hz or kHz]",
          "Top-up injection rate",
          "35 mJ laser → ~100 Hz typical"),
    ( 8, "Beam lifetime in ring τ [ms–s]",
          "Determines required top-up rate",
          "IBS + Compton recoil + gas"),
    ( 9, "Laser cavity finesse F",
          "Stored photon number N_ph",
          "F ~ 1000–100000; E_stored ∝ F"),
    (10, "Laser cavity mode-locking / sync",
          "Must be synchronised to ring",
          "Timing jitter < 100 fs required"),
    (11, "IP beta functions β*_x, β*_y [m]",
          "Determines beam waist at IP",
          "β* ~ σ²/ε ~ 0.01–0.1 m"),
    (12, "Nb3Sn cavity gradient [MV/m]",
          "Determines number of cells",
          "Nb3Sn at 4K: ~10–20 MV/m"),
    (13, "Cryocooler capacity [W at 4K]",
          "Feasibility of mobile cooling",
          "Need ~1–5 W at 4 K; ≈100 W plug"),
    (14, "Compton recoil energy loss/turn",
          "Limits ring current / lifetime",
          "ΔE ≈ E_γ per scatter event"),
    (15, "Intrabeam scattering (IBS) time",
          "Emittance growth → beam loss",
          "Especially critical at low energy"),
]

for idx, param, why, estimate in missing:
    print(f"  [{idx:02d}]  {param:<35} {why:<28} {estimate}")


# ═════════════════════════════════════════════════════════════════════════════
# 9.  SUMMARY
# ═════════════════════════════════════════════════════════════════════════════
section("9.  SUMMARY")

print()
print(f"  {'Quantity':<40} {'Value':<20} {'Confidence'}")
print(f"  {'─'*70}")

summary = [
    ("Compton photon energy (exact)",    f"{E_gamma_eV_exact:.0f} eV",       "✅ High — physics well-known"),
    ("Compton wavelength (exact)",       f"{lambda_gamma_nm_exact:.1f} nm",   "✅ High — in water window"),
    ("Electron γ at 3 MeV",             f"{gamma:.2f}",                      "✅ High — directly stated"),
    ("Electrons in ring at 1 A",        f"{N_e_ring:.1e}",                   "⚠  Medium — C assumed = 2 m"),
    ("Stored laser photons",            f"{N_photons_stored:.1e}",            "⚠  Low — finesse unknown"),
    ("Photon flux (central estimate)",  f"~{flux_total:.0e} ph/s",           "⚠  Low — 2 orders uncertain"),
    ("Flux range (pessimistic)",        f"~1e13 ph/s",                        "⚠  100 mA, poor overlap"),
    ("Flux range (optimistic)",         f"~1e17 ph/s",                        "⚠  1 A, tight focus, hi-F"),
    ("Water-window compatibility",      "Yes (2.3–4.4 nm)",                  "✅ High"),
    ("Compact ring feasibility",        f"Bρ = {Brho_Tm:.3f} T·m → C~2 m",  "✅ High"),
]

for q, v, c_ in summary:
    print(f"  {q:<40} {v:<20} {c_}")

print()
print(f"{'═'*W}")
print(f"  Next recommended action: Ring lattice design (item #1–3 above)")
print(f"  followed by LWFA injection simulation (items #6–8).")
print(f"{'═'*W}")
print()
