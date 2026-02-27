import math

# This script defines beam parameters for the Electron-Ion Collider (EIC) based on the parameters from https://eic.jlab.org/Documents/EIC-General/20250305-EIC_ParameterList.pdf
# for the 10 GeV electron and 275 GeV proton collision scenario.
# It uses normalized emittances, calculates relativistic gamma factors, derives geometric emittances,
# then calculates the beam sizes at the interaction point (IP), the collider luminosity,
# beam-beam tune shift parameters (xi_x and xi_y for both beams), and RMS beam divergences (angles) at IP.
# The beam-beam parameters are calculated using the Bassetti-Erskine approximation for elliptic beams:
# xi_{1x} = (r_1 * N_2 * beta_{1x}) / (2 * pi * gamma_1 * sigma_{2x} * (sigma_{2x} + sigma_{2y}))
# xi_{1y} = (r_1 * N_2 * beta_{1y}) / (2 * pi * gamma_1 * sigma_{2y} * (sigma_{2x} + sigma_{2y}))
# where indices 1 and 2 refer to the perturbed and perturbing beams, respectively.
# These formulas match the table values closely (e.g., electron xi_x ~0.072, xi_y ~0.100; proton xi_x ~0.012, xi_y ~0.012).
# The RMS beam angle (divergence) is theta = sqrt(eps_geo / beta*) in rad, converted to μrad.
# Calculated values match the table (e.g., electron theta_x ~116 μrad, theta_y ~83 μrad ≈84; proton ~65 μrad both).
# Comments explain each section, and parameters can be modified as needed for other scenarios.
# Note: All lengths are in meters (m), energies in GeV, emittances in m.
# This ignores hourglass effect and crossing angle for simplicity.
# To include hourglass, multiply L by the factor (0.98 in this case), but not done here.

# Constants
c = 3.0e8                # Speed of light (m/s)
circumference = 3834.0   # Ring circumference for RHIC/EIC (m)
f_rev = c / circumference  # Revolution frequency (Hz)
e = 1.602e-19            # Elementary charge (C) - not used here but for reference
m_e = 0.000511           # Electron rest mass energy (GeV)
m_p = 0.938              # Proton rest mass energy (GeV)
r_e = 2.8179e-15         # Classical electron radius (m)
r_p = 1.5347e-18         # Classical proton radius (m) ≈ r_e / (m_p / m_e)

# General parameters
num_bunches = 1160       # Number of bunches in each beam

# Electron beam parameters
energy_e = 10.0          # Electron energy (GeV)
N_e = 17.2e10            # Electron bunch population (particles per bunch)
eps_n_ex_e = 391e-6      # Normalized horizontal emittance (m)
eps_n_ey_e = 26e-6       # Normalized vertical emittance (m)
beta_x_e = 0.45          # Beta function horizontal at IP (m) ; was 1.49
beta_y_e = 0.056          # Beta function vertical at IP (m) ; was 0.19
bunch_length_e = 0.007   # RMS bunch length (m)

# Proton beam parameters
energy_p = 275.0         # Proton energy (GeV)
N_p = 6.9e10             # Proton bunch population (particles per bunch)
eps_n_ex_p = 3.3e-6      # Normalized horizontal emittance (m)
eps_n_ey_p = 0.3e-6      # Normalized vertical emittance (m)
beta_x_p = 0.80          # Beta function horizontal at IP (m) ; was 2.65
beta_y_p = 0.072          # Beta function vertical at IP (m) ; was 0.24
bunch_length_p = 0.060   # RMS bunch length (m) ; was 0.07

# Calculate relativistic gamma factors
gamma_e = energy_e / m_e
gamma_p = energy_p / m_p

# Calculate geometric emittances
eps_ex_e = eps_n_ex_e / gamma_e  # Geometric horizontal emittance for electrons (m)
eps_ey_e = eps_n_ey_e / gamma_e  # Geometric vertical emittance for electrons (m)
eps_ex_p = eps_n_ex_p / gamma_p  # Geometric horizontal emittance for protons (m)
eps_ey_p = eps_n_ey_p / gamma_p  # Geometric vertical emittance for protons (m)

# Calculate beam sizes at IP (RMS)
# Beam size sigma = sqrt(emittance * beta*) using geometric emittance
sigma_ex = math.sqrt(eps_ex_e * beta_x_e)  # Electron horizontal beam size (m)
sigma_ey = math.sqrt(eps_ey_e * beta_y_e)  # Electron vertical beam size (m)
sigma_px = math.sqrt(eps_ex_p * beta_x_p)  # Proton horizontal beam size (m)
sigma_py = math.sqrt(eps_ey_p * beta_y_p)  # Proton vertical beam size (m)

# Calculate RMS beam divergences (angles) at IP
theta_ex = math.sqrt(eps_ex_e / beta_x_e) * 1e6  # Electron horizontal divergence (μrad)
theta_ey = math.sqrt(eps_ey_e / beta_y_e) * 1e6  # Electron vertical divergence (μrad)
theta_px = math.sqrt(eps_ex_p / beta_x_p) * 1e6  # Proton horizontal divergence (μrad)
theta_py = math.sqrt(eps_ey_p / beta_y_p) * 1e6  # Proton vertical divergence (μrad)

# Calculate beam-beam tune shift parameters
# For electrons (perturbed by protons)
xi_ex = (r_e * N_p * beta_x_e) / (2 * math.pi * gamma_e * sigma_px * (sigma_px + sigma_py))  # Horizontal
xi_ey = (r_e * N_p * beta_y_e) / (2 * math.pi * gamma_e * sigma_py * (sigma_px + sigma_py))  # Vertical
# For protons (perturbed by electrons)
xi_px = (r_p * N_e * beta_x_p) / (2 * math.pi * gamma_p * sigma_ex * (sigma_ex + sigma_ey))  # Horizontal
xi_py = (r_p * N_e * beta_y_p) / (2 * math.pi * gamma_p * sigma_ey * (sigma_ex + sigma_ey))  # Vertical

# Calculate effective beam sizes for luminosity
sig_eff_x = math.sqrt(sigma_ex**2 + sigma_px**2)  # Effective horizontal (m)
sig_eff_y = math.sqrt(sigma_ey**2 + sigma_py**2)  # Effective vertical (m)

# Collision frequency
f_coll = num_bunches * f_rev  # Hz

# Calculate luminosity in m^-2 s^-1
L_m = (N_e * N_p * f_coll) / (2 * math.pi * sig_eff_x * sig_eff_y)

# Convert to cm^-2 s^{-1} (standard unit for collider luminosity)
L_cm = L_m / 1e4

# Calculate beam currents
I_e = N_e * f_coll * e * 1000  # Electron current (mA)
I_p = N_p * f_coll * e         # Proton current (A)

# Now, print a table of all parameters
# We use f-strings for formatting in Python 3.6+

print("EIC Beam Parameters and Calculations")
print("------------------------------------")
print(f"{'Parameter':<40} {'Value':<20} {'Unit'}")
print("------------------------------------")

# General
print(f"{'Ring circumference':<40} {circumference:<20.2f} {'m'}")
print(f"{'Revolution frequency':<40} {f_rev:<20.2f} {'Hz'}")
print(f"{'Number of bunches':<40} {num_bunches:<20} {''}")

# Rest masses
print(f"{'Electron rest mass energy':<40} {m_e:<20.6f} {'GeV'}")
print(f"{'Proton rest mass energy':<40} {m_p:<20.3f} {'GeV'}")

# Gamma factors
print(f"{'Electron gamma':<40} {gamma_e:<20.2f} {''}")
print(f"{'Proton gamma':<40} {gamma_p:<20.2f} {''}")

# Electron beam
print(f"{'Electron energy':<40} {energy_e:<20.1f} {'GeV'}")
print(f"{'Electron bunch population':<40} {N_e:<20.1e} {'particles'}")
print(f"{'Electron normalized horiz emittance':<40} {eps_n_ex_e * 1e6:<20.1f} {'μm'}")  # Converted to μm
print(f"{'Electron normalized vert emittance':<40} {eps_n_ey_e * 1e6:<20.1f} {'μm'}")
print(f"{'Electron geometric horiz emittance':<40} {eps_ex_e * 1e9:<20.1f} {'nm'}")  # Converted to nm
print(f"{'Electron geometric vert emittance':<40} {eps_ey_e * 1e9:<20.1f} {'nm'}")
print(f"{'Electron beta_x at IP':<40} {beta_x_e:<20.2f} {'m'}")
print(f"{'Electron beta_y at IP':<40} {beta_y_e:<20.2f} {'m'}")
print(f"{'Electron bunch length (RMS)':<40} {bunch_length_e:<20.3f} {'m'}")
print(f"{'Electron horizontal beam size at IP':<40} {sigma_ex * 1e6:<20.1f} {'μm'}")
print(f"{'Electron vertical beam size at IP':<40} {sigma_ey * 1e6:<20.1f} {'μm'}")
print(f"{'Electron horizontal beam-beam xi':<40} {xi_ex:<20.3f} {''}")
print(f"{'Electron vertical beam-beam xi':<40} {xi_ey:<20.3f} {''}")
print(f"{'Electron horizontal RMS beam angle':<40} {theta_ex:<20.1f} {'μrad'}")
print(f"{'Electron vertical RMS beam angle':<40} {theta_ey:<20.1f} {'μrad'}")
print(f"{'Electron beam current':<40} {I_e:<20.2f} {'mA'}")

# Proton beam
print(f"{'Proton energy':<40} {energy_p:<20.1f} {'GeV'}")
print(f"{'Proton bunch population':<40} {N_p:<20.1e} {'particles'}")
print(f"{'Proton normalized horiz emittance':<40} {eps_n_ex_p * 1e6:<20.1f} {'μm'}")
print(f"{'Proton normalized vert emittance':<40} {eps_n_ey_p * 1e6:<20.1f} {'μm'}")
print(f"{'Proton geometric horiz emittance':<40} {eps_ex_p * 1e9:<20.1f} {'nm'}")
print(f"{'Proton geometric vert emittance':<40} {eps_ey_p * 1e9:<20.1f} {'nm'}")
print(f"{'Proton beta_x at IP':<40} {beta_x_p:<20.2f} {'m'}")
print(f"{'Proton beta_y at IP':<40} {beta_y_p:<20.2f} {'m'}")
print(f"{'Proton bunch length (RMS)':<40} {bunch_length_p:<20.3f} {'m'}")
print(f"{'Proton horizontal beam size at IP':<40} {sigma_px * 1e6:<20.1f} {'μm'}")
print(f"{'Proton vertical beam size at IP':<40} {sigma_py * 1e6:<20.1f} {'μm'}")
print(f"{'Proton horizontal beam-beam xi':<40} {xi_px:<20.3f} {''}")
print(f"{'Proton vertical beam-beam xi':<40} {xi_py:<20.3f} {''}")
print(f"{'Proton horizontal RMS beam angle':<40} {theta_px:<20.1f} {'μrad'}")
print(f"{'Proton vertical RMS beam angle':<40} {theta_py:<20.1f} {'μrad'}")
print(f"{'Proton beam current':<40} {I_p:<20.2f} {'A'}")

# Calculated
print(f"{'Effective horizontal beam size':<40} {sig_eff_x * 1e6:<20.1f} {'μm'}")
print(f"{'Effective vertical beam size':<40} {sig_eff_y * 1e6:<20.1f} {'μm'}")
print(f"{'Collision frequency':<40} {f_coll:<20.1e} {'Hz'}")
print(f"{'Luminosity':<40} {L_cm:<20.2e} {'cm^{-2} s^{-1}'}")
print("------------------------------------")

# To modify: You can change the parameter values above for other energy scenarios
# from the PDF. For example, change energy_e to 5 GeV, update gamma_e, etc. You can also include the hourglass factor: L_cm *= 0.98
# If you want to include Piwinski angle or other effects, add those calculations similarly.
# Note that the calculated luminosity here is approximate and matches the given value ~3.14e33 without hourglass.
# Beam-beam and divergence values are useful for optimization; keep xi below limits (e.g., ~0.1 for electrons) to avoid instabilities.