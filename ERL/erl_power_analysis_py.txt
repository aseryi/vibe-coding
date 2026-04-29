#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-Pass Energy Recovery Linac (ERL) — Comprehensive Wall-Plug Power Analysis

Parameters:
  Bunch charge:   1.5 nC
  Bunch length:   1.5 mm (rms)
  Rep rate:       40.08 MHz  → 60 mA circulating current
  Injection:      7 MeV
  SRF gain:       10 GeV  (801.58 MHz, 4.5 K)
  Return arc:     180° phase shift → deceleration pass
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend; remove this line to get an interactive window
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LogLocator

# ── Physical constants ────────────────────────────────────────────────────────
C_GAMMA = 8.85e-5   # m/GeV³  (synchrotron radiation constant for electrons)
# Exact meaning: U0 [GeV] = C_GAMMA × E[GeV]^4 / ρ[m]  for a full isomagnetic turn
# For a partial arc of total angle θ [rad]:
#   ΔE [GeV] = C_GAMMA × E^4 / ρ × (θ / 2π)

# ── Machine parameters ────────────────────────────────────────────────────────
BUNCH_CHARGE   = 1.5e-9     # C
BUNCH_LENGTH   = 1.5e-3     # m  rms
REP_RATE       = 40.08e6    # Hz
I_BEAM         = 60e-3      # A  (= BUNCH_CHARGE × REP_RATE)
I_LINAC        = 120e-3     # A  (acc + dec beams coexist in linac)

E_INJ          = 7e-3       # GeV  (7 MeV)
E_GAIN         = 10.0       # GeV
E_PEAK         = E_INJ + E_GAIN  # GeV

F_RF           = 801.58e6   # Hz  SRF frequency
T_SRF          = 4.5        # K   operating temperature
T_ROOM         = 300.0      # K
Q0             = 2e10       # SRF intrinsic quality factor
R_OVER_Q       = 100.0      # Ω   per cavity (circuit convention)
V_CAVITY       = 25e6       # V   per cavity → 400 cavities = 10 GV
N_CAVITIES     = int(round(E_GAIN * 1e9 / V_CAVITY))   # 400

RHO_ARC        = 100.0      # m   bending radius in return arc
THETA_ARC      = np.pi      # rad 180° total bending

CRYO_FRAC      = 0.25       # fraction of Carnot efficiency achieved
ETA_KLYSTRON   = 0.65
ETA_XMISSION   = 0.90
ETA_RF         = ETA_KLYSTRON * ETA_XMISSION   # 0.585

SR_ENERGY      = 10.0       # GeV  reference storage ring
SR_RHO         = 80.0       # m    average bending radius
SR_CURRENT     = 200e-3     # A
ETA_RF_SR      = 0.55       # RF efficiency for normal-conducting SR cavities

# ── Derived quantities ────────────────────────────────────────────────────────
I_CALC = BUNCH_CHARGE * REP_RATE   # should be ≈ 60 mA

# SR energy loss in return arc (single pass, 180°)
U_ARC   = C_GAMMA * E_PEAK**4 / RHO_ARC * (THETA_ARC / (2 * np.pi))   # GeV

E_BEFORE_DECEL  = E_PEAK - U_ARC           # GeV  (energy entering decel section)
E_RECOVERED     = E_BEFORE_DECEL - E_INJ   # GeV  (returned to RF cavities)
E_NOT_RECOVERED = E_GAIN - E_RECOVERED     # GeV  (≈ U_ARC, must be supplied by RF)
ETA_RECOVERY    = E_RECOVERED / E_GAIN

# Power flow  [P in W = E_eV × I_A,  E_eV = E_GeV × 1e9]
P_INJ_BEAM      = E_INJ   * 1e9 * I_BEAM   # W  injection beam power
P_PEAK_BEAM     = E_PEAK  * 1e9 * I_BEAM   # W  beam power at full energy
P_RECOVERED_W   = E_RECOVERED * 1e9 * I_BEAM
P_SR_ARC        = U_ARC   * 1e9 * I_BEAM   # W  SR power lost in arc
P_RF_BEAM_NET   = E_NOT_RECOVERED * 1e9 * I_BEAM   # W  ≈ P_SR_ARC
P_DUMP          = E_INJ   * 1e9 * I_BEAM   # W  ≈ injection power
P_WITHOUT_ER    = E_GAIN  * 1e9 * I_BEAM   # W  beam power WITHOUT recovery (600 MW!)

# SRF cavity ohmic (wall) losses
P_WALL_CAV  = V_CAVITY**2 / (2 * R_OVER_Q * Q0)   # W per cavity
P_WALL_TOT  = P_WALL_CAV * N_CAVITIES              # W total

# RF generators must supply: SR compensation + wall losses
P_RF_NEEDED = P_RF_BEAM_NET + P_WALL_TOT
P_RF_WP     = P_RF_NEEDED / ETA_RF   # W  wall-plug

# Cryogenic system
COP_CARNOT    = T_SRF / (T_ROOM - T_SRF)
COP_ACTUAL    = COP_CARNOT * CRYO_FRAC
P_STATIC_CRYO = 0.3 * N_CAVITIES              # W  static heat load (~0.3 W/cavity)
P_CRYO_LOAD   = P_WALL_TOT + P_STATIC_CRYO   # W  total heat to remove at 4.5 K
P_CRYO_WP     = P_CRYO_LOAD / COP_ACTUAL     # W  cryogenic wall-plug

# Supporting systems
P_INJECTOR = P_INJ_BEAM / 0.10   # W  injector (10% wall-plug → beam efficiency)
P_MAGNETS  = 2.0e6               # W  arc magnets estimate
P_INFRA    = 1.0e6               # W  cooling water, controls, vacuum

P_ERL_TOTAL = P_RF_WP + P_CRYO_WP + P_INJECTOR + P_MAGNETS + P_INFRA

# Storage ring
U0_SR      = C_GAMMA * SR_ENERGY**4 / SR_RHO   # GeV per turn
P_SR_BEAM  = U0_SR * 1e9 * SR_CURRENT          # W  SR power to be compensated
P_RF_SR_WP = P_SR_BEAM / ETA_RF_SR             # W  RF wall-plug
P_MAG_SR   = 15.0e6                            # W  larger magnet system
P_INF_SR   = 5.0e6                             # W
P_SR_TOTAL = P_RF_SR_WP + P_MAG_SR + P_INF_SR

# ── Print report ──────────────────────────────────────────────────────────────
SEP = "═" * 68

def mw(w): return w / 1e6

print(SEP)
print("   SINGLE-PASS ERL — COMPREHENSIVE WALL-PLUG POWER ANALYSIS")
print(SEP)

print(f"""
BEAM PARAMETERS
  Bunch charge               {BUNCH_CHARGE*1e9:.2f} nC
  Bunch length (rms)         {BUNCH_LENGTH*1e3:.1f} mm
  Repetition rate            {REP_RATE/1e6:.3f} MHz
  Beam current (calculated)  {I_CALC*1e3:.2f} mA
  Circulating current        {I_BEAM*1e3:.1f} mA
  Within linac               {I_LINAC*1e3:.1f} mA  (acc + dec beams)
  SRF frequency              {F_RF/1e6:.2f} MHz
  SRF temperature            {T_SRF} K
  Harmonic ratio             {F_RF/REP_RATE:.1f}  (RF / rep-rate)
  Cavities                   {N_CAVITIES}  ×  {V_CAVITY/1e6:.0f} MV  =  {E_GAIN:.0f} GeV total gain

ENERGY PROFILE
  Injection energy           {E_INJ*1e3:.1f} MeV
  RF gain                    {E_GAIN:.3f} GeV
  Peak energy                {E_PEAK*1e3:.1f} MeV  ({E_PEAK:.4f} GeV)
  SR loss in 180° return arc {U_ARC*1e3:.3f} MeV   (ρ = {RHO_ARC:.0f} m)
  Energy before decel        {E_BEFORE_DECEL*1e3:.3f} MeV
  Energy recovered to RF     {E_RECOVERED*1e3:.4f} MeV
  Energy NOT recovered       {E_NOT_RECOVERED*1e3:.4f} MeV   (= SR loss; must be supplied)
  Beam dump energy           ≈ {E_INJ*1e3:.1f} MeV  (≈ injection energy)

ENERGY RECOVERY
  η_recovery                 {ETA_RECOVERY*100:.4f} %
  Peak beam power            {mw(P_PEAK_BEAM):.2f} MW
  Power returned to RF       {mw(P_RECOVERED_W):.3f} MW
  WITHOUT energy recovery    {mw(P_WITHOUT_ER):.1f} MW  ← beam power alone!
  Net RF needed for beam     {P_RF_BEAM_NET/1e3:.2f} kW   (compensates only SR losses)
""")

print("RF SYSTEM")
print(f"  SRF wall loss per cavity   {P_WALL_CAV:.1f} W")
print(f"  Total SRF wall losses      {P_WALL_TOT/1e3:.1f} kW")
print(f"  Net RF for SR compensation {P_RF_BEAM_NET/1e3:.2f} kW")
print(f"  Total RF needed            {P_RF_NEEDED/1e3:.1f} kW")
print(f"  Klystron efficiency        {ETA_KLYSTRON*100:.0f}%   transmission {ETA_XMISSION*100:.0f}%")
print(f"  RF wall-plug               {mw(P_RF_WP):.3f} MW")

print(f"""
CRYOGENIC SYSTEM  (4.5 K SRF)
  Carnot COP                 {COP_CARNOT*100:.3f} %
  Actual COP ({CRYO_FRAC*100:.0f}% of Carnot)  {COP_ACTUAL*100:.4f} %
  Dynamic heat load          {P_WALL_TOT:.1f} W
  Static heat load           {P_STATIC_CRYO:.1f} W
  Total cryo load at 4.5 K   {P_CRYO_LOAD:.1f} W  ({P_CRYO_LOAD/1e3:.2f} kW)
  Cryo wall-plug             {mw(P_CRYO_WP):.3f} MW
  Specific power             {P_CRYO_WP/P_CRYO_LOAD:.0f} W_wall / W_cryo
""")

print("ERL TOTAL WALL-PLUG BUDGET")
erl_items = [
    ("RF generators",    P_RF_WP),
    ("Cryogenic system", P_CRYO_WP),
    ("Injector",         P_INJECTOR),
    ("Arc magnets",      P_MAGNETS),
    ("Infrastructure",   P_INFRA),
]
for name, pwr in erl_items:
    print(f"  {name:<22}  {mw(pwr):7.3f} MW   ({pwr/P_ERL_TOTAL*100:5.1f}%)")
print(f"  {'TOTAL':<22}  {mw(P_ERL_TOTAL):7.3f} MW")

print(f"""
10 GeV STORAGE RING REFERENCE  (200 mA, ρ = {SR_RHO:.0f} m)
  SR energy loss per turn    {U0_SR*1e3:.2f} MeV
  SR beam power              {mw(P_SR_BEAM):.3f} MW""")
sr_items = [
    ("RF for SR losses",  P_RF_SR_WP),
    ("Magnets",           P_MAG_SR),
    ("Infrastructure",    P_INF_SR),
]
for name, pwr in sr_items:
    print(f"  {name:<22}  {mw(pwr):7.3f} MW   ({pwr/P_SR_TOTAL*100:5.1f}%)")
print(f"  {'TOTAL':<22}  {mw(P_SR_TOTAL):7.3f} MW")

print(f"""
SUMMARY
  ERL total wall-plug        {mw(P_ERL_TOTAL):.2f} MW
  Storage ring wall-plug     {mw(P_SR_TOTAL):.2f} MW
  Power ratio  (ERL / SR)    {P_ERL_TOTAL/P_SR_TOTAL:.2f}×
  Without ERL beam power     {mw(P_WITHOUT_ER/ETA_RF):.0f} MW  (RF alone, no recovery)
""")
print(SEP)

# ── Visualisation ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 15))
fig.suptitle(
    "Single-Pass Energy Recovery Linac — Comprehensive Wall-Plug Power Analysis",
    fontsize=15, fontweight="bold", y=0.99,
)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38)

BLUE   = "#1565C0"
RED    = "#C62828"
GREEN  = "#2E7D32"
ORANGE = "#E65100"
PURPLE = "#6A1B9A"
TEAL   = "#00695C"
GREY   = "#546E7A"
BROWN  = "#4E342E"

# ── 1. Energy profile along ERL ───────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])

# Seven representative stations
x_pos = [0, 0.8, 1.6, 2.2, 2.8, 3.6, 4.4]
e_mev = [
    E_INJ * 1e3,                                           # injector
    E_INJ * 1e3 + E_GAIN * 1e3 * 0.5,                     # linac mid (acc)
    E_PEAK * 1e3,                                          # end of acc
    (E_PEAK - U_ARC * 0.5) * 1e3,                         # arc midpoint
    E_BEFORE_DECEL * 1e3,                                  # end of arc
    E_BEFORE_DECEL * 1e3 - E_GAIN * 1e3 * 0.5,            # linac mid (dec)
    E_INJ * 1e3,                                           # dump
]
xlabels = ["Injector\n7 MeV", "Linac\n(acc) ½", "End of\nlinac\n10007 MeV",
           "Arc\nmidpoint", "Post-arc\n(pre-decel)", "Linac\n(dec) ½", "Beam\ndump"]

ax1.plot(x_pos, e_mev, "o-", color=BLUE, lw=2.5, ms=7, zorder=5)
ax1.fill_between(x_pos, 0, e_mev, alpha=0.10, color=BLUE)
ax1.axhline(E_INJ * 1e3, color=GREY, ls="--", lw=1.3, alpha=0.7,
            label=f"Injection energy  {E_INJ*1e3:.0f} MeV")
ax1.annotate(
    f"SR loss\n{U_ARC*1e3:.2f} MeV\n(ρ = {RHO_ARC:.0f} m)",
    xy=(2.8, E_BEFORE_DECEL * 1e3),
    xytext=(3.0, E_PEAK * 1e3 * 0.70),
    arrowprops=dict(arrowstyle="->", color=RED, lw=1.5),
    color=RED, fontsize=9, va="top",
)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(xlabels, fontsize=8.5)
ax1.set_ylabel("Beam energy (MeV)", fontsize=11)
ax1.set_title("ERL Energy Profile Along Beamline", fontsize=12, fontweight="bold")
ax1.set_ylim(-200, E_PEAK * 1e3 * 1.18)
ax1.legend(fontsize=9, loc="upper right")
ax1.grid(True, alpha=0.3)

for span, color, label in [
    ((0, 1.6), "#1565C0", "Acceleration"),
    ((1.6, 2.8), "#E65100", "Return arc (180°)"),
    ((2.8, 4.4), "#2E7D32", "Deceleration"),
]:
    ax1.axvspan(*span, alpha=0.06, color=color)
    ax1.text(
        (span[0] + span[1]) / 2, E_PEAK * 1e3 * 1.10,
        label, ha="center", fontsize=8, color=color, style="italic",
    )

# ── 2. Power-flow bar chart ───────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
flow_names = ["Injection\nbeam", "Net RF\nsupplied", "Peak beam\npower",
              "Recovered\nto RF", "SR loss\n(arc)", "Dump\npower"]
flow_vals  = [P_INJ_BEAM, P_RF_NEEDED, P_PEAK_BEAM,
              P_RECOVERED_W, P_SR_ARC, P_DUMP]
flow_colors = [BLUE, RED, "#0D47A1", GREEN, ORANGE, TEAL]
bars = ax2.bar(range(len(flow_names)), [v / 1e6 for v in flow_vals],
               color=flow_colors, edgecolor="black", lw=0.6, width=0.62)
ax2.set_xticks(range(len(flow_names)))
ax2.set_xticklabels(flow_names, fontsize=7.8)
ax2.set_ylabel("Power (MW)", fontsize=10)
ax2.set_title("ERL Power Flow", fontsize=11, fontweight="bold")
ax2.grid(True, axis="y", alpha=0.3)
for bar, val in zip(bars, flow_vals):
    ax2.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + max(flow_vals) * 0.01,
             f"{val/1e6:.2f}", ha="center", va="bottom", fontsize=8)

# ── 3. ERL wall-plug pie ──────────────────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
erl_labels = ["RF generators", "Cryogenics", "Injector", "Arc magnets", "Infrastructure"]
erl_vals   = [P_RF_WP, P_CRYO_WP, P_INJECTOR, P_MAGNETS, P_INFRA]
erl_colors = [RED, PURPLE, BLUE, GREY, BROWN]
wedges, _, autotexts = ax3.pie(
    erl_vals, autopct="%1.1f%%", colors=erl_colors,
    startangle=90, pctdistance=0.70,
)
ax3.legend(wedges, erl_labels, loc="upper left", bbox_to_anchor=(-0.12, -0.02),
           fontsize=7.5, ncol=1)
ax3.set_title(f"ERL Wall-Plug Breakdown\n(Total: {mw(P_ERL_TOTAL):.2f} MW)",
              fontsize=10, fontweight="bold")

# ── 4. Storage-ring wall-plug pie ─────────────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
sr_labels = ["RF for SR losses", "Magnets", "Infrastructure"]
sr_vals   = [P_RF_SR_WP, P_MAG_SR, P_INF_SR]
sr_colors = [RED, GREY, BROWN]
wedges2, _, _ = ax4.pie(
    sr_vals, autopct="%1.1f%%", colors=sr_colors,
    startangle=90, pctdistance=0.70,
)
ax4.legend(wedges2, sr_labels, loc="upper left", bbox_to_anchor=(-0.12, -0.02),
           fontsize=7.5, ncol=1)
ax4.set_title(f"Storage Ring Wall-Plug Breakdown\n(Total: {mw(P_SR_TOTAL):.2f} MW)",
              fontsize=10, fontweight="bold")

# ── 5. Three-way comparison bar ───────────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2])
comp_labels = [f"ERL\n({I_BEAM*1e3:.0f} mA)", f"Storage ring\n({SR_CURRENT*1e3:.0f} mA)",
               "Linac w/o\nrecovery"]
comp_vals   = [mw(P_ERL_TOTAL), mw(P_SR_TOTAL), mw(P_WITHOUT_ER / ETA_RF)]
comp_colors = [BLUE, RED, ORANGE]
bars5 = ax5.bar(comp_labels, comp_vals, color=comp_colors,
                edgecolor="black", lw=1, width=0.45)
for bar, val in zip(bars5, comp_vals):
    ax5.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 4,
             f"{val:.1f} MW", ha="center", va="bottom",
             fontsize=10, fontweight="bold")
ax5.set_ylabel("Total wall-plug power (MW)", fontsize=10)
ax5.set_title("Architecture Comparison\n(same 10 GeV beam energy)",
              fontsize=11, fontweight="bold")
ax5.set_ylim(0, max(comp_vals) * 1.22)
ax5.grid(True, axis="y", alpha=0.3)
ax5.text(0.5, 0.91,
         f"ERL saves {comp_vals[2] - comp_vals[0]:.0f} MW vs no-recovery",
         transform=ax5.transAxes, ha="center", fontsize=8.5,
         color=GREEN, fontweight="bold")

# ── 6. Cryogenic power vs temperature ────────────────────────────────────────
ax6 = fig.add_subplot(gs[2, 0])
T_arr   = np.linspace(1.9, 10, 400)
COP_c   = T_arr / (T_ROOM - T_arr)
COP_r   = COP_c * CRYO_FRAC
P_cryo_arr = P_CRYO_LOAD / COP_r / 1e6   # MW

ax6.semilogy(T_arr, P_cryo_arr, color=PURPLE, lw=2.5)
ax6.axvline(T_SRF, color=RED, ls="--", lw=1.8,
            label=f"Design: {T_SRF} K → {mw(P_CRYO_WP):.1f} MW")
ax6.scatter([T_SRF], [mw(P_CRYO_WP)], color=RED, s=80, zorder=6)
ax6.set_xlabel("SRF operating temperature (K)", fontsize=10)
ax6.set_ylabel("Cryogenic wall-plug (MW)", fontsize=10)
ax6.set_title("Cryogenic Power vs Temperature\n"
              f"(heat load = {P_CRYO_LOAD:.0f} W, {CRYO_FRAC*100:.0f}% Carnot)",
              fontsize=10, fontweight="bold")
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3, which="both")

# Annotate 2 K point (possible upgrade)
T2 = 2.0
COP2 = (T2 / (T_ROOM - T2)) * CRYO_FRAC
P2 = mw(P_CRYO_LOAD / COP2)
ax6.scatter([T2], [P2], color=BLUE, s=60, zorder=6)
ax6.annotate(f"2 K → {P2:.0f} MW", xy=(T2, P2),
             xytext=(3.5, P2 * 1.5), fontsize=8, color=BLUE,
             arrowprops=dict(arrowstyle="->", color=BLUE))

# ── 7. Recovery efficiency vs arc bending radius ──────────────────────────────
ax7  = fig.add_subplot(gs[2, 1])
rho_arr = np.linspace(5, 600, 500)
U_arr   = C_GAMMA * E_PEAK**4 / rho_arr * (np.pi / (2 * np.pi))   # GeV, 180° arc
eta_arr = (1 - U_arr / E_GAIN) * 100   # %
Pnet_arr = U_arr * 1e9 * I_BEAM / 1e3  # kW  net RF needed

ax7b = ax7.twinx()
l1, = ax7.plot(rho_arr, eta_arr, color=BLUE, lw=2.5, label="Recovery efficiency")
l2, = ax7b.plot(rho_arr, Pnet_arr, color=RED, lw=2, ls="--", label="Net RF power (kW)")
ax7.axvline(RHO_ARC, color=GREEN, ls=":", lw=2,
            label=f"Design ρ = {RHO_ARC:.0f} m")
ax7.set_xlabel("Arc bending radius ρ (m)", fontsize=10)
ax7.set_ylabel("Energy recovery efficiency (%)", fontsize=10, color=BLUE)
ax7b.set_ylabel("Net RF power needed (kW)", fontsize=10, color=RED)
ax7.set_ylim(80, 100.5)
ax7.set_title("Recovery Efficiency & Net RF vs Arc Radius",
              fontsize=10, fontweight="bold")
lines = [l1, l2]
ax7.legend(lines, [l.get_label() for l in lines] + [f"Design ρ = {RHO_ARC:.0f} m"],
           fontsize=8, loc="lower right")
ax7.grid(True, alpha=0.3)

# Design point marker
ax7.scatter([RHO_ARC], [ETA_RECOVERY * 100], color=GREEN, s=80, zorder=6)
ax7b.scatter([RHO_ARC], [P_RF_BEAM_NET / 1e3], color=RED, s=80, zorder=6)

# ── 8. Summary table ──────────────────────────────────────────────────────────
ax8 = fig.add_subplot(gs[2, 2])
ax8.axis("off")
table_data = [
    ["Parameter", "ERL", "Storage Ring"],
    ["Beam energy",         "10 GeV",                          "10 GeV"],
    ["Beam current",        f"{I_BEAM*1e3:.0f} mA",            f"{SR_CURRENT*1e3:.0f} mA"],
    ["Energy recovery η",   f"{ETA_RECOVERY*100:.3f} %",       "N/A"],
    ["SR loss / pass",      f"{U_ARC*1e3:.2f} MeV*",          f"{U0_SR*1e3:.1f} MeV/turn"],
    ["SR power",            f"{mw(P_SR_ARC):.2f} MW",          f"{mw(P_SR_BEAM):.2f} MW"],
    ["RF wall-plug",        f"{mw(P_RF_WP):.3f} MW",           f"{mw(P_RF_SR_WP):.3f} MW"],
    ["Cryo wall-plug",      f"{mw(P_CRYO_WP):.2f} MW",        "—"],
    ["Injector wall-plug",  f"{mw(P_INJECTOR):.2f} MW",        "—"],
    ["Magnets",             f"{mw(P_MAGNETS):.2f} MW",          f"{mw(P_MAG_SR):.2f} MW"],
    ["TOTAL",               f"{mw(P_ERL_TOTAL):.2f} MW",       f"{mw(P_SR_TOTAL):.2f} MW"],
]

tbl = ax8.table(cellText=table_data, cellLoc="center", loc="center",
                bbox=[0, 0.04, 1.0, 0.96])
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.5)

# Header row
for j in range(3):
    tbl[0, j].set_facecolor(BLUE)
    tbl[0, j].set_text_props(color="white", fontweight="bold")
# Last (total) row
for j in range(3):
    tbl[len(table_data) - 1, j].set_facecolor("#C8E6C9")
    tbl[len(table_data) - 1, j].set_text_props(fontweight="bold")
# Alternating rows
for i in range(1, len(table_data) - 1):
    if i % 2 == 0:
        for j in range(3):
            tbl[i, j].set_facecolor("#E3F2FD")

ax8.set_title("Summary Comparison", fontsize=11, fontweight="bold", pad=4)
ax8.text(0.5, 0.01, "* single pass through 180° return arc",
         transform=ax8.transAxes, ha="center", fontsize=7, style="italic")

out_path = "erl_power_analysis.png"
plt.tight_layout(rect=[0, 0, 1, 0.975])
plt.savefig(out_path, dpi=150, facecolor="white", edgecolor="none")
print(f"\nFigure saved -> {out_path}")
