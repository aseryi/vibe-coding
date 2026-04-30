# ERL Simulation — Model Assumptions & Discussion Guide

**File:** `erl_simulation.py` | **Repository:** [github.com/aseryi/vibe-coding](https://github.com/aseryi/vibe-coding) (EIC-LHeC folder)  
**Purpose:** Educational simulation of single-pass and two-pass Energy Recovery Linacs (ERLs) for electron–proton collisions at EIC-like parameters.  
**Audience:** Graduate students and early-career accelerator physicists.

---

**Legend:**  
✅ **GIVEN** — specified in the original prompt  
🟡 **ASSUMED** — chosen by the author / AI; could be varied  
🔵 **DERIVED** — follows from other quantities by formula  
🔴 **MODEL** — a deliberate simplification of the physics

---

## 1. Electron Beam Parameters

| Parameter | Value | Status | Assumption / Justification | Ref |
|-----------|-------|--------|---------------------------|-----|
| `Q_bunch` | 1.5 nC | ✅ GIVEN | Specified in the prompt. Representative of high-current ERL designs (e.g. Cornell ERL, eRHIC). Sets the average current together with the repetition rate. | [1, 2] |
| `sigma_z` | 1.5 mm (rms) | ✅ GIVEN | Specified in the prompt. Short bunches are needed for high luminosity (hourglass factor) and to limit coherent SR in arcs. Achievable with velocity bunching or a magnetic chicane compressor at low energy. | [1, 3] |
| `E_inj` | 7 MeV | ✅ GIVEN | Typical photoinjector output energy: high enough for the beam to be fully relativistic (γ ≈ 14) and for space-charge forces to be manageable, yet low enough to limit the dump power P = I·E_inj. | [4, 5] |
| `delta_in_rms` | 4×10⁻³ | 🟡 ASSUMED | Initial fractional energy spread at injection. Typical for a DC photogun followed by a bunching section and booster cavity. **Discussion:** What happens to the final dump energy spread if this is 10× larger? Smaller? | [4] |
| `N_e` | 9.36×10⁹ | 🔵 DERIVED | N_e = Q / e. No loss is assumed between injection and the IP (100% transmission). In reality 1–5% beam loss budgets are typical. | — |
| `eps_n_e` (1-pass) | 20 μm·rad | 🟡 ASSUMED | Normalised transverse emittance at the collision point. Set by the photocathode quantum emittance and solenoid transport, not by SR damping. Deliberately conservative; advanced guns achieve < 1 μm·rad at lower charge. | [4, 6] |
| `eps_n_e` (2-pass) | 10 μm·rad | 🟡 ASSUMED | Assigned a factor-of-two improvement for the two-pass case, motivated by the longer acceleration path allowing better phase-space control. This is speculative; both cases use the same injector. | — |

> **Student exercise:** The ERL emittance is gun-limited (~10–20 μm·rad), while a 10 GeV storage ring has SR equilibrium emittance ~0.1 μm·rad. Why does the ERL still achieve competitive or higher luminosity? What compensates for the larger emittance?

---

## 2. RF and Timing Parameters

| Parameter | Value | Status | Assumption / Justification | Ref |
|-----------|-------|--------|---------------------------|-----|
| `f_RF` | 801.58 MHz | ✅ GIVEN | Specified in the prompt. Close to 800 MHz, a common choice for high-current ERLs (e.g. eRHIC design studies). Lower than 1.3 GHz TESLA cavities — better HOM damping at high current, larger bore radius. | [2, 7] |
| `f_rep` | 40.08 MHz | ✅ GIVEN | Specified. Sets average current I = Q × f_rep = 60 mA. Chosen as 1/20 of f_RF, so every 20th RF bucket carries a bunch (harmonic = 20). This sparse filling reduces HOM excitation. | [2] |
| `harmonic` | 20 | 🔵 DERIVED | h = f_RF / f_rep = 801.58 / 40.08 ≈ 20.00. Exact integer implies synchronous timing between RF and bunch arrival — a necessary operating condition. | — |
| `I_arc` | 60 mA | 🔵 DERIVED | Current in the arcs (one-way beam). Each arc carries only the accelerating *or* decelerating bunch train — not both. | — |
| `I_linac` | 120 mA | 🔵 DERIVED | The SRF linac is traversed by both the accelerating and decelerating bunch trains simultaneously. Their beam-loading currents are anti-phased and largely cancel — the key energy-recovery mechanism. Any mismatch drives residual generator power requirements. | [8] |
| `phi_acc` | 0° (crest) | 🟡 ASSUMED | Maximum energy gain. In practice a small off-crest angle (~5°) provides longitudinal focusing to prevent phase-space growth in multi-cell cavities. The code does not model this correction. | [9] |

---

## 3. SRF Linac Parameters

| Parameter | Value | Status | Assumption / Justification | Ref |
|-----------|-------|--------|---------------------------|-----|
| `G_acc` | 15 MV/m | 🟡 ASSUMED | State-of-the-art for production SRF cavities at 800 MHz. Current CEBAF cavities achieve 19–22 MV/m; 15 MV/m is a conservative, reliable operating gradient leaving headroom for quench margin. **Discussion:** Increasing to 20 MV/m would reduce cavity count by 25% and cryoplant power proportionally. | [9, 10] |
| `L_cell` | 0.187 m | 🔵 DERIVED | L_cell = λ/2 = c / (2 f_RF). Standard pillbox-cell geometry for a β = 1 cavity. | [9] |
| `N_cells_cav` | 5 | 🟡 ASSUMED | Five-cell cavity, a common choice balancing HOM spectrum (fewer modes than 7-cell) against cost per unit voltage. CEBAF uses 5-cell; TESLA uses 9-cell at 1.3 GHz. | [9, 11] |
| `V_cav` | ≈ 14 MV | 🔵 DERIVED | V_cav = G × L_cav = 15 MV/m × 0.935 m. | — |
| `Q_0` | 2×10¹⁰ | 🟡 ASSUMED | Intrinsic quality factor. Routinely achieved in modern Nb cavities at 2 K. Nb₃Sn and Nb/Cu coatings are targeting Q₀ > 10¹¹ for future machines; that would reduce cryoplant power by > 5×. | [9, 12] |
| `R/Q` | 400 Ω | 🟡 ASSUMED | Geometric shunt impedance per cavity (linac convention: R/Q = V² / (ω U)). Estimated for a 5-cell 802 MHz cavity by scaling from CEBAF and TESLA designs. Determines dynamic heat load: P_2K = V² / (R/Q · Q₀). | [9, 10] |
| `eta_cryo` | 250 W_plug/W_2K | 🟡 ASSUMED | Cryoplant coefficient of performance (inverse). Carnot limit at 2 K is T_hot/ΔT ≈ 150; realistic plants are ~1.5–2× Carnot, giving 200–300 W/W. 250 W/W is a standard engineering estimate. | [13] |
| `eta_RFPA` | 65% | 🟡 ASSUMED | RF power amplifier (klystron or IOT) wall-plug efficiency. Modern klystrons reach 65–70%. Solid-state amplifiers can reach > 80% but are not yet commercially available at the required power levels for SRF. | [13] |
| `eta_inj` | 35% | 🟡 ASSUMED | Injector wall-plug efficiency (DC gun + bunching + booster). Normal-conducting guns are typically 30–40% efficient. This is the dominant power cost in the ERL budget. Improving this (e.g. using a superconducting gun) would directly reduce total wall-plug power. | [4, 5] |
| static heat load | neglected | 🔴 MODEL | The code includes only dynamic (RF-driven) heat load. Static heat from thermal radiation, support conduction, and HOM deposited power can add 20–50% to the cryogenic load in a real system. | [9] |

> ⚠️ **Key limitation:** The cryoplant power (4.4 MW for 1-pass ERL) is the largest single budget item, driven by the assumed Q₀ and R/Q. Students should examine how this scales: P_cryo ∝ V² / (R/Q · Q₀ · η_cryo). Doubling Q₀ halves the cryoplant cost.

---

## 4. Longitudinal Dynamics Model

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Thin-lens RF kick | 🔴 MODEL | Each accelerating pass is treated as a single instantaneous kick: ΔE = eV₀[cos(φ_s + φ) − cos(φ_s)]. This is exact in the limit where the bunch is much shorter than the cavity length, which holds here (σ_z = 1.5 mm ≪ L_cav = 935 mm). The full linac (hundreds of cavities) is collapsed into a single kick per pass. | [14, 15] |
| No longitudinal focusing | 🔴 MODEL | A linac has no restoring force in the longitudinal plane (unlike a storage ring where synchrotron oscillations confine the bunch). The RF phase is fixed; particles with different phases receive different kicks but are not returned to the reference phase. This is correct for a linac; the model therefore does not include synchrotron oscillations. | [14] |
| β ≈ 1 throughout | 🔴 MODEL | Electrons are ultra-relativistic at all energies considered (γ ≥ 14 at 7 MeV injection). Phase slippage between different-energy particles is therefore negligible and is not modelled. This would break down below ~2 MeV. | [14] |
| No wake fields / HOMs | 🔴 MODEL | Resistive-wall and cavity higher-order-mode (HOM) wake fields can drive emittance growth and the beam-breakup (BBU) instability in ERLs. These are entirely absent from the simulation. BBU is a critical design constraint for high-current ERLs and is the subject of dedicated codes (e.g. TDBBU, BMAD). | [16, 17] |
| No space charge | 🔴 MODEL | Space-charge forces scale as 1/γ² and are negligible for γ ≫ 1. They are important only at very low energy (< few MeV) near the cathode and gun exit, which is before injection into the simulation. | [14] |
| No coherent SR in arcs | 🔴 MODEL | Coherent synchrotron radiation (CSR) in the bending magnets of the return arc can cause longitudinal emittance growth for short bunches. With σ_z = 1.5 mm at 10 GeV, CSR is moderate but non-negligible. The arc SR power is included in the power budget (incoherent average), but CSR-induced emittance growth is not tracked. | [18] |
| Gaussian initial distribution | 🔴 MODEL | Both the longitudinal phase (φ) and energy deviation (δE) at injection are drawn from independent Gaussian distributions. Real photocathode guns produce distributions closer to a truncated Gaussian with non-Gaussian tails. Tails matter for halo formation and machine protection. | [4] |
| 50,000 macro-particles | 🟡 ASSUMED | Sufficient to resolve the core of the distribution at the 10⁻⁴ level. The simulation is not attempting to track beam halo (would need > 10⁶ particles and importance sampling). | — |
| Fixed random seed (42) | 🟡 ASSUMED | Ensures reproducibility. With 50,000 particles statistical noise is small (< 0.5% on rms quantities) and results are seed-independent. | — |

---

## 5. Return Arc and 180° Phase Reversal

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Ideal phase flip | 🔴 MODEL | The return arc is assumed to impart *exactly* a 180° phase advance relative to the RF. In practice this is achieved by tuning the arc path length (and hence its R₅₆) so that the net phase at the linac entrance shifts by (2n+1)π for some integer n. Any deviation is the "phase error" studied in Figure 4. | [8, 19] |
| Linear arc transport | 🔴 MODEL | The arc is modelled as a perfect phase shifter with no emittance growth, chromatic aberrations, or dispersion mismatch. In a real design the arc optics must be carefully matched to suppress second-order dispersion (T₅₆₆) which would smear the phase distribution of off-momentum particles. | [19] |
| Single arc (1-pass) / two arcs (2-pass) | 🔴 MODEL | The single-pass ERL has one return arc from the IP back to the linac entrance. The two-pass ERL has two arcs of different energies on each side (acceleration and deceleration), for a total of four arc segments. Only the high-energy arc SR is included in the power budget; lower-energy arcs contribute negligibly. | [1, 2] |
| Correlated phase error | 🔴 MODEL | In the efficiency-vs-phase-error analysis (Fig. 4), the same phase error Δφ is applied to every deceleration pass. In a two-pass ERL this may be conservative (systematic error) or optimistic (if errors are independent per pass). | — |
| `rho_arc` for SR | 300 m | 🟡 ASSUMED | Bending radius used for the SR energy-loss calculation in the arc. A larger radius reduces SR losses but increases arc length and cost. For EIC-scale machines, radii of 100–300 m are plausible. **Discussion:** How does the SR arc power scale with energy and radius? (P ∝ E⁴ / ρ). | [20] |
| Arc SR fraction | 0.5 (half-turn) | 🔴 MODEL | SR energy loss is approximated as U₀/2, treating the return arc as sweeping half of a full ring. This is a rough geometric estimate; the actual loss depends on the arc bending-magnet layout. | — |

---

## 6. Energy Recovery Efficiency

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Analytical formula η = cos(Δφ) | 🔴 MODEL | For a reference particle accelerated at the RF crest (φ_acc = 0) and decelerated at π + Δφ: E_recovered = E_gain × cos(Δφ), so η = cos(Δφ). This is exact for a mono-energetic beam; the finite energy spread of the bunch broadens the effective efficiency curve slightly. | [8] |
| 100% beam transmission | 🔴 MODEL | All electrons that are accelerated are assumed to be decelerated and reach the dump. Beam loss along the return path (halo, Touschek scattering, gas scattering) would reduce the recovered power and deposit radioactivating energy in the tunnel. Typical ERL loss budgets target < 1 W total beam loss. | [1] |
| Residual RF = 0.1% of beam power | 🟡 ASSUMED | Generator RF power required to compensate imperfect beam-loading cancellation between the accelerating and decelerating bunches. In an ideal ERL these cancel exactly. The 0.1% residual is a rough engineering estimate; the actual value depends on phase and amplitude stability of the RF control system (LLRF). | [8, 21] |

> **Student exercise:** The efficiency formula η = cos(Δφ) is independent of the number of passes. Why? Sketch an argument. Would this change if the phase error were different in each pass?

---

## 7. Proton Beam Parameters (EIC-like)

| Parameter | Value | Status | Assumption / Justification | Ref |
|-----------|-------|--------|---------------------------|-----|
| `E_p` | 275 GeV | 🟡 ASSUMED | Maximum proton energy in the EIC Hadron Storage Ring (RHIC-derived). Sets √s ≈ 105 GeV for 10 GeV electrons. | [22] |
| `N_p` | 6.88×10¹⁰ | 🟡 ASSUMED | EIC design bunch population for 275 GeV protons. This drives the large beam-beam parameter ξ_e ~ 3.5 for the electron beam, which would be catastrophic in a storage ring but is acceptable in a single-pass ERL because there is no resonance accumulation. | [22] |
| `eps_n_p` | 1 μm·rad | 🟡 ASSUMED | EIC design normalised proton emittance. Smaller than the electron emittance in this simulation, so the proton beam is the narrower of the two at the IP. Proton emittance can be reduced further by coherent electron cooling (CeC). | [22, 23] |
| `beta_star` | 5 cm (round) | 🟡 ASSUMED | Equal β* in x and y, and equal for electrons and protons. This gives the simplest round-beam IP geometry. Real EIC designs use flat beams (β*_x ≠ β*_y) to optimise luminosity vs. aperture constraints. | [22] |
| `f_rep_p = f_rep_e` | — | 🟡 ASSUMED | Both beams are assumed to have the same bunch repetition rate (40.08 MHz), ensuring every electron bunch collides with exactly one proton bunch. In practice the proton ring circumference may not be an exact integer multiple of the electron bunch spacing, requiring careful RF gymnastics. | [22] |

---

## 8. Beam-Beam Interaction Model

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Round-beam, strong-weak ξ | 🔴 MODEL | The beam-beam tune-shift parameter is computed in the strong-weak (test-particle) approximation for round beams: ξ = r_e · N_opp / (4π · γ · σ²). The convolved beam size σ_coll = √(σ_e² + σ_p²) is used. The weak-strong approximation ignores the back-action on the strong (proton) beam from the weaker (electron) beam. | [24, 25] |
| ξ_e unconstrained in ERL | 🔴 MODEL | In a storage ring, ξ > 0.05–0.10 drives resonant emittance growth over many turns and is catastrophic. In a single-pass ERL, each electron makes only *one collision*, so there is no turn-by-turn resonance. Large ξ is therefore tolerable up to the point where the beam is geometrically disrupted. The relevant limit is then the disruption parameter D. | [24, 26] |
| Hourglass factor | 🔴 MODEL | F_hg ≈ 1/√(1 + (σ_z/β*)²). This is an approximation to the exact Gaussian integral. With σ_z = 1.5 mm and β* = 5 cm, the correction is < 0.05% (F ≈ 0.9996) — negligible here. It would matter for σ_z ~ β*. | [27] |
| Pinch enhancement | 🔴 MODEL | H_pinch ≈ 1 + 0.5 D^0.45 (empirical fit from Yokoya & Chen). When D > 1 the opposing beam focuses ("pinches") the other beam during the collision, enhancing instantaneous luminosity. The exponent 0.45 and coefficient 0.5 are fits to simulation data and have significant uncertainty. | [24, 25] |
| Beamstrahlung energy spread | 🔴 MODEL | The energy spread added to the electron beam by beamstrahlung (synchrotron radiation in the proton field) is approximated as a Gaussian smearing with width δ_bs ≈ ξ_e × (σ_z/β*) × E_beam/2. This is a very rough linear parametrisation. A full treatment requires the Yokoya-Chen formalism integrated over the bunch longitudinal profile. | [24, 28] |
| No crossing angle | 🔴 MODEL | Head-on collisions are assumed. A finite crossing angle (needed to separate the beams after the IP) reduces the geometric luminosity by the Piwinski factor F = 1/√(1 + (φ σ_z / 2σ*)²) and is not included. | [27] |
| Luminosity formula | 🔴 MODEL | L = H · N_e · N_p · f_rep / (4π σ_e σ_p). This uses the product of individual beam sizes in the denominator rather than the sum of squares. For an asymmetric collision the exact formula uses (σ_e² + σ_p²); here σ_e σ_p ≠ (σ_e² + σ_p²) by a factor of order 2 when beams have comparable sizes. Luminosity values are therefore approximate to within a factor ~2. | [27] |

---

## 9. Power Budget Assumptions

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Dump energy = E_inj (ideal) | 🔴 MODEL | In the ideal ERL all accelerated energy is returned to the RF and only the injection energy (7 MeV) is lost at the dump. P_dump = I × E_inj = 421 kW. Any phase error or beam loss increases this. | [8] |
| Beam loss = 0 | 🔴 MODEL | All electrons reach the dump. Even 10 W of lost beam at 10 GeV activates materials significantly and is a hard engineering limit, but the *power* impact on the budget is negligible at the 0.01% level. | [1] |
| 2-pass cryoplant ∝ V² | 🔵 DERIVED | The two-pass ERL uses one linac at half the voltage (5 GeV instead of 10 GeV). Dynamic SRF losses scale as V², so the cryoplant power is 4× smaller: a major advantage of multi-pass operation. | [9] |
| Storage ring τ_beam = 2 hr | 🟡 ASSUMED | Beam lifetime used to estimate injector power for the storage ring comparison. A 10 GeV ring will be Touschek-lifetime-limited; 2 hours is plausible for moderate bunch intensity but depends strongly on the ring optics (momentum aperture) and vacuum pressure. | [20] |
| Ring magnet power ~1 kW/m | 🟡 ASSUMED | Rough estimate for normal-conducting dipole magnets. For a compact ring (C ~ 2π × 50 × 4 ≈ 1260 m) this gives ~1.3 MW, but is a very coarse model. Superconducting arc magnets would have negligible resistive loss. | — |
| Ring RF efficiency η_RF = 50% | 🟡 ASSUMED | Wall-plug efficiency of the entire ring RF system (klystrons + waveguides + cavities). Comparable to the ERL RFPA assumption (65%), reduced by the additional normal-conducting cavity losses in a ring RF system. | [13] |

---

## 10. Two-Pass ERL Specific Assumptions

| Assumption | Status | Details and Limitations | Ref |
|------------|--------|------------------------|-----|
| Same single linac, 4 passes | 🔴 MODEL | One SRF linac section is traversed four times per electron: acceleration pass 1, acceleration pass 2, deceleration pass 1, deceleration pass 2. At any instant the linac contains bunches at four different energies. Net beam loading cancels between the two accelerating and two decelerating passes. The linac cavity count is therefore half that of the 1-pass case (357 vs. 714). | [1, 26] |
| Independent arc phases | 🔴 MODEL | In a two-pass ERL each recirculation arc must deliver a distinct phase advance to route the beam to the correct linac entrance. In practice a spreader/recombiner magnet system separates beams of different energies at the linac ends. This complexity is not modelled; only the energy and phase at each linac traversal are tracked. | [1, 29] |
| eps_n_e (2-pass) = 10 μm·rad | 🟡 ASSUMED | Assigned 2× smaller than the 1-pass case. In reality both start from the same gun, and any emittance improvement in the 2-pass case would come from better phase-space manipulation, not from the recirculation itself. This assumption is optimistic and should be treated as an upper-bound estimate of 2-pass luminosity gain. | — |

---

## 11. What the Model Does NOT Include

The following physical effects are absent from the simulation and represent the primary areas for future improvement:

- **Beam-breakup (BBU) instability** — driven by HOM wake fields, this is the primary current limit for multi-pass ERLs. References: [16, 17].
- **Touschek scattering and intra-beam scattering (IBS)** — relevant to emittance growth and beam lifetime in the arcs. Reference: [20].
- **Coherent synchrotron radiation (CSR) emittance growth** — particularly important for short (< 1 mm) bunches in the arcs. Reference: [18].
- **Nonlinear optics and dynamic aperture** — chromatic aberrations in the arcs limit the energy acceptance (typically ±0.5–2%).
- **Spin polarisation** — required for many EIC physics programmes; spin rotation and Siberian snakes are not included.
- **Crab cavities** — needed to restore head-on geometry when a crossing angle is used.
- **Beam-gas scattering and residual-gas ionisation** — contribute to emittance growth and beam loss at the 10⁻⁹–10⁻¹⁰ Torr vacuum required.
- **Radiation damage and activation** — from any beam loss, not modelled.

---

## 12. Suggested Student Exercises

1. Vary the SRF quality factor Q₀ from 10⁹ to 10¹¹ and plot the total wall-plug power. At what Q₀ does the cryoplant cease to dominate?
2. Repeat the efficiency analysis with a phase error that grows linearly with pass number (independent errors). How does two-pass compare to one-pass?
3. Replace the round-beam hourglass formula with the exact integral and compare to the approximate result as a function of σ_z/β*.
4. Add a simple beam-loss term (e.g. 0.1% loss per pass at 10 GeV) and recalculate the dump energy and total wall-plug power.
5. Model the SRF static heat load as a fixed fraction (30%) of the dynamic load and determine the new total cryoplant cost.
6. Implement the full Yokoya-Chen beamstrahlung energy-loss distribution (instead of the Gaussian smearing) and compare the post-collision energy spread.
7. Estimate the BBU threshold current using the simplified single-HOM model and compare it to the operating current of 60 mA. (See [16, 17].)
8. Design a two-pass ERL with 6 GeV/pass (reaching 12 GeV) and compare the power budget to the 5 GeV/pass case modelled here.

---

## References

1. Hoffstaetter, G.H. et al., "CBETA: Cornell-BNL ERL Test Accelerator," *IPAC 2017*, WEYAA1 (2017). [WEYAA1](https://accelconf.web.cern.ch/ipac2017/papers/weyaa1.pdf)
2. Ben-Zvi, I. et al., "eRHIC Design Study," BNL-94022-2011 (2011).
3. Gruner, S.M. et al., "Considerations for a Soft X-Ray Free Electron Laser," *Rev. Sci. Instrum.* **73**, 1402 (2002). DOI: [10.1063/1.1445799](https://doi.org/10.1063/1.1445799)
4. Bazarov, I.V. et al., "Benchmarking of codes for electron beam brightness and emittance growth from a photoinjector," *PRST-AB* **11**, 100703 (2008). DOI: [10.1103/PhysRevSTAB.11.100703](https://doi.org/10.1103/PhysRevSTAB.11.100703)
5. Poelker, M., "DC photoguns for ERLs," *ERL Workshop* 2007.
6. Karkare, S. et al., "Ultracold Electrons via Near-Threshold Photoemission from Single-Crystal Cu(100)," *PRL* **125**, 054801 (2020). DOI: [10.1103/PhysRevLett.125.054801](https://doi.org/10.1103/PhysRevLett.125.054801)
7. Calaga, R. et al., "High-current SRF cavity design at 802 MHz for the eRHIC ERL," *SRF 2009*, THPPO077.
8. Neil, G.R. and Merminga, L., "Technical approaches for high-average-power free-electron lasers," *Rev. Mod. Phys.* **74**, 685 (2002). DOI: [10.1103/RevModPhys.74.685](https://doi.org/10.1103/RevModPhys.74.685)
9. Padamsee, H., Knobloch, J. and Hays, T., *RF Superconductivity for Accelerators*, 2nd ed., Wiley-VCH (2008). ISBN 978-3-527-40842-9.
10. Liepe, M. et al., "SRF experience at CEBAF," *SRF 2011*, MOIOA05.
11. Aune, B. et al., "Superconducting TESLA cavities," *PRST-AB* **3**, 092001 (2000). DOI: [10.1103/PhysRevSTAB.3.092001](https://doi.org/10.1103/PhysRevSTAB.3.092001)
12. Posen, S. and Hall, D.L., "Nb₃Sn superconducting radiofrequency cavities," *Supercond. Sci. Technol.* **30**, 033004 (2017). DOI: [10.1088/1361-6668/30/3/033004](https://doi.org/10.1088/1361-6668/30/3/033004)
13. Faircloth, D. and Hitit, S., "RF power sources for particle accelerators," *CERN Yellow Report* CERN-2011-007, pp. 211–247 (2011).
14. Wiedemann, H., *Particle Accelerator Physics*, 4th ed., Springer (2015). ISBN 978-3-319-18317-6.
15. Reiser, M., *Theory and Design of Charged Particle Beams*, 2nd ed., Wiley-VCH (2008). ISBN 978-3-527-40741-5.
16. Hoffstaetter, G.H. and Bazarov, I.V., "Beam-breakup instability theory for energy recovery linacs," *PRST-AB* **7**, 054401 (2004). DOI: [10.1103/PhysRevSTAB.7.054401](https://doi.org/10.1103/PhysRevSTAB.7.054401)
17. Tennant, C. et al., "Beam breakup in the Jefferson Laboratory Free Electron Laser upgrade driver," *PRST-AB* **8**, 074403 (2005). DOI: [10.1103/PhysRevSTAB.8.074403](https://doi.org/10.1103/PhysRevSTAB.8.074403)
18. Derbenev, Y.S. and Shiltsev, V., "Transverse effects of microbunch radiative interaction," FERMILAB-TM-1974 (1996).
19. Tigner, M., "A possible apparatus for electron-clashing-beam experiments," *Nuovo Cimento* **37**, 1228 (1965). [First ERL proposal.] DOI: [10.1007/BF02773204](https://doi.org/10.1007/BF02773204)
20. Sands, M., *The Physics of Electron Storage Rings*, SLAC-121 (1970). [SLAC-121](https://www.slac.stanford.edu/pubs/slacreports/slac-r-121.html)
21. Lonza, M., "Multi-bunch feedback systems," *CAS RF for Accelerators*, CERN-2011-007, pp. 467–510 (2011).
22. Khalek, R.A. et al. (EIC Yellow Report), "Science Requirements and Detector Concepts for the Electron-Ion Collider," *Nuclear Physics A* **1026**, 122447 (2022). DOI: [10.1016/j.nuclphysa.2022.122447](https://doi.org/10.1016/j.nuclphysa.2022.122447)
23. Litvinenko, V.N. and Derbenev, Y.S., "Coherent electron cooling," *PRL* **102**, 114801 (2009). DOI: [10.1103/PhysRevLett.102.114801](https://doi.org/10.1103/PhysRevLett.102.114801)
24. Yokoya, K. and Chen, P., "Beam-beam phenomena in linear colliders," in *Frontiers of Particle Beams: Intensity Limitations*, Lecture Notes in Physics 400, Springer (1992), pp. 415–445.
25. Bassetti, M. and Erskine, G.A., "Closed expression for the electrical field of a two-dimensional Gaussian charge," CERN-ISR-TH/80-06 (1980).
26. Litvinenko, V.N., "ERL-based electron-ion colliders," *IPAC 2010*, WEXRA02 (2010).
27. Hirata, K., "Analysis of beam–beam interaction with a large Piwinski angle," *PRST-AB* **67**, 1901 (1987).
28. Chen, P., "An introduction to beamstrahlung and disruption," in *Frontiers of Particle Beams*, Lecture Notes in Physics 296, Springer (1988), pp. 495–532.
29. Trbojevic, D. et al., "FFAG lattice for the multi-pass ERL," *IPAC 2012*, TUPPR006 (2012).

---

*Generated by **Claude Sonnet 4.6** (Anthropic) as part of the EIC-LHeC mini-project educational exercise — vibe-coding with AI.*  
*Repository: [github.com/aseryi/vibe-coding](https://github.com/aseryi/vibe-coding) | Code file: `EIC-LHeC/erl_simulation.py`*
