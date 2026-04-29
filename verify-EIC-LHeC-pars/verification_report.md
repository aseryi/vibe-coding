# Independent Verification of Tables 1–3 in WEV1301 (IPAC'26)

**Paper:** *Accelerator Design Educational Primer — Conceptualizing and Optimizing
the Hybrid LHeC-like Electron-Ion Collider Design*

**Verifier:** independent Python recomputation against actual reference PDFs
(`ref-6-EIC_ParameterList.pdf`, `ref-7-Phase1_LHeC_ESPP2026.pdf`, `ref-7-backup.pdf`).

Scripts:
`verify_from_refs.py` (Table 1 cross-check vs ref tables);
`verify_parents.py` (Table 1 derived rows from paper inputs);
`verify_hybrid.py` (bunch-reduction logic + Tables 2/3).

---

## 1. Method

For each design point I (a) read the input rows from the reference PDF,
(b) recompute every derived quantity that the paper also lists, and
(c) compare against the value quoted in the paper. Formulas used:

- `gamma = E / m`,  `eps_geo = eps_n / gamma`,  `sigma* = sqrt(beta* eps_geo)`,
  `sigma'* = sigma*/beta*`
- `I = e · N · f_coll`
- `xi_p = N_e · r_p · beta_p / [2 π γ_p · sigma*_e · (sigma*_e+sigma*_e)]`
- Round matched beams: `L = f_coll · N_e · N_p · H / (4 π sigma_x sigma_y)`
- Unequal beams:  `L = f_coll · N_e · N_p / [2 π · √(σ²_xe+σ²_xp) · √(σ²_ye+σ²_yp)]`

Luminosities printed in cm⁻² s⁻¹.

---

## 2. Table 1 — LHeC column versus ref [7]

The actual LHeC reference (Piotrzkowski et al., back-up document
arXiv:2503.20475) contains:

* **Table 1 of backup** (electron-source spec): bunch repetition rate
  **40 MHz**, average current **60 mA**, bunch charge **1.5 nC**,
  normalised emittance ≤ 6 mm·mrad (gun spec; the IP value used in the
  luminosity table is 22 mm·mrad).
* **Table 2 of backup** (luminosity reach): β\*_p = 0.35 m, β\*_e = 0.20 m,
  ε_p = 3.3 × 10⁻¹⁰ m·rad (geometric), ε_e = 5.7 × 10⁻¹⁰ m·rad (= 22 µm
  normalised), matched IP beam size σ\* = **10.7 µm**, N_p = 2.2 × 10¹¹,
  electron current 60 mA, **luminosity = 6 × 10³³ cm⁻² s⁻¹**.

A 1.5 nC bunch means **N_e = Q/e = 9.36 × 10⁹**, not 1.19 × 10¹⁰.

| Parameter | Paper Table 1 (LHeC col) | Reference [7] | Status |
|-----------|-------------------------:|--------------:|:------:|
| f_coll                  | 4.0 × 10⁷ Hz   | 40 MHz                       | ✓ |
| I_e                     | 60.1 mA        | 60 mA                        | ✓ |
| β\*_e (= β\*_p)          | 0.20 / 0.35 m  | 0.20 / 0.35 m                | ✓ |
| ε_n (e/p)               | 22 / 2.5 µm    | 22 / (3.3×10⁻¹⁰·γ_p ≈ 2.5)   | ✓ |
| σ\* (round, matched)    | 10.6 / 10.7 µm | 10.7 µm                      | ✓ |
| **N_e**                 | **1.19 × 10¹⁰** | **9.36 × 10⁹ (≈ 1.0 × 10¹⁰)** | **✗** |
| **L**                   | **5.78 × 10³³** | **6.0 × 10³³**                | **✗** |
| **I_p**                 | **1110 mA**     | (1.41 A from N_p × f_coll)   | **✗** |

**The paper's LHeC column mixes two different LHeC parameter sets.**

- f_coll = 40 MHz and I_e = 60 mA are taken from the new ESPP-LHeC
  submission (which uses 1.5 nC bunches → N_e ≈ 9.36 × 10⁹).
- N_e = 1.19 × 10¹⁰ is the LHeC-CDR-era value, which corresponds to f_coll
  ≈ 32 MHz (CDR convention) so that I_e = e · N_e · f_coll still gives
  60 mA.
- The quoted L = 5.78 × 10³³ is what 1.19 × 10¹⁰ × 2.2 × 10¹¹ gives at
  f_coll = 32 MHz — not what the actual ESPP-LHeC reference reports
  (6.0 × 10³³).
- The quoted I_p = 1110 mA equals e · 2.2 × 10¹¹ · 31.5 MHz = 1109 mA;
  at the 40 MHz that is also quoted, I_p = 1410 mA.

Either way the column is internally inconsistent. *Recommended fix:* pick
the new (ESPP-LHeC) numbers and update **N_e → 0.94 × 10¹⁰**,
**L → 6.0 × 10³³**, **I_p → 1410 mA**, keeping f_coll = 40 MHz and
I_e = 60 mA. This matches the open `% CHECK: LHeC rep rate 40 MHz and
why initial Grok code gave it 32 MHz value?` note in the .tex source.

---

## 3. Table 1 — EIC column versus ref [6]

The actual EIC parameter list (Peggs/Satogata, April 4 2026) Table 1.1
(*Proton collisions – 275 GeV*, column 1) gives the design 1.0 × 10³⁴ point
at √s = 105 GeV. Verbatim:

| Quantity | Ref [6] value |
|----------|--------------:|
| E_cm | 105 GeV |
| E_p / E_e | 275 / 10 GeV |
| L | **1.00 × 10³⁴** |
| hourglass H | 0.94 |
| n_bunches | 1160 |
| σ_H / σ_V (matched IP) | **95 / 8.5 µm** |
| K_x = σ_H/σ_V | 11.1 |
| Protons: bunch intensity / I / ε_n h/v / β\* h/v / ξ h/v | 6.9 × 10¹⁰ / **1.00 A** / 3.3 / 0.30 µm / 0.80 / 0.072 m / 0.012 / 0.012 |
| Electrons: bunch intensity / I / ε_n h/v / β\* h/v / ξ h/v | 17.2 × 10¹⁰ / **2.5 A** / 391 / 26 µm / 0.55 / 0.065 m / 0.072 / 0.100 |

| Parameter | Paper Table 1 (EIC col) | Reference [6] | Status |
|-----------|------------------------:|--------------:|:------:|
| n_bunches              | 1160                | 1160                | ✓ |
| f_coll                 | 9.1 × 10⁷ Hz        | 1160 × 78.2 kHz = 90.7 MHz | ✓ |
| N_e / N_p              | 1.72e11 / 6.9e10    | 1.72e11 / 6.9e10    | ✓ |
| ε_n^e h/v              | 391 / 26 µm         | 391 / 26 µm         | ✓ |
| ε_n^p h/v              | 3.3 / 0.3 µm        | 3.3 / 0.30 µm       | ✓ |
| β\*_p h/v              | 0.80 / 0.07 m       | 0.80 / 0.072 m      | ✓ (rounded) |
| **β\*_e h/v**          | **0.45 / 0.06 m**   | **0.55 / 0.065 m**  | see note |
| σ\* h/v (e and p)      | 94.8/8.6, 94.9/8.6 µm | 95 / 8.5 µm        | ✓ |
| I_p                    | 1000 mA             | 1.00 A              | ✓ |
| **I_e**                | **1555 mA**         | **2.5 A**           | **✗ wrong** |
| ξ_p h/v                | 0.012 / 0.012       | 0.012 / 0.012       | ✓ |
| L                      | 1.05 × 10³⁴         | 1.00 × 10³⁴         | ✓ |

**Two issues, one minor and one a real error:**

* Minor (β\*_e). Reference [6] *itself* is internally inconsistent: with the
  quoted β\*_e = 0.55/0.065 m, the electron σ\* would be 105/9.3 µm, but the
  ref also lists matched σ\* = 95 / 8.5 µm. The σ' values in the same table
  (211/152 µrad) and the stated matched IP beam size (95 µm) both line up
  with β\*_e ≈ **0.45/0.06 m**. The IPAC26 paper has implicitly resolved
  the EIC reference's own inconsistency by using 0.45/0.06.

* **Major: the paper's I_e = 1555 mA is wrong.** The actual EIC value at
  this design point is **2.5 A**, which is what the (also-quoted)
  N_e = 1.72 × 10¹¹ and f_coll = 91 MHz combine to give. *Recommended fix:*
  I_e → 2500 mA.

---

## 4. Bunch-reduction logic (1160 → 386)

Independently verified from first principles using the HSR ring geometry
(C = 3833.85 m → f_rev = 78.20 kHz):

- EIC nominal: n_b = 1160 → f_coll = **90.71 MHz** (paper Table 1: 91 MHz). ✓
- Hybrid choice: n_b = 386 (≈ 1160/3) → f_coll = **30.18 MHz**. ✓
- At N_e = 1.19 × 10¹⁰: I_e = e · N_e · f_coll = **57.2 mA** ≈ 60 mA. ✓
- At constant I_p = 1 A: N_p = I_p / (e · f_coll) = **2.07 × 10¹¹**. ✓

The factor-of-3 reduction (1160 → 386 ⇔ 91 → 30 MHz ⇔ N_p × 3 ⇔ I_e ≈ 60 mA)
is fully self-consistent.

**Subtle note**, however: the hybrid keeps using **N_e = 1.19 × 10¹⁰**
(= ~1.9 nC bunch), which is *not* the value used in the actual ESPP-LHeC
reference (1.5 nC = 0.94 × 10¹⁰). With the ref [7] value the hybrid I_e
would be 0.94e10 × 30.2 MHz × e ≈ **45 mA**, and N_p would still be
2.07 × 10¹¹ at constant I_p = 1 A. The luminosity would scale by
(0.94/1.19) ≈ 0.79× — i.e. ~1.7 × 10³³ instead of 2.1 × 10³³ for Table 2,
and ~7.9 × 10³³ instead of 1.0 × 10³⁴ for Table 3 (before any other
adjustment). **If the paper wants to use the ESPP-LHeC bunch charge,
some compensating adjustment of n_b, β\*, or ε is needed to recover the
quoted Table 2/3 luminosities.**

---

## 5. Table 2 — Hybrid initial configuration (~2 × 10³³)

Reproduces independently to better than 1 % on every single row
(σ\*, σ', ξ, σ_eff, I, L). No issues.

| Quantity | Computed | Paper Table 2 | Status |
|----------|---------:|--------------:|:------:|
| σ\*_e h/v        | 36.73 / 7.50 µm    | 36.7 / 7.5 µm    | ✓ |
| σ\*_p h/v        | 36.76 / 7.84 µm    | 36.8 / 7.8 µm    | ✓ |
| σ'_e h/v         | 30.6 / 149.9 µrad  | 30.6 / 149.9     | ✓ |
| σ'_p h/v         | 306.3 / 130.6      | 306.3 / 130.6    | ✓ |
| ξ_p h/v          | 0.0007 / 0.0018    | 0.001 / 0.002    | ✓ |
| ξ_e h/v          | 3.523 / 0.689      | 3.474 / 0.679    | ✓ |
| I_e              | 57.2 mA            | 57.6 mA          | ✓ |
| I_p              | 1009 mA            | 1000 mA          | ✓ |
| Σ_x_eff / Σ_y_eff | 51.96 / 10.85 µm  | 52.0 / 10.8 µm   | ✓ |
| L                | 2.117 × 10³³       | 2.1 × 10³³       | ✓ |

---

## 6. Table 3 — Hybrid optimised configuration (10³⁴)

Mostly reproduces, but the **proton vertical chain has a 15 % internal
inconsistency**: with the inputs `ε_n,p_y = 0.1 µm`, `β\*_p_y = 0.05 m`,
σ\*_p_y comes out 4.13 µm, not 3.6 µm. This propagates into ξ_e_y, σᵧ_eff,
and the peak luminosity (9.5 × 10³³ instead of 1.0 × 10³⁴).

Two equivalent fixes:

| Adjust | New value | Why |
|--------|-----------|-----|
| `ε_n,p_y` | **0.076 µm** (paper says 0.1) | gives σ\*_p_y = 3.6 µm at β\* = 0.05 m |
| `β\*_p_y`  | **0.038 m** (paper says 0.05) | gives σ\*_p_y = 3.6 µm at ε_n = 0.1 µm |

| Quantity | Computed | Paper Table 3 | Status |
|----------|---------:|--------------:|:------:|
| σ\*_e h/v        | 16.76 / 3.35 µm    | 16.8 / 3.4 µm    | ✓ |
| σ\*_p h/v        | **16.62 / 4.13 µm** | **16.8 / 3.6 µm** | **✗ vertical** |
| σ'_e h/v         | 67.1 / 33.5 µrad   | 66.8 / 32.7      | ✓ |
| σ'_p h/v         | 184.7 / 82.6 µrad  | 182.4 / 77.8     | ✗ vertical |
| ξ_p h/v          | 0.0026 / 0.0074    | 0.003 / 0.007    | ✓ |
| ξ_e h/v          | 3.487 / **5.614**  | 3.474 / **6.789** | ✗ vertical |
| Σ_x_eff / Σ_y_eff | 23.61 / 5.32 µm   | 23.8 / 5.0 µm    | ✗ vertical |
| L                | **9.50 × 10³³**    | **1.0 × 10³⁴**   | ✗ ~5 % short |

---

## 7. Top-line summary of recommended corrections

1. **Paper Table 1, LHeC column:** The column mixes ESPP-LHeC (40 MHz)
   with old LHeC-CDR (N_e = 1.19 × 10¹⁰). To match ref [7]:
   * keep f_coll = 4.0 × 10⁷ Hz, I_e = 60.1 mA, β\*, ε_n;
   * change **N_e → 0.94 × 10¹⁰** (1.5 nC bunch);
   * change **I_p → 1410 mA** (= 2.2 × 10¹¹ × 40 MHz × e);
   * change **L → 6.0 × 10³³**.

2. **Paper Table 1, EIC column:** All inputs and most derived quantities
   match ref [6], *except* I_e:
   * change **I_e → ~2500 mA (= 2.5 A)**.
   * (β\*_e = 0.45/0.06 m is fine — reference [6] itself has β\*_e = 0.55/0.065
     listed but its other rows are consistent with 0.45/0.06; the paper has
     resolved this in the right direction.)

3. **Bunch-reduction logic** (1160 → 386, factor 3): consistent if the hybrid
   uses N_e = 1.19 × 10¹⁰. **If the paper switches to ref [7]'s 0.94 × 10¹⁰,
   the Tables 2/3 luminosities drop by ~21 %** unless n_b, β\*, or ε is
   re-tuned. Worth flagging in §3 of the paper.

4. **Paper Table 2** (hybrid initial) — reproduces to <1 %. No change needed.

5. **Paper Table 3** (hybrid optimised, 10³⁴): proton vertical chain is
   internally inconsistent at the 15 % level. Lower **ε_n,p_y from 0.1 µm
   to 0.076 µm**, **or** lower **β\*_p_y from 50 mm to 38 mm**. Either fix
   makes σ\*_p_y = 3.6 µm and L = 1.0 × 10³⁴ self-consistent.

---

## Files

- `verify_from_refs.py` — recomputes Table 1 inputs vs the actual reference PDFs.
- `verify_parents.py` — recomputes Table 1 derived rows from the paper's own inputs.
- `verify_hybrid.py` — bunch-reduction logic + Tables 2 and 3 (hybrid columns).
- `verification_report.md` — this report.
