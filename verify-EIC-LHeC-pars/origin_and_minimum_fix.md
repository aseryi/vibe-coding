# Origins of inconsistencies and minimum fix proposal

Based on the Grok prompt history (`Prompts to Grok AI July 2025-1.pdf`), the
two reference PDFs, and the verification work in `verify_from_refs.py`,
`verify_parents.py`, `verify_hybrid.py`.

---

## 1. Where each inconsistency came from

### 1.1 LHeC column in Table 1 — three linked issues

**What Grok was given.** In summer 2025 the user pointed Grok at the LHeC
backup PDF (the same arXiv:2503.20475 we now have as `ref-7-backup.pdf`).
That document, in Table 1 of its electron-source spec, lists:

> *bunch repetition rate 40 MHz, average current 60 mA, bunch charge 1.5 nC*

so an honest reading of this row gives **N_e = Q/e = 9.36 × 10⁹ ≈ 1.0 × 10¹⁰**.
The luminosity Table 2 of the same document then quotes **L = 6 × 10³³**.

**What ended up in Table 1 of the IPAC paper.** The numbers actually printed
(N_e = 1.19 × 10¹⁰, I_e = 60.1 mA, I_p = 1110 mA, L = 5.78 × 10³³) form a
*different* internally consistent set: they correspond to
**N_e = 1.19 × 10¹⁰ at f_coll ≃ 32 MHz** (the original LHeC CDR Phase-1
convention from arXiv:1206.2913), not to the 40 MHz now quoted in the
collision-frequency row.

Independent check (from `verify_from_refs.py`):

| f_coll | × N_e (1.19×10¹⁰) | gives I_e |
|--------|-------------------|-----------|
| 32 MHz | × 1.19 × 10¹⁰     | **61 mA**  → matches Table 1's 60.1 mA |
| 40 MHz | × 1.19 × 10¹⁰     | 76 mA  → does *not* match  |
| 40 MHz | × 0.94 × 10¹⁰     | **60 mA** → matches if N_e is the ref [7] value |

So **Grok mixed two LHeC conventions**: the new ESPP submission
(40 MHz / 1.5 nC = 0.94 × 10¹⁰) supplied the f_coll, while the older LHeC
CDR Phase-1 set (32 MHz / 1.91 nC = 1.19 × 10¹⁰) supplied N_e and L.

This is the same bug that the original `% CHECK: LHeC rep rate 40 MHz and
why initial Grok code gave it 32 MHz value?` comment in the .tex source
already suspected.

### 1.2 §3 text — "~1.45 A electron current"

The paper's §3 says:

> *"naively combining the two implies ∼1.45 A electron current, far beyond
> the 60 mA LHeC design value"*

Recompute (`verify_from_refs.py` again):

|                   | × EIC 91 MHz | × LHeC 40 MHz |
|-------------------|--------------|---------------|
| N_e = 1.0 × 10¹⁰   | **136 mA**   | 60 mA         |
| N_e = 1.19 × 10¹⁰  | 173 mA       | 76 mA         |
| N_e = 1.0 × 10¹¹   | **1.45 A**   | 0.6 A         |

So **the "1.45 A" claim has a 10× factor error**: it would require
N_e ≈ 10¹¹, but the LHeC bunch population is only ≈ 10¹⁰. The correct
sentence should read "**~145 mA, more than 2× the 60 mA LHeC design**" — the
qualitative point still stands but the number is off by 10×.

### 1.3 EIC column in Table 1 — I_e = 1555 mA

EIC reference [6] Table 1.1, column 1, lists I_e = **2.5 A** at 10 GeV with
1160 bunches (consistent with N_e × f_coll × e). With f_coll = 91 MHz and
N_e = 1.72 × 10¹¹ (correctly quoted in the IPAC paper), the only consistent
current is 2.5 A.

**Where 1555 mA came from.** With f_coll = 91 MHz, 1555 mA implies
N_e ≈ 1.07 × 10¹¹ (≈ 17 nC bunches). 17 nC is not a standard EIC value;
it's about 60 % of the 28 nC value quoted in ref [6]. The most plausible
explanation is that the March 2025 version of the EIC parameter list (which
Grok was pointed at — that's the version in the prompt history, and the
current April 2026 version is marked "UNDER CONSTRUCTION") had different
electron-bunch numbers than the current version, **or** the value was
manually transcribed in error. Either way the row is now contradicted by the
rest of the same row.

### 1.4 Table 3 — proton vertical chain

`σ*_p_y = 3.6 µm`, `β*_p_y = 0.05 m`, `ε_n,p_y = 0.1 µm` cannot all be true
simultaneously. With β* = 0.05 m and ε_n = 0.1 µm we get σ* = 4.13 µm. Most
likely Grok started from σ_p_y = 3.6 µm (set by the matched-beam target),
and *separately* declared "factor-of-3 emittance reduction → 0.1 µm",
without re-deriving σ from ε. The actual ε that yields σ = 3.6 µm at
β* = 0.05 m is **0.076 µm** (call it 0.08).

---

## 2. Minimum fix to make the paper self-consistent

I recommend touching exactly **5 numbers and one sentence**. All four
tables, the hybrid story in §3, and the bunch-reduction logic survive.

### Fix #1 (Table 1, LHeC column) — three numbers

Anchor the column to ref [7] and keep f_coll = 4.0 × 10⁷ Hz (which is
already there).

| Row | Currently | Change to | Reason |
|-----|-----------|-----------|--------|
| `e` bunch population | 1.19 × 10¹⁰ | **1.0 × 10¹⁰** | matches ref [7]'s 1.5 nC |
| `p` beam current     | 1110 mA     | **1410 mA**    | = e · 2.2 × 10¹¹ · 40 MHz |
| Peak ℒ                | 5.78 × 10³³ | **6.0 × 10³³** | matches ref [7] |

`I_e = 60.1 mA`, `√s = 749 GeV`, σ\*, ξ, β\*, ε all stay the same and are
correct.

### Fix #2 (Table 1, EIC column) — one number

| Row | Currently | Change to | Reason |
|-----|-----------|-----------|--------|
| `e` beam current | 1555 mA | **2500 mA** | = e · 1.72 × 10¹¹ · 91 MHz; matches ref [6] (2.5 A) |

Everything else in the column reproduces.

### Fix #3 (Table 3, optimised hybrid) — one number

Either change ε or change β\*; ε is the cleaner edit because it lives in
a row that's also already different from Table 2.

| Row | Currently | Change to | Reason |
|-----|-----------|-----------|--------|
| `p` Norm. ε h/v | 0.9 / 0.1 µm | **0.9 / 0.08 µm** | gives σ\*_p_y = 3.6 µm at β\* = 0.05 m |

After this single edit Table 3 becomes self-consistent: σ\*_p_y = 3.60 µm,
σ'\_p_y = 80 µrad, ξ_e_y ≈ 6.6, σᵧ_eff = 5.0 µm, peak ℒ = 1.0 × 10³⁴.

### Fix #4 (§3 prose) — one number

Replace

> *"naively combining the two implies ~1.45 A electron current, far beyond
> the 60 mA LHeC design value"*

with

> *"naively combining the two implies ~145 mA electron current — more
> than twice the 60 mA LHeC design value."*

(Same paragraph, qualitative argument unchanged.)

### Fix #5 (§3 prose) — explicit hybrid bunch-charge note

Add one sentence so that the hybrid Tables 2–3 are not silently using a
different N_e from Table 1:

> *"The hybrid adopts a slightly larger electron bunch charge,
> N_e = 1.19 × 10¹⁰ (1.91 nC), than the ESPP-LHeC baseline 1.5 nC, in order
> to keep I_e ≈ 60 mA at the reduced 30 MHz collision frequency."*

This reconciles the two N_e values that currently appear unannotated.

---

## 3. Does Table 2 still deliver 2 × 10³³?

**Yes, with no further edits.** Table 2 is already self-consistent: my
recomputation gives σ\*_e = 36.7/7.5, σ\*_p = 36.8/7.8, σᵧ_eff = 10.85 µm,
σₓ_eff = 51.96 µm, I_e = 57 mA, I_p = 1009 mA, **peak ℒ = 2.117 × 10³³ cm⁻²
s⁻¹** (paper quotes 2.1 × 10³³). The 2 × 10³³ target is met.

The reason it's robust: the hybrid uses N_e = 1.19 × 10¹⁰, n_b = 386,
f_coll = 30 MHz, N_p = 2.1 × 10¹¹, all of which combine to
ℒ ≈ 2.1 × 10³³ at the matched (52 / 11 µm) effective beam sizes. Fix #5
above just makes that choice explicit.

---

## 4. After-fix predictions (for sanity)

| Quantity | Before | After Fix #1 | After Fix #2 | After Fix #3 |
|---|---|---|---|---|
| Table 1 LHeC L | 5.78 × 10³³ (inconsistent with f) | **6.0 × 10³³** | — | — |
| Table 1 EIC I_e | 1555 mA (off by 1.6×) | — | **2500 mA** | — |
| Table 3 ℒ | 1 × 10³⁴ (claimed; recomputes to 9.5 × 10³³) | — | — | **1.0 × 10³⁴** (consistent) |
| Table 2 ℒ | 2.1 × 10³³ ✓ | unchanged | unchanged | unchanged |
| Hybrid I_e | 57.6 mA ✓ | unchanged | unchanged | unchanged |
| Bunch-reduction logic | self-consistent ✓ | unchanged | unchanged | unchanged |

After Fix #1–5, every published number can be re-derived from the input
rows of the same table. The 2 × 10³³ initial-hybrid target and the 1 × 10³⁴
optimised-hybrid target are both achieved.
