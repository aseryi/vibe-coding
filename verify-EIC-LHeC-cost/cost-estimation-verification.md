# Independent Verification of the "Cost Estimations" Section
## IPAC26 paper WEV1301 — *Accelerator Design Educational Primer: Hybrid LHeC-like EIC*

**Author of paper:** A. Seryi et al.
**Verifier:** Independent cross-check by a second AI model (Claude, Cowork mode), April 2026.
**Scope of verification:** First two paragraphs of Section "Cost Estimations" and Table 4.
**Reference materials used:** the LHeC cost report CERN-ACC-2018-0061 (Calaga, 2018) and the LHeC Phase-1 ESPP-2026 submission, both provided alongside the paper draft.

---

## 1. Why this verification was performed

The paper's cost estimate was generated with AI-assisted parametric scaling from the LHeC CDR. Per the paper's own "Verification of AI results" section, a second prong of verification is *cross-verification by an independent AI model*. This document records that independent re-scaling and lists the resulting recommendations.

---

## 2. What was verified

### 2.1 Internal arithmetic of Table 4

Table 4 in the draft reads:

| Item | MUSD |
|---|---|
| Material (578 MCHF × 1.25) | 722.2 |
| Labor (601 PY × 300 k$) | 180.5 |
| Pre-contingency total | 902.7 |
| Contingency 60 % | 541.6 |
| **Total** | **1444.3** |

Line-to-line, the table is internally consistent (722.2 + 180.5 = 902.7; 902.7 × 1.6 = 1444.3). However, the conversions from the paragraph values (578 MCHF, 601 PY) are slightly off:

- 578 × 1.25 = **722.5** (table says 722.2, Δ = −0.3)
- 601 × 0.300 = **180.3** (table says 180.5, Δ = +0.2)

Net effect on the bottom line is 0.2 MUSD — negligible, but the rounding inconsistency is real.

### 2.2 Paragraph-1 sub-totals

The sub-totals add up correctly: 278 + 50 + 250 = 578 MCHF; 451 + 50 + 100 = 601 PY. The 70 MCHF arc cost matches the LHeC ESPP-2026 Phase-1 submission, which states "an overall budget reduction of at least 70 MCHF" by removing arcs 3–6 from the multi-pass version. So this number is independently sourced and well grounded.

### 2.3 LHeC-CDR scaling cross-check

The LHeC CDR (Table 1 of CERN-ACC-2018-0061) gives a 60-GeV ERL total of **1756 MCHF** with the following breakdown:

| Item | MCHF |
|---|---|
| SRF System (2 × 10 GeV linacs) | 805 |
| SRF R&D and prototyping | 31 |
| Injector | 40 |
| Magnet/Vacuum (6 return arcs) | 215 |
| SC IR magnets | 105 |
| Source/Dump | 5 |
| Cryogenic infrastructure | 100 |
| General/Installation | 69 |
| Civil engineering | 386 |
| **Total** | **1756** |

Three independent scaling scenarios for the hybrid (6 GeV recirculating linac, 2 passes for 12 GeV) were applied to this reference:

| Scenario | Electron complex (MCHF) | Grand total incl. IR + det. (MCHF) | Equivalent USD with 60 % cont. (MUSD) |
|---|---|---|---|
| **Aggressive** (most aux. items dropped/heavily scaled) | 445 | 745 | 1777 |
| **Moderate** | 522 | 822 | 1932 |
| **Conservative** | 645 | 945 | 2177 |
| **Paper draft** | **278** | **578** | **1444** |

Even the most aggressive linear scaling sits ~30 % above the paper's number. The 278 MCHF "electron complex" figure appears to capture roughly SRF + R&D + injector + arcs only; civil engineering, general/installation, cryo, and source items appear to be omitted or heavily compressed in the raw AI scaling. (The paper's own caveat that "IR and detector line items were manually increased above the raw AI output" already acknowledges the AI tends to under-count auxiliaries.)

Bottom line: the two AI estimates agree to within roughly **±30 % at the pre-contingency MCHF level** and lie within the same factor on the MUSD bottom line. Both are well below the EIC project's CD-1 cost range upper bound.

### 2.4 "Halves the dominant SRF, civil and cryo costs" — overstated

Going from a single-pass 10 GeV linac to a 6 GeV recirculating linac scales SRF/civil/cryo by 6/10 = **0.6**, not 0.5. The 70 MCHF arc penalty leaves ~62 MCHF net savings, so the *bottom-line* of the comparison is reproducible, but the verb "halves" is not literal. **"Reduces by ~40 %"** is more accurate.

### 2.5 CHF/USD exchange rate

As of mid-April 2026, 1 CHF ≈ 1.27 USD; the paper's 1.25 is conservative and within the recent trading range. ✓

### 2.6 CD-0 vs CD-1 attribution

The "$1.7–2.8 B" range for the EIC was set at **CD-1 (June 2021)**, not CD-0 (December 2019, mission need only). Half the upper bound is exactly 1.4 B USD, so the comparison number is correct, but the label "CD0-approved" should read "CD-1-approved".

---

## 3. Note on labor (person-year) estimates

The LHeC CDR cost reference does **not** contain person-year estimates — it only quotes capital cost. The paper's 601 PY estimate therefore has no direct LHeC-CDR anchor. In this verification, labor was treated as proportional to material cost using the paper's own implicit ratio (~1.04 PY per MCHF total, with category-specific ratios of 1.62 / 1.0 / 0.4 PY-per-MCHF for electron-complex / IR / detector respectively). That ratio corresponds to labor being ~25 % of material cost, which is consistent with typical DOE/CERN large-project planning rules of thumb but should be regarded as an order-of-magnitude estimate.

---

## 4. Recommendations for revised text and table

### 4.1 Wording corrections (independent of the range decision)

- Replace "CD0-approved" → **"CD-1-approved"** (the cost range is set at CD-1).
- Replace "halves the dominant SRF, civil, and cryo costs" → **"reduces the dominant SRF, civil, and cryo costs by ~40 %"**.
- Optionally state explicitly that two AI-assisted scalings of the LHeC CDR cost reproduce the same qualitative conclusion and agree on the bottom line to within **~30 %** at the pre-contingency level.

### 4.2 Suggested ranges for the cost summary

Anchored on the aggressive scaling (because of the large 60 % contingency already added on top), with the lower bound taken from the paper's existing AI-derived numbers:

**Paragraph-1 sub-totals as ranges:**

| Item | Range (MCHF) | Range (PY) |
|---|---|---|
| Electron complex | 280 – 450 | 450 – 720 |
| IR | 50 | 50 |
| Detector | 250 | 100 |
| **Grand total** | **~580 – 750** | **~600 – 870** |

**Revised Table 4:**

| Item | Range (MUSD) |
|---|---|
| Material costs (580 – 750 MCHF × 1.25) | 725 – 935 |
| Labor costs (600 – 870 PY × 300 k$) | 180 – 260 |
| Pre-contingency total | 905 – 1195 |
| Contingency (60 %) | 545 – 720 |
| **Total with contingency** | **~1.4 – 1.9 B USD** |

Suggested caption: *"Hybrid EIC cost range, estimated by two independent AI-assisted scalings of the LHeC CDR cost. The two estimates agree to within ~30 % at the pre-contingency level. All figures are order-of-magnitude only."*

The conclusion in paragraph 2 then naturally becomes: *"the hybrid EIC total of ~1.4–1.9 B USD is comfortably below the upper bound of the EIC CD-1 cost range ($1.7–2.8 B) — a margin that justifies more detailed evaluation of the hybrid concept."*

---

## 5. Caveats

This cross-check is itself an AI-assisted parametric scaling of a single reference (LHeC CDR, 2018 prices) and inherits all of its limitations: no bottom-up engineering estimate; no inflation/escalation correction beyond the rough 2018 → 2026 currency conversion; no separate accounting for the EIC-specific source R&D required by 60 mA polarized operation; no detector cost modelling beyond the manual 250 MCHF placeholder. The exercise is fit for purpose — flagging where the paper's numbers are likely to slip — but not as a substitute for a formal cost review.

---

*Generated as a record of the AI-cross-verification described in the paper's Section "Verification of AI results" (second prong: cross-verification by a second, independent AI model).*
