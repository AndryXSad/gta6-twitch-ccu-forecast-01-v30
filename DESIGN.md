# v30 — design record

## What v30 is

v30 starts from the v29 packaged cell and moves it in the direction the
last seven rounds all pointed: the median is the bottleneck, and the
disagreement shift is the only lever that keeps buying practice AND LOO
at the same time. v30 deepens that shift (0.40 + 0.10·size, late band
0.35), steepens the size pull-down (−0.02), the mid ramp (−0.12), the
early ramp (+0.03), the disagreement factor (k = 0.80), and the late
spread ramp (→ 0.90). Nothing about the members changed.

## How the baseline was re-established (and a bug found)

`main.py` was written from the v29 spec before any harness existed. The
first direct grade came out at 0.120935 instead of 0.119155. Comparing
the computed members against the cached v29 rows showed Y5, Z0 and clone
exact to 0.0000 on medians and spreads — but Y1's median was off by
−0.18 (p00) to +0.10 (p02). Testing candidate shape constructions
against the cached rows proved Y1's shape (level, amplitude) regressions
run on **[lnV1, lnPREM]** — the v1 original features — not [lnV1, g12]
as the v29 brief's shorthand ("median template") had suggested. After
that fix, all four members match the cached rows to 0.0000 and the
direct grade reproduces v29 exactly: 0.119155 (per-case 0.1086 / 0.1353
/ 0.0792 / 0.1535), LOO 0.209806 / 0.211492. The combination formula was
also validated separately from the cached rows, reproducing the same
three numbers exactly.

## The search

All sweeps ran through `main.combine` on the cached v29 member rows
(practice + LOO from the same rows), never by duplicating predictor code
in the harness. A cell passes only if it beats v29 on practice AND both
LOO gates.

* **Phase 1 — idea 1 (deeper shift / slopes / ramps / k).** 432 cells
  (shift_base 0.35/0.38/0.40 × shift_slope 0.05/0.075/0.10 × k 0.7/0.75/0.8
  × mid ramp −0.115/−0.12 × early ramp 0.025/0.03 × late shift None/0.35 ×
  late spread end 0.92/0.90). 104 cells passed all three gates. Best
  practice: 0.118929 (0.40/0.10/0.8/−0.12/0.03, no late override).
* **Refinement around the winner.** shift_late 0.30/0.35/0.40 × late
  spread end 0.90/0.92 × size slope −0.015/−0.02. The packaged cell
  (shift_late 0.35, late spread end 0.90, size slope −0.02) scored
  practice **0.118710**, LOO-large **0.209744**, LOO-mega **0.211400** —
  better than the phase-1 best on practice and both LOO gates. Two
  near-tied alternatives: (0.30/0.90/−0.02) 0.118771 / 0.209619 /
  0.211322 and (0.30/0.92/−0.02) 0.118900 / 0.209413 / 0.211139.
  Note the near-tie (0.40/0.90/−0.02) reached practice 0.118651 but
  LOO-mega 0.211493 — one millionth above the v29 gate — and was
  rejected for failing the mega gate.
* **Phase 2 — idea 2 (median sources).** On the best phase-1 cell,
  all 27 combinations of early_src (y1/med3/y1z0) × mid_src
  (y5clone/med4/y5clonez0) × late_src (med3/med4/mean3). Including
  clone in the late median (med4) improved practice to 0.118596 but
  wrecked LOO (0.210970 / 0.212616) — clone is the worst LOO member
  (0.2198) and does not belong in the late median. No median-source
  change passed both LOO gates. Honest conclusion: **there is no
  further median-source lever**; the per-band sources from v24–v29
  stand.
* **Phase 3 — idea 3, wild (spread grows with launch size).**
  `spread *= 1 + slope * max(0, lnV1 - 11.5)`. Slope 0.05 → LOO
  0.209207/0.210809 but practice 0.119236; slope 0.20 → LOO
  0.208199/0.209738 but practice 0.120193. The wild idea trades
  practice for LOO monotonically and fails the practice gate at every
  nonzero slope, so it was **not packaged** (slope 0.0). It does,
  however, expose a real tension: the LOO launches say bigger launches
  need wider spreads, while the practice cases — including p02, the
  biggest and smoothest practice launch (0.0786) — say narrower.
  The final case (lnV1 ~ 14+) is beyond both, so neither gate can
  settle which side is right for the final.

## Results (direct-graded, exact RULES scorer, "Average viewers")

* Practice **0.118710** (p00 0.1083 / p01 0.1355 / p02 0.0786 / p03 0.1525)
  vs v29 0.119155 — **−0.000445**.
* LOO-large **0.209744** vs 0.209806 — **−0.000062**.
* LOO-mega **0.211400** vs 0.211492 — **−0.000092**.
* Median-only practice **0.184501** (p00 0.1675 / p01 0.2277 / p02 0.1050 /
  p03 0.2378) vs v29 0.186892 — **−0.002391**, the seventh round in a row
  the median-only loss moves. Median-only LOO-large 0.317142 (v29
  0.318193), LOO-mega 0.317647.

One direct grade was used (this candidate). The harness sweep numbers
matched the direct grade to six decimals, because the harness scores
through the same `combine`.

## Honest caveats

1. **Selection bias.** Every parameter was chosen on the same four
   practice launches and the same 37-launch LOO that all previous
   rounds selected on. The margins (−0.000445 practice, −0.00006 LOO)
   are thin, and part of them is selection noise. The practice board
   is only a proxy: placements are held until finalize, and the final
   is graded on GTA VI's real release day.
2. **The practice board is flooded.** As of 17 Nov ~21:08 WIB, 198
   entries score below 0.122315 — the top is dominated by replay /
   memorization entries (67 entries at exactly 1e-06). v30 at 0.118710
   would sit around rank #199 on the practice board, nowhere near the
   top 20, let alone the top 10 the owner now requires for submission.
   The final rewards honest generalization (v30's LOO), not memorized
   practice curves — but whether prizes follow the practice rank or the
   final rank is decided by the platform's finalize step, and that is
   not in this repo's control.
3. **Extrapolation.** The final launch (lnV1 ~ 14+) sits beyond the
   library's largest launch (13.68). The size shift (−0.02·size) and
   the size-indexed disagreement shift (0.40 + 0.10·size) both
   extrapolate hardest exactly there, and no gate in this repo can
   validate them at that size.
4. **0.100 is practically unreachable honestly.** The median is the
   bottleneck (best median-only practice ever is now 0.184501), the
   per-hour OLS residual sd is 0.35–0.45 log, and the practice answers
   are 2–2.5× rougher than the comparables. 0.100 would need a near-
   perfect smooth model. v30 gains −0.000445, in line with the last
   rounds' shrinking steps (v28 −0.0009, v29 −0.0005).
5. **Not submitted.** The owner's rule (17 Nov ~21:21 WIB) is to submit
   only a candidate that reaches the top 10. v30 does not. It is a
   verified fallback, better than v27 on all three gates.
