# Changelog

## 1.0.0 — 2026-10-18

First public release.

- Shiny reimplementation of the "Power Calc – Interactive" Stata .do file
  (SCUS_IDEAL_HAEC, Interactive V4): two-sample means and proportions tests,
  individual-level and cluster designs, vary-MDES and vary-sample-size
  methods, with the do-file's prompts, defaults, validation messages, and
  summary wording.
- Interactive plotly graph with hover tooltips, click-to-pick, zoom/pan.
- Scrollable results table with a sticky header row.
- Effects in either direction: treatment ranges may sit below the control,
  and Step 7 accepts negative MDES values.
- Bonferroni alpha adjustment for multiple primary outcomes.
- Real-world adjustments: expected attrition (1/(1 − a)) and take-up (1/c²,
  intention-to-treat) applied to the final result.
- Power-vs-N curve at the selected MDES.
- PNG (graph) and CSV (table) downloads; copyable "parameters used" memo.
- Stale-results warning when inputs change after a run; full reset button.
- Fixes relative to the .do file: input errors allow retry instead of
  exiting; cluster Method 1 uses both the control (K1) and treatment (K2)
  cluster counts (the .do file passes the treatment count for both); the
  MDES-interval check works for decrease ranges; very large effects no
  longer crash the cluster-proportions solver.
