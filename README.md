# Power Calculation R Shiny App

An interactive [Shiny](https://shiny.posit.co/) app for running statistical power
calculations for impact evaluations and research studies. It walks you through
choosing your parameters, shows the trade-off between the minimum detectable
effect size (MDES) and sample size as an interactive table and graph, and
produces a plain-language summary you can paste into a proposal.

Created by **Rigzom Wangchuk**. Based on the
[*"Power Calc – Interactive"* Stata .do file](https://www.fsnnetwork.org/resource/power-calculation-stata-code)
(SCUS_IDEAL_HAEC, Interactive V4, April 2023) by Ananda Young and
Rigzom Wangchuk, written to help applicants to HAEC's research grant run power
calculations in Stata — reimplemented here as a point-and-click web app that
needs no Stata license.

## Features

- **Two-sample means and proportions tests**, for individual-level and
  cluster-randomized designs (cluster designs via the design effect
  1 + (m − 1)ρ, with your choice of fixing cluster number or cluster size)
- **Vary MDES or vary sample size**, in either direction — effects can be
  increases *or* decreases relative to the control group
- **Interactive graph** (plotly): hover for exact values, click a point to
  select it, zoom and pan
- **Scrollable results table** with a sticky header row
- **Bonferroni correction** for multiple primary outcomes
- **Real-world adjustments**: expected attrition (inflates required N by
  1/(1 − a)) and imperfect take-up (by 1/c² under intention-to-treat)
- **Power curve** at your selected MDES (power vs. N)
- **Exports**: graph as a 300-dpi PNG, table as CSV, plus a copyable
  "parameters used" memo for proposal appendices
- **Guardrails**: the do-file's input validation (retry instead of exit), a
  stale-results warning when inputs change after a run, and a one-click reset

## Getting started

Requires R ≥ 4.1 (tested on 4.5). The three packages the app uses — shiny,
ggplot2, and plotly — are installed automatically on first run if missing,
so no separate install step is needed.

Open `Power Calculation R Shiny App.R` in RStudio and click **Source**
(or **Run App**, if shown), or run:

```r
source("Power Calculation R Shiny App.R")
```

The script is self-launching — sourcing it starts the app.

## How to use it

1. **Set standard parameters** — power (default 0.8), alpha (default 0.05),
   and the number of primary outcomes (alpha is Bonferroni-adjusted when > 1).
2. **Choose the test** — means vs. proportions, individual vs. cluster design.
3. **Enter study parameters** — control mean/SD or control proportion, nratio
   or ICC and cluster structure, and the MDES or sample-size range to explore.
   For an expected *decrease*, set the treatment range below the control.
4. **(Optional) real-world adjustments** — expected attrition and take-up.
5. **Run** — review the table and the MDES-vs-N curve.
6. **Pick the point where the curve flattens** (type it, or click the graph),
   using a negative MDES for decreases.
7. **Read the result** — the summary sentence, the attrition/take-up-adjusted
   enrollment target, the parameters memo, and the power curve.

## Methodology

Formulas match Stata's power commands under the normal approximation:

- `power twomeans` — two-sample z-test with known SD, supporting `nratio`
- `power twoproportions` — Pearson's chi-squared (pooled variance), the Stata
  default
- Cluster designs apply the design effect 1 + (m − 1)ρ; cluster size or
  cluster number is solved analytically (means) or by root-finding
  (proportions)

### Distributions and assumptions

- **All critical values and power come from the standard normal (z)
  distribution** (`qnorm()` / `pnorm()`). No t-distribution is used: the
  means test treats the standard deviation as known, matching the default
  behavior of Stata's `power twomeans`. For small samples (roughly under
  ~30 per arm) a t-based calculation would require slightly more
  observations than shown here.
- **Means test** — the outcome is assumed (approximately) normally
  distributed with equal standard deviation in both arms, or the sample
  large enough for the central limit theorem to make the difference in
  means approximately normal.
- **Proportions test** — the outcome is binary (Bernoulli/binomial); the
  test statistic uses the normal approximation to the binomial with the
  pooled variance under the null, i.e. the large-sample version of
  Pearson's chi-squared test. The approximation is weakest when a
  proportion is near 0 or 1 or cell counts are small (a common rule of
  thumb: n·p and n·(1 − p) should both exceed ~5 in each arm).
- **Cluster designs** — outcomes are correlated within clusters; this is
  captured through a single intra-cluster correlation ρ applied as the
  design effect 1 + (m − 1)ρ, assuming a common ρ across arms and equal
  cluster sizes within each arm. The test statistic is still referred to
  the normal distribution (no degrees-of-freedom correction for a small
  number of clusters, which Stata's default likewise omits).

Results may differ from Stata by a unit or two per group due to rounding
conventions. Two-sided tests throughout, so the required N depends only on the
magnitude of the effect, not its direction.

**Note on the original .do file:** in its cluster Method 1 branches the Stata
code passes `k1($treat_cluster_num) k2($treat_cluster_num)`, so the control
cluster count entered by the user is never used. This app uses both K1 and K2
correctly; results will differ from the .do file whenever the two counts
differ.

## Repository structure

```
.
├── Power Calculation R Shiny App.R   # the entire app: stats, UI, and server
├── README.md
├── LICENSE
├── CHANGELOG.md
└── .gitignore
```

## Limitations

- Normal-approximation formulas only (no exact/finite-population corrections,
  no t-distribution small-sample adjustment)
- Cluster calculations assume a common ICC across arms and, where solved,
  a common cluster size or cluster number across arms
- The take-up adjustment uses the standard 1/c² intention-to-treat
  approximation

## Learn more

- The original Stata code: [Power Calculation Stata Code](https://www.fsnnetwork.org/resource/power-calculation-stata-code) (FSN Network)
- J-PAL, [Power calculations](https://www.povertyactionlab.org/resource/power-calculations)
- J-PAL North America, [Quick Guide to Power Calculations](https://www.povertyactionlab.org/resource/quick-guide-power-calculations)
- EGAP, [10 Things to Know About Statistical Power](https://egap.org/resource/10-things-to-know-about-statistical-power)

## License

Licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
(Creative Commons Attribution-NonCommercial): you are free to use, share, and
adapt this app for any **non-commercial** purpose, as long as you credit
Rigzom Wangchuk and link to the license. Commercial use is not permitted.
Provided "as is", without warranty of any kind. See [LICENSE](LICENSE).
