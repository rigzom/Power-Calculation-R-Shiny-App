# ============================================================================
# Power Calculation R Shiny App — interactive power calculations for studies
#
# Created by:   Rigzom Wangchuk
# Date:         October 18, 2026
# Based on:     "Power Calc - Interactive" Stata .do file (SCUS_IDEAL_HAEC,
#               Interactive V4, April 2023) by Ananda Young & Rigzom Wangchuk
#               https://www.fsnnetwork.org/resource/power-calculation-stata-code
# License:      CC BY-NC 4.0 — free to use, share, and adapt for
#               non-commercial purposes with attribution; no warranty.
#               https://creativecommons.org/licenses/by-nc/4.0/
#
# Mirrors the do-file's flow: power & alpha defaults -> test type -> design ->
# method -> parameter entry -> table + MDES-vs-N graph -> pick a point ->
# summary sentence in the do-file's wording.
#
# Extras beyond the do-file:
#   - hover tooltips + click-to-pick on the graph (plotly)
#   - PNG / CSV downloads of the graph and table
#   - Bonferroni alpha adjustment for multiple primary outcomes
#   - attrition and take-up adjustments in the final result
#   - a power-vs-N curve at the selected MDES
#   - copyable "parameters used" memo, stale-results warning, full reset
#
# Formulas match Stata's `power twomeans` (z-test, known sd, nratio) and
# `power twoproportions` (Pearson chi-squared, pooled). Cluster designs use
# the design effect 1 + (m - 1) * rho. Rounding conventions may differ from
# Stata by a unit or two per group.
#
# To run:  source this file (it is self-launching). Required packages are
#          installed automatically on first run if missing.
# ============================================================================

# Install required packages if not already installed, then load them
packages <- c("shiny", "ggplot2", "plotly")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(packages, library, character.only = TRUE))

# Small helper: fall back to b when a is NULL/NA (input not rendered yet)
orv <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

# Two-sample means, individual level: required N given delta
means_n <- function(m1, m2, sd, alpha, power, r) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  delta <- m2 - m1
  n1 <- ceiling((za + zb)^2 * sd^2 * (1 + 1 / r) / delta^2 - 1e-9)
  n2 <- ceiling(r * n1 - 1e-9)
  list(N = n1 + n2, n1 = n1, n2 = n2, delta = delta)
}

# Two-sample means, individual level: detectable delta given total N
means_delta <- function(n, sd, alpha, power, r) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  n1 <- n / (1 + r); n2 <- n - n1
  (za + zb) * sd * sqrt(1 / n1 + 1 / n2)
}

# Cluster means, Method 1: K1, K2 fixed -> solve common cluster size M
means_cluster_m <- function(m1, m2, sd, alpha, power, rho, K1, K2) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  delta <- m2 - m1
  C <- (za + zb)^2 * sd^2 * (1 / K1 + 1 / K2)
  X <- delta^2 / C
  if (X <= rho) return(NULL)  # infeasible: ICC puts a floor on the SE
  M <- ceiling((1 - rho) / (X - rho) - 1e-9)
  list(M = M, N1 = K1 * M, N2 = K2 * M, N = (K1 + K2) * M, delta = delta)
}

# Cluster means, Method 2: M1, M2 fixed -> solve common cluster number K
means_cluster_k <- function(m1, m2, sd, alpha, power, rho, M1, M2) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  delta <- m2 - m1
  D1 <- 1 + (M1 - 1) * rho; D2 <- 1 + (M2 - 1) * rho
  K <- ceiling((za + zb)^2 * sd^2 * (D1 / M1 + D2 / M2) / delta^2 - 1e-9)
  list(K = K, N1 = K * M1, N2 = K * M2, N = K * (M1 + M2), delta = delta)
}

# Two-sample proportions, individual: required N given p1, p2 (pooled chi2)
props_n <- function(p1, p2, alpha, power, r) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  d <- abs(p2 - p1)
  pbar <- (p1 + r * p2) / (1 + r); qbar <- 1 - pbar
  num <- za * sqrt((1 + 1 / r) * pbar * qbar) +
         zb * sqrt(p1 * (1 - p1) + p2 * (1 - p2) / r)
  n1 <- ceiling(num^2 / d^2 - 1e-9)
  n2 <- ceiling(r * n1 - 1e-9)
  list(N = n1 + n2, n1 = n1, n2 = n2, delta = p2 - p1)
}

# Two-sample proportions, individual: detectable delta given total N
props_delta <- function(n, p1, alpha, power, r) {
  za <- qnorm(1 - alpha / 2)
  n1 <- n / (1 + r)
  pw <- function(p2) {
    d <- p2 - p1
    pbar <- (p1 + r * p2) / (1 + r); qbar <- 1 - pbar
    se0 <- sqrt((1 + 1 / r) * pbar * qbar)
    se1 <- sqrt(p1 * (1 - p1) + p2 * (1 - p2) / r)
    pnorm((d * sqrt(n1) - za * se0) / se1)
  }
  upper <- 1 - 1e-8
  if (pw(upper) < power) return(NA_real_)
  root <- uniroot(function(p2) pw(p2) - power,
                  lower = p1 + 1e-8, upper = upper, tol = 1e-10)
  root$root - p1
}

# Cluster proportions: power of the design-effect-adjusted pooled z-test
props_cluster_power <- function(p1, p2, alpha, K1, K2, M1, M2, rho) {
  za <- qnorm(1 - alpha / 2)
  d <- abs(p2 - p1)
  D1 <- 1 + (M1 - 1) * rho; D2 <- 1 + (M2 - 1) * rho
  n1 <- K1 * M1; n2 <- K2 * M2
  pbar <- (n1 * p1 + n2 * p2) / (n1 + n2); qbar <- 1 - pbar
  se0 <- sqrt(pbar * qbar * (D1 / n1 + D2 / n2))
  se1 <- sqrt(D1 * p1 * (1 - p1) / n1 + D2 * p2 * (1 - p2) / n2)
  pnorm((d - za * se0) / se1)
}

# Cluster proportions, Method 1: K1, K2 fixed -> solve common cluster size M
props_cluster_m <- function(p1, p2, alpha, power, rho, K1, K2) {
  f <- function(M) props_cluster_power(p1, p2, alpha, K1, K2, M, M, rho) - power
  if (f(1e7) < 0) return(NULL)   # target power unreachable (ICC floor)
  # If even cluster size 1 already reaches the target power (very large
  # effects), there is no root to find -> the answer is simply M = 1.
  M <- if (f(1) >= 0) 1 else
    ceiling(uniroot(f, lower = 1, upper = 1e7, tol = 1e-8)$root - 1e-9)
  list(M = M, N1 = K1 * M, N2 = K2 * M, N = (K1 + K2) * M, delta = p2 - p1)
}

# Cluster proportions, Method 2: M1, M2 fixed -> solve common cluster number K
props_cluster_k <- function(p1, p2, alpha, power, rho, M1, M2) {
  f <- function(K) props_cluster_power(p1, p2, alpha, K, K, M1, M2, rho) - power
  if (f(1e7) < 0) return(NULL)
  # Same guard: one cluster per arm may already be enough for huge effects.
  K <- if (f(1) >= 0) 1 else
    ceiling(uniroot(f, lower = 1, upper = 1e7, tol = 1e-8)$root - 1e-9)
  list(K = K, N1 = K * M1, N2 = K * M2, N = K * (M1 + M2), delta = p2 - p1)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

stata_css <- "
  body { background-color: #EEF1EE; }
  .stata-red   { color: #8B1E2D; font-family: monospace; }
  .stata-green { color: #1F6E43; font-family: monospace; font-size: 12.5px; }
  .well { background-color: #FFFFFF; border: 1px solid #D4DBD4; }
  .summary-box { background: #22302B; color: #F2F5F1; border-radius: 8px;
                 padding: 18px 20px; font-size: 15px; line-height: 1.6; }
  .summary-adjust { margin: 10px 0 0; color: #9FD8B4; font-size: 14px; }
  .flag-box { background: #FFF7E6; border: 1px solid #E8D9A8; border-radius: 6px;
              padding: 10px 14px; color: #6B5A20; font-size: 12.5px; margin-top: 12px; }
  h2 { font-weight: 800; letter-spacing: -0.02em; }
  /* Scrollable results table with a sticky header row */
  .table-scroll { max-height: 300px; overflow-y: auto; border: 1px solid #D4DBD4;
                  border-radius: 6px; background: #FFFFFF; margin-bottom: 12px; }
  .table-scroll table { margin-bottom: 0; }
  .table-scroll thead th { position: sticky; top: 0; z-index: 2; background: #22302B;
                           color: #F2F5F1; border-bottom: none;
                           box-shadow: inset 0 -2px 0 #8B1E2D; }
  .kicker { font-family: monospace; font-size: 11px; letter-spacing: .14em;
            text-transform: uppercase; color: #5C6B63; }
  /* Input labels: main topic stays bold, the (range / suggestion) hint does not */
  .label-hint { font-weight: 400; color: #5C6B63; }
  /* Copyable parameters memo under the result */
  .memo { font-family: monospace; font-size: 12px; background: #F6F8F5;
          border: 1px solid #D4DBD4; border-radius: 6px; padding: 12px 14px;
          margin-top: 10px; white-space: pre-wrap; color: #22302B; }
  .dl-row { margin: 0 0 18px; }
  .dl-row .btn { margin-right: 8px; }
"

# Input-label helper: everything before the first "(" is the main topic and
# keeps Shiny's default bold; the parenthetical range/suggestion (and trailing
# colon) renders in normal weight via the .label-hint CSS class.
lbl <- function(txt) {
  i <- regexpr("(", txt, fixed = TRUE)
  if (i < 0) return(txt)
  HTML(paste0(htmltools::htmlEscape(trimws(substr(txt, 1, i - 1))),
              " <span class=\"label-hint\">",
              htmltools::htmlEscape(substr(txt, i, nchar(txt))),
              "</span>"))
}

ui <- fluidPage(
  title = "Power Calculation R Shiny App",
  tags$head(tags$style(HTML(stata_css))),
  br(),
  div(class = "kicker", "Power Calculations · Interactive"),
  h2("Power Calculation R Shiny App"),
  p(class = "stata-green",
    "Statistical power is the probability of correctly rejecting the null hypothesis when it is false. "),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      # Step 1: standard parameters -------------------------------------
      h4("Step 1: Standard parameters"),
      numericInput("power", lbl("Power (between 0 and 1; typically 0.8 or 0.9):"), value = 0.8,
                   min = 0.01, max = 0.99, step = 0.05),
      numericInput("alpha", lbl("Alpha (between 0 and 1; typically 0.01, 0.05, or 0.10):"),
                   value = 0.05, min = 0.001, max = 0.5, step = 0.01),
      numericInput("n_outcomes",
                   lbl("Number of primary outcomes (whole number ≥ 1; alpha is divided by this — Bonferroni):"),
                   value = 1, min = 1, step = 1),
      p(class = "stata-green",
        "With several primary outcomes, the calculations use alpha / number of outcomes ",
        "(Bonferroni correction) to keep the overall false-positive risk at your chosen alpha."),

      # Step 2: test type -------------------------------------------------
      h4("Step 2: Test type"),
      radioButtons("test", "Two-sample means or proportions test?",
                   choices = c("Means" = "means", "Proportions" = "props")),
      radioButtons("design", "Individual-level or cluster design?",
                   choices = c("Individual-level" = "ind", "Cluster" = "cl")),

      # Step 3-4: branch-specific parameters -------------------------------
      h4("Step 3: Parameters"),
      uiOutput("params"),

      # Step 5: real-world adjustments --------------------------------------
      h4("Step 5: Real-world adjustments (optional)"),
      p(class = "stata-green",
        "Attrition inflates the required sample by 1/(1−a). ",
        "Imperfect take-up inflates it by 1/c² under an intention-to-treat analysis. ",
        "These adjustments are applied to the final Step 8 result only."),
      numericInput("attrition",
                   lbl("Expected attrition rate (0–0.99; share of the enrolled sample lost by endline; 0 = none):"),
                   value = 0, min = 0, max = 0.99, step = 0.05),
      numericInput("takeup",
                   lbl("Expected take-up rate (0–1; share of the treatment group actually treated; 1 = full):"),
                   value = 1, min = 0.01, max = 1, step = 0.05),

      actionButton("run", "Run power calculations", class = "btn-primary"),
      br(), br(),
      actionButton("reset_all", "Clear and start again"),
      uiOutput("bug_flag")
    ),

    mainPanel(
      width = 8,
      uiOutput("stale_banner"),
      uiOutput("results_header"),
      plotlyOutput("power_plot", height = "340px"),
      div(class = "table-scroll", tableOutput("power_table")),
      uiOutput("dl_ui"),
      uiOutput("pick_ui"),
      uiOutput("summary_ui"),
      uiOutput("power_curve_block")
    )
  ),
  br()
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  is_means   <- reactive(input$test == "means")
  is_cluster <- reactive(input$design == "cl")

  # Bonferroni-adjusted alpha used in every calculation
  alpha_eff <- reactive({
    k <- max(1, round(orv(input$n_outcomes, 1)))
    input$alpha / k
  })

  # --- Full reset ----------------------------------------------------------
  observeEvent(input$reset_all, session$reload())

  # --- Dynamic parameter inputs, mirroring the do-file's prompts ----------
  output$params <- renderUI({
    common_means <- tagList(
      numericInput("control_mean", lbl("Control mean (usually 100, or baseline mean if known):"), 100),
      numericInput("sd", lbl("Standard deviation (> 0; usually 100, or baseline SD if known):"), 100, min = 1e-9)
    )
    mde_range_means <- tagList(
      p(class = "stata-green",
        "The difference between the treatment and control means is the assumed MDES. ",
        "Expecting an increase? Set the range above the control mean (e.g. 100 to 130). ",
        "Expecting a decrease? Set it below (e.g. 70 to 100)."),
      numericInput("tm_lower", lbl("Lower bound of treatment mean (for a decrease, set below the control mean, e.g. 70):"), 100),
      numericInput("tm_upper", lbl("Upper bound of treatment mean (must exceed the lower bound; for an increase, above the control, e.g. 130):"), 130),
      numericInput("tm_int", lbl("MDES interval to graph over (smaller than upper − lower bound; e.g. 2–10):"), 5, min = 1e-9)
    )
    mde_range_props <- tagList(
      p(class = "stata-green",
        "Intervals are in percentage points: 0.01 = 1 pp, 0.1 = 10 pp. ",
        "Expecting an increase? Set the range above the control proportion (e.g. 0.5 to 0.65). ",
        "Expecting a decrease? Set it below (e.g. 0.35 to 0.5)."),
      numericInput("tp_lower", lbl("Lower bound of treatment proportion (between 0 and 1; for a decrease, below the control, e.g. 0.35):"),
                   0.5, min = 0.001, max = 0.999, step = 0.01),
      numericInput("tp_upper", lbl("Upper bound of treatment proportion (between 0 and 1; must exceed the lower bound; for an increase, above the control, e.g. 0.65):"), 0.65,
                   min = 0.001, max = 0.999, step = 0.01),
      numericInput("tp_int", lbl("MDES interval to graph over (between 0 and 1; smaller than upper − lower; typically 0.01–0.05):"),
                   0.03, min = 0.001, max = 0.999, step = 0.01)
    )
    nratio_input <- tagList(
      p(class = "stata-green",
        "nratio = treatment size / control size. 1 = equal groups; control ×2 of treatment = 0.5; ",
        "control ½ of treatment = 2."),
      numericInput("nratio", lbl("nratio (> 0; typically between 0.5 and 2):"), 1, min = 1e-9, step = 0.1)
    )
    icc_input <- tagList(
      p(class = "stata-green",
        "The ICC compares within-cluster to between-cluster variance. Usually 0.2 if no baseline data. ",
        "With baseline data: loneway outcome cluster_id, then r(rho)."),
      numericInput("icc", lbl("ICC, ρ (between 0 and 1; typically 0.05–0.20 for field studies):"),
                   0.2, min = 0.001, max = 0.999, step = 0.01)
    )
    n_range <- function(up_default, up_hint, int_default, int_hint) tagList(
      numericInput("n_lower", lbl("Lower bound of sample size (positive whole number; e.g. 100):"), 100, min = 1),
      numericInput("n_upper", lbl(paste0("Upper bound of sample size (must exceed the lower bound; e.g. ", up_hint, "):")), up_default, min = 2),
      numericInput("n_int", lbl(paste0("Sample size interval to graph over (smaller than upper − lower; e.g. ", int_hint, "):")), int_default, min = 1)
    )
    cluster_counts <- tagList(
      numericInput("K1", lbl("Number of control clusters (positive whole number; typically 10–60 per arm):"), 20, min = 1),
      numericInput("K2", lbl("Number of treatment clusters (positive whole number; typically 10–60 per arm):"), 20, min = 1))
    cluster_sizes <- tagList(
      numericInput("M1", lbl("Control cluster size (positive whole number; typically 10–50 units per cluster):"), 25, min = 1),
      numericInput("M2", lbl("Treatment cluster size (positive whole number; typically 10–50 units per cluster):"), 25, min = 1))

    if (is_means() && !is_cluster()) {
      tagList(
        nratio_input, common_means,
        radioButtons("method", "Step 4: vary MDES or sample size?",
                     choices = c("1 — Vary MDES (get sample sizes)" = "mde",
                                 "2 — Vary sample size (get MDEs)" = "n")),
        conditionalPanel("input.method == 'mde'", mde_range_means),
        conditionalPanel("input.method == 'n'", n_range(1000, "500–2,000", 100, "50–200"))
      )
    } else if (is_means() && is_cluster()) {
      tagList(
        common_means, mde_range_means, icc_input,
        radioButtons("method", "Step 4: what varies alongside MDES?",
                     choices = c("1 — Cluster size varies (cluster number fixed)" = "mde",
                                 "2 — Cluster number varies (cluster size fixed)" = "n")),
        conditionalPanel("input.method == 'mde'", cluster_counts),
        conditionalPanel("input.method == 'n'", cluster_sizes)
      )
    } else if (!is_means() && !is_cluster()) {
      tagList(
        nratio_input,
        numericInput("control_prop", lbl("Control proportion (0–1; usually 0.5, or baseline if known):"),
                     0.5, min = 0.001, max = 0.999, step = 0.01),
        radioButtons("method", "Step 4: vary MDES or sample size?",
                     choices = c("1 — Vary MDES (get sample sizes)" = "mde",
                                 "2 — Vary sample size (get MDEs)" = "n")),
        conditionalPanel("input.method == 'mde'", mde_range_props),
        conditionalPanel("input.method == 'n'", n_range(2000, "1,000–5,000", 200, "100–500"))
      )
    } else {
      tagList(
        numericInput("control_prop", lbl("Control proportion (0–1; usually 0.5, or baseline if known):"),
                     0.5, min = 0.001, max = 0.999, step = 0.01),
        mde_range_props, icc_input,
        radioButtons("method", "Step 4: what varies alongside MDES?",
                     choices = c("1 — Cluster size varies (cluster number fixed)" = "mde",
                                 "2 — Cluster number varies (cluster size fixed)" = "n")),
        conditionalPanel("input.method == 'mde'", cluster_counts),
        conditionalPanel("input.method == 'n'", cluster_sizes)
      )
    }
  })

  output$bug_flag <- renderUI({
    req(is_cluster(), input$method == "mde")
    div(class = "flag-box",
        HTML("<b>Note on the original .do file:</b> in cluster Method 1 the Stata code passes
              <code>k1($treat_cluster_num) k2($treat_cluster_num)</code>, so the control cluster
              count is never used. This app uses K1 and K2 correctly."))
  })

  # --- Stale-results guard -------------------------------------------------
  has_run <- reactiveVal(FALSE)
  stale   <- reactiveVal(FALSE)
  observeEvent(input$run, { has_run(TRUE); stale(FALSE) }, priority = 10)
  observe({
    # Touch every input that feeds results(); any change after a run marks
    # the displayed table/graph as out of date. (Attrition/take-up are not
    # listed: they only affect Step 8, which recomputes live.)
    list(input$power, input$alpha, input$n_outcomes, input$test, input$design,
         input$method, input$nratio, input$control_mean, input$sd,
         input$tm_lower, input$tm_upper, input$tm_int,
         input$n_lower, input$n_upper, input$n_int,
         input$control_prop, input$tp_lower, input$tp_upper, input$tp_int,
         input$icc, input$K1, input$K2, input$M1, input$M2)
    if (isolate(has_run())) stale(TRUE)
  })
  output$stale_banner <- renderUI({
    req(has_run(), stale())
    div(class = "flag-box", style = "margin: 0 0 12px;",
        strong("Inputs changed since the last run. "),
        "The table, graph, and summary below still reflect the previous inputs — ",
        "click ‘Run power calculations’ to refresh.")
  })

  # --- Validation (do-file error codes + the new optional inputs) ----------
  validate_inputs <- function() {
    validate(
      need(input$power > 0 && input$power < 1,
           "ERROR: Power must be between 0 and 1. Start again!"),
      need(input$alpha > 0 && input$alpha < 1,
           "ERROR: Alpha must be between 0 and 1. Start again!"),
      need(!is.na(input$n_outcomes) && input$n_outcomes >= 1 &&
             abs(input$n_outcomes - round(input$n_outcomes)) < 1e-9,
           "ERROR: Number of primary outcomes must be a whole number of at least 1."),
      need(!is.na(input$attrition) && input$attrition >= 0 && input$attrition < 1,
           "ERROR: Attrition must be at least 0 and below 1. Start again!"),
      need(!is.na(input$takeup) && input$takeup > 0 && input$takeup <= 1,
           "ERROR: Take-up must be above 0 and at most 1. Start again!")
    )
    if (!is_cluster()) {
      validate(need(input$nratio > 0,
                    "ERROR: nratio must be a positive number. Start again!"))
    }
    if (is_means()) {
      validate(need(input$sd > 0, "ERROR: Standard deviation must be positive."))
      if (is_cluster() || input$method == "mde") {
        validate(
          need(input$tm_upper > input$tm_lower,
               "ERROR: Upper bound of the treatment mean must exceed the lower bound. Start again!"),
          need(input$tm_int > 0 && input$tm_int < (input$tm_upper - input$tm_lower),
               "ERROR: MDES interval is equal to or more than the width of the treatment-mean range. Start again!")
        )
      }
    } else {
      validate(need(input$control_prop > 0 && input$control_prop < 1,
                    "ERROR: Control proportion must be between 0 and 1. Start again!"))
      if (is_cluster() || input$method == "mde") {
        validate(
          need(input$tp_lower > 0 && input$tp_lower < 1 && input$tp_upper > 0 && input$tp_upper < 1,
               "ERROR: Treatment proportions must be between 0 and 1. Start again!"),
          need(input$tp_upper > input$tp_lower,
               "ERROR: Upper bound of the treatment proportion must exceed the lower bound. Start again!"),
          need(input$tp_int > 0 && input$tp_int < (input$tp_upper - input$tp_lower),
               "ERROR: MDES interval is equal to or more than the width of the treatment-proportion range. Start again!")
        )
      }
    }
    if (!is_cluster() && input$method == "n") {
      validate(
        need(input$n_lower > 0, "ERROR: Lower bound of sample size must be more than 0. Start again!"),
        need(input$n_upper > input$n_lower,
             "ERROR: Upper bound of sample size must exceed the lower bound. Start again!"),
        need(input$n_int > 0 && input$n_int < abs(input$n_upper - input$n_lower),
             "ERROR: The sample size interval is equal to or more than the difference between the lower and upper limits of the sample size. Start again!")
      )
    }
    if (is_cluster()) {
      validate(need(input$icc > 0 && input$icc < 1,
                    "ERROR: ICC must be between 0 and 1. Start again!"))
      if (input$method == "mde") {
        validate(need(input$K1 > 0 && input$K2 > 0, "ERROR: Cluster numbers must be positive."))
      } else {
        validate(need(input$M1 > 0 && input$M2 > 0, "ERROR: Cluster sizes must be positive."))
      }
    }
  }

  # Build a grid of values, refusing absurdly fine grids (they make the app slow)
  make_grid <- function(lo, hi, by) {
    g <- seq(lo, hi, by = by)
    validate(need(length(g) <= 300,
      sprintf("ERROR: this grid has %d points, which will be slow to compute and plot. Increase the interval so the grid has at most 300 points.", length(g))))
    g
  }

  # --- Results grid (the do-file's "table") --------------------------------
  results <- eventReactive(input$run, {
    validate_inputs()
    a <- alpha_eff(); p <- input$power

    if (is_means() && !is_cluster() && input$method == "mde") {
      grid <- make_grid(input$tm_lower, input$tm_upper, input$tm_int)
      rows <- lapply(grid, function(tm) {
        if (isTRUE(all.equal(tm, input$control_mean))) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = 0))
        r <- means_n(input$control_mean, tm, input$sd, a, p, input$nratio)
        data.frame(N = r$N, N1 = r$n1, N2 = r$n2, delta = r$delta)
      })
      df <- do.call(rbind, rows); df$extra <- NA; extra_lab <- NULL
    } else if (is_means() && !is_cluster() && input$method == "n") {
      grid <- make_grid(input$n_lower, input$n_upper, input$n_int)
      df <- data.frame(
        N = grid,
        N1 = floor(grid / (1 + input$nratio)),
        N2 = grid - floor(grid / (1 + input$nratio)),
        delta = sapply(grid, means_delta, sd = input$sd, alpha = a, power = p, r = input$nratio),
        extra = NA
      ); extra_lab <- NULL
    } else if (is_means() && is_cluster()) {
      grid <- make_grid(input$tm_lower, input$tm_upper, input$tm_int)
      rows <- lapply(grid, function(tm) {
        if (isTRUE(all.equal(tm, input$control_mean))) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = 0, extra = NA))
        r <- if (input$method == "mde")
          means_cluster_m(input$control_mean, tm, input$sd, a, p, input$icc, input$K1, input$K2)
        else
          means_cluster_k(input$control_mean, tm, input$sd, a, p, input$icc, input$M1, input$M2)
        if (is.null(r)) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = tm - input$control_mean, extra = NA))
        data.frame(N = r$N, N1 = r$N1, N2 = r$N2, delta = r$delta,
                   extra = if (input$method == "mde") r$M else r$K)
      })
      df <- do.call(rbind, rows)
      extra_lab <- if (input$method == "mde") "Cluster size (M)" else "Clusters per arm (K)"
    } else if (!is_means() && !is_cluster() && input$method == "mde") {
      grid <- make_grid(input$tp_lower, input$tp_upper, input$tp_int)
      rows <- lapply(grid, function(tp) {
        if (isTRUE(all.equal(tp, input$control_prop))) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = 0))
        r <- props_n(input$control_prop, tp, a, p, input$nratio)
        data.frame(N = r$N, N1 = r$n1, N2 = r$n2, delta = r$delta)
      })
      df <- do.call(rbind, rows); df$extra <- NA; extra_lab <- NULL
    } else if (!is_means() && !is_cluster() && input$method == "n") {
      grid <- make_grid(input$n_lower, input$n_upper, input$n_int)
      df <- data.frame(
        N = grid,
        N1 = floor(grid / (1 + input$nratio)),
        N2 = grid - floor(grid / (1 + input$nratio)),
        delta = sapply(grid, props_delta, p1 = input$control_prop, alpha = a, power = p, r = input$nratio),
        extra = NA
      ); extra_lab <- NULL
    } else {
      grid <- make_grid(input$tp_lower, input$tp_upper, input$tp_int)
      rows <- lapply(grid, function(tp) {
        if (isTRUE(all.equal(tp, input$control_prop))) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = 0, extra = NA))
        r <- if (input$method == "mde")
          props_cluster_m(input$control_prop, tp, a, p, input$icc, input$K1, input$K2)
        else
          props_cluster_k(input$control_prop, tp, a, p, input$icc, input$M1, input$M2)
        if (is.null(r)) return(data.frame(N = NA, N1 = NA, N2 = NA, delta = tp - input$control_prop, extra = NA))
        data.frame(N = r$N, N1 = r$N1, N2 = r$N2, delta = r$delta,
                   extra = if (input$method == "mde") r$M else r$K)
      })
      df <- do.call(rbind, rows)
      extra_lab <- if (input$method == "mde") "Cluster size (M)" else "Clusters per arm (K)"
    }

    df$alpha <- a; df$power <- p
    list(df = df, extra_lab = extra_lab,
         test = input$test, design = input$design, method = input$method,
         vary_n = (!is_cluster() && input$method == "n"))
  })

  output$results_header <- renderUI({
    req(results())
    tagList(
      h4("Step 6: Table & graph"),
      p(class = "stata-green",
        "See where the trade-off between MDES and sample size is, and what value is acceptable — ",
        "pick the point after which the curve flattens.",
        if (results()$test == "props")
          " Note: the MDES (y-axis) is in proportions — 0.01 = 1 percentage point." else "")
    )
  })

  # --- Graph (shared between on-screen plotly and the PNG download) --------
  build_gg <- reactive({
    r <- results()
    d <- r$df[!is.na(r$df$N) & !is.na(r$df$delta) & r$df$delta != 0, ]
    validate(need(nrow(d) > 0, "No feasible points to plot — adjust the parameters (e.g., more clusters or a larger MDES)."))
    dig <- if (r$test == "props") 3 else 2
    d$hover <- sprintf("Total sample size (N): %s<br>MDES: %s",
                       formatC(d$N, format = "d", big.mark = ","),
                       formatC(d$delta, format = "f", digits = dig))
    g <- ggplot(d, aes(x = N, y = delta, group = 1, text = hover)) +
      geom_line(color = "#8B1E2D", linewidth = 1.1) +
      geom_point(color = "#8B1E2D", size = 2.3) +
      labs(x = "Total sample size (N)",
           y = if (r$test == "props") "MDES (proportion)" else "MDES (delta)",
           title = "Trade-off between MDES and sample size") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
    sel <- pick_result()
    if (!is.null(sel) && !is.na(sel$N)) {
      seldf <- data.frame(
        N = sel$N, delta = sel$delta,
        hover = sprintf("Your selection<br>Total sample size (N): %s<br>MDES: %s",
                        formatC(sel$N, format = "d", big.mark = ","),
                        formatC(sel$delta, format = "f", digits = dig)))
      g <- g + geom_point(data = seldf, aes(x = N, y = delta, text = hover),
                          color = "#1F6E43", size = 4.5, inherit.aes = FALSE)
    }
    g
  })

  output$power_plot <- renderPlotly({
    p <- ggplotly(build_gg(), tooltip = "text", source = "powerplot")
    p <- event_register(p, "plotly_click")
    config(p, displayModeBar = FALSE)
  })

  # Click-to-pick: clicking a point on the graph fills the Step 7 box
  observeEvent(event_data("plotly_click", source = "powerplot"), {
    ed <- event_data("plotly_click", source = "powerplot")
    req(ed)
    r <- results()
    val <- if (r$vary_n) round(ed$x[1])
           else round(ed$y[1], if (r$test == "props") 3 else 2)
    updateNumericInput(session, "pick", value = val)
  }, ignoreInit = TRUE)

  output$power_table <- renderTable({
    r <- results()
    d <- r$df
    out <- data.frame(
      alpha = formatC(d$alpha, format = "f", digits = 3),
      power = formatC(d$power, format = "f", digits = 2),
      N  = ifelse(is.na(d$N), "not feasible", formatC(d$N, format = "d", big.mark = ",")),
      `N1 (control)`   = ifelse(is.na(d$N1), "", formatC(d$N1, format = "d", big.mark = ",")),
      `N2 (treatment)` = ifelse(is.na(d$N2), "", formatC(d$N2, format = "d", big.mark = ",")),
      delta = ifelse(d$delta == 0, "0 (no effect)",
                     formatC(d$delta, format = "f", digits = if (r$test == "props") 3 else 2)),
      check.names = FALSE
    )
    if (!is.null(r$extra_lab)) out[[r$extra_lab]] <-
      ifelse(is.na(d$extra), "", formatC(d$extra, format = "d", big.mark = ","))
    out
  }, striped = TRUE, align = "r")

  # --- Downloads -------------------------------------------------------------
  output$dl_ui <- renderUI({
    req(results())
    div(class = "dl-row",
        downloadButton("dl_plot", "Download graph (PNG)"),
        downloadButton("dl_table", "Download table (CSV)"))
  })
  output$dl_table <- downloadHandler(
    filename = function() paste0("power_table_", Sys.Date(), ".csv"),
    content = function(file) {
      r <- results(); d <- r$df
      out <- data.frame(alpha = d$alpha, power = d$power, N_total = d$N,
                        N_control = d$N1, N_treatment = d$N2, MDES = d$delta)
      if (!is.null(r$extra_lab)) out[[make.names(r$extra_lab)]] <- d$extra
      write.csv(out, file, row.names = FALSE, na = "")
    })
  output$dl_plot <- downloadHandler(
    filename = function() paste0("power_graph_", Sys.Date(), ".png"),
    content = function(file) suppressWarnings(
      ggsave(file, plot = build_gg(), width = 8, height = 5, dpi = 300, bg = "white")))

  # --- Step 7: pick from the graph -----------------------------------------
  output$pick_ui <- renderUI({
    r <- results()
    tagList(
      h4("Step 7: Select from the graph"),
      if (r$vary_n)
        p(class = "stata-red", "Enter sample size (on x-axis) where efficiency gain is maximized, i.e. after which the curve flattens:")
      else
        tagList(
          p(class = "stata-red", "Enter MDES (on y-axis) where efficiency gain is maximized, i.e. after which the curve flattens:"),
          p(class = "stata-green",
            "Use a negative MDES if the treatment value is below the control — e.g. −10 for means, −0.10 for proportions.")
        ),
      p(class = "stata-green", "Tip: click a point on the graph to fill this box, then press Enter."),
      fluidRow(
        column(4, numericInput("pick", NULL, value = NA,
                               step = if (r$test == "props" && !r$vary_n) 0.01 else 1)),
        column(3, actionButton("go_pick", "Enter"))
      )
    )
  })

  compute_pick <- function() {
    r <- results(); a <- alpha_eff(); p <- input$power
    # For proportion MDES picks (positive or negative), the implied treatment
    # proportion must stay strictly between 0 and 1.
    if (r$test == "props" && !r$vary_n) {
      p2 <- input$control_prop + input$pick
      if (p2 <= 0 || p2 >= 1) return(list(N = NA, delta = NA,
        text = "ERROR: control proportion + MDES must lie strictly between 0 and 1. Enter a smaller MDES."))
    }
    if (r$test == "means" && r$design == "ind" && r$method == "mde") {
      res <- means_n(input$control_mean, input$control_mean + input$pick, input$sd, a, p, input$nratio)
      list(N = res$N, n1 = res$n1, n2 = res$n2, delta = res$delta,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s with a probability of %s if the effect is true and the ratio of units in treatment and control is %s.",
                          format(res$N, big.mark = ","), format(res$n2, big.mark = ","), format(res$n1, big.mark = ","),
                          round(res$delta, 3), p, input$nratio))
    } else if (r$test == "means" && r$design == "ind" && r$method == "n") {
      d <- means_delta(input$pick, input$sd, a, p, input$nratio)
      n1 <- floor(input$pick / (1 + input$nratio)); n2 <- input$pick - n1
      list(N = input$pick, n1 = n1, n2 = n2, delta = d,
           text = sprintf("A sample size of %s (treatment = %s and control = %s) is able to detect an effect size of %s with a probability of %s if the effect is true and the ratio of units in treatment and control is %s.",
                          format(input$pick, big.mark = ","), format(n2, big.mark = ","), format(n1, big.mark = ","),
                          round(d, 2), p, input$nratio))
    } else if (r$test == "means" && r$design == "cl" && r$method == "mde") {
      res <- means_cluster_m(input$control_mean, input$control_mean + input$pick, input$sd, a, p, input$icc, input$K1, input$K2)
      if (is.null(res)) return(list(N = NA, delta = NA,
        text = "This MDES is not detectable with the given number of clusters and ICC — increase the number of clusters or the MDES."))
      list(N = res$N, n1 = res$N1, n2 = res$N2, delta = res$delta, M = res$M,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s with a probability of %s if the effect is true. There are %s control clusters each of size %s. There are %s treatment clusters each of size %s.",
                          format(res$N, big.mark = ","), format(res$N2, big.mark = ","), format(res$N1, big.mark = ","),
                          round(res$delta, 3), p, input$K1, res$M, input$K2, res$M))
    } else if (r$test == "means" && r$design == "cl" && r$method == "n") {
      res <- means_cluster_k(input$control_mean, input$control_mean + input$pick, input$sd, a, p, input$icc, input$M1, input$M2)
      list(N = res$N, n1 = res$N1, n2 = res$N2, delta = res$delta, K = res$K,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s with a probability of %s if the effect is true. There are %s control clusters each of size %s. There are %s treatment clusters each of size %s.",
                          format(res$N, big.mark = ","), format(res$N2, big.mark = ","), format(res$N1, big.mark = ","),
                          round(res$delta, 3), p, res$K, input$M1, res$K, input$M2))
    } else if (r$test == "props" && r$design == "ind" && r$method == "mde") {
      res <- props_n(input$control_prop, input$control_prop + input$pick, a, p, input$nratio)
      list(N = res$N, n1 = res$n1, n2 = res$n2, delta = res$delta,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s percentage points with a probability of %s if the effect is true and the ratio of units in treatment and control is %s.",
                          format(res$N, big.mark = ","), format(res$n2, big.mark = ","), format(res$n1, big.mark = ","),
                          round(res$delta * 100, 1), p, input$nratio))
    } else if (r$test == "props" && r$design == "ind" && r$method == "n") {
      d <- props_delta(input$pick, input$control_prop, a, p, input$nratio)
      if (is.na(d)) return(list(N = NA, delta = NA,
        text = "No detectable effect below 100% at this sample size — increase the sample size."))
      n1 <- floor(input$pick / (1 + input$nratio)); n2 <- input$pick - n1
      list(N = input$pick, n1 = n1, n2 = n2, delta = d,
           text = sprintf("A sample size of %s (treatment = %s and control = %s) is able to detect an effect size of %s percentage points with a probability of %s if the effect is true and the ratio of units in treatment and control is %s.",
                          format(input$pick, big.mark = ","), format(n2, big.mark = ","), format(n1, big.mark = ","),
                          round(d * 100, 2), p, input$nratio))
    } else if (r$test == "props" && r$design == "cl" && r$method == "mde") {
      res <- props_cluster_m(input$control_prop, input$control_prop + input$pick, a, p, input$icc, input$K1, input$K2)
      if (is.null(res)) return(list(N = NA, delta = NA,
        text = "This MDES is not detectable with the given number of clusters and ICC — increase the number of clusters or the MDES."))
      list(N = res$N, n1 = res$N1, n2 = res$N2, delta = res$delta, M = res$M,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s percentage points with a probability of %s if the effect is true. There are %s control clusters each of size %s. There are %s treatment clusters each of size %s.",
                          format(res$N, big.mark = ","), format(res$N2, big.mark = ","), format(res$N1, big.mark = ","),
                          round(res$delta * 100, 1), p, input$K1, res$M, input$K2, res$M))
    } else {
      res <- props_cluster_k(input$control_prop, input$control_prop + input$pick, a, p, input$icc, input$M1, input$M2)
      if (is.null(res)) return(list(N = NA, delta = NA,
        text = "This MDES is not detectable with the given cluster sizes and ICC."))
      list(N = res$N, n1 = res$N1, n2 = res$N2, delta = res$delta, K = res$K,
           text = sprintf("The minimum sample size needed is %s (treatment = %s and control = %s) to detect an effect size of %s percentage points with a probability of %s if the effect is true. There are %s control clusters each of size %s. There are %s treatment clusters each of size %s.",
                          format(res$N, big.mark = ","), format(res$N2, big.mark = ","), format(res$N1, big.mark = ","),
                          round(res$delta * 100, 1), p, res$K, input$M1, res$K, input$M2))
    }
  }

  # Holds the Step 7 selection; NULL until the user clicks Enter.
  pick_result <- reactiveVal(NULL)
  observeEvent(input$run, { pick_result(NULL) })
  observeEvent(input$go_pick, {
    req(!is.na(input$pick))
    if (stale()) {
      pick_result(list(N = NA, delta = NA,
        text = "Inputs changed since the last run — click ‘Run power calculations’ before selecting from the graph."))
      return()
    }
    # A sample-size pick must be positive; an MDES pick may be negative
    # (treatment below control) but not zero.
    if (results()$vary_n) req(input$pick > 0) else req(input$pick != 0)
    pick_result(compute_pick())
  })

  # --- Real-world adjustment of the final result ----------------------------
  adjust_text <- function(sel, r) {
    attr <- orv(input$attrition, 0); tk <- orv(input$takeup, 1)
    if ((attr <= 0 && tk >= 1) || is.na(sel$N)) return(NULL)
    bits <- c(if (attr > 0) paste0(round(attr * 100, 1), "% expected attrition"),
              if (tk < 1)  paste0("take-up of ", round(tk * 100, 1), "%"))
    lead <- paste0("Real-world adjustment (", paste(bits, collapse = " and "), "): ")
    if (r$vary_n) {
      n_eff <- max(4, floor(sel$N * (1 - attr)))
      d0 <- if (r$test == "means")
        means_delta(n_eff, input$sd, alpha_eff(), input$power, input$nratio)
      else
        props_delta(n_eff, input$control_prop, alpha_eff(), input$power, input$nratio)
      if (is.na(d0)) return(paste0(lead, "the effective sample is too small to detect any effect — increase the sample size."))
      d_eff <- d0 / tk
      paste0(lead, "the smallest effect this design can reliably detect grows to about ",
             if (r$test == "props") paste0(round(d_eff * 100, 1), " percentage points")
             else round(d_eff, 2),
             " (intention-to-treat approximation).")
    } else {
      fac <- 1 / ((1 - attr) * tk^2)
      if (r$design == "cl" && r$method == "mde") {
        M_adj <- ceiling(sel$M * fac)
        N1a <- input$K1 * M_adj; N2a <- input$K2 * M_adj
        paste0(lead, "plan to enroll ", format(N1a + N2a, big.mark = ","),
               " in total (treatment = ", format(N2a, big.mark = ","),
               ", control = ", format(N1a, big.mark = ","),
               "), i.e. about ", M_adj, " units per cluster.")
      } else if (r$design == "cl") {
        K_adj <- ceiling(sel$K * fac)
        N1a <- K_adj * input$M1; N2a <- K_adj * input$M2
        paste0(lead, "plan to enroll ", format(N1a + N2a, big.mark = ","),
               " in total (treatment = ", format(N2a, big.mark = ","),
               ", control = ", format(N1a, big.mark = ","),
               "), i.e. about ", K_adj, " clusters per arm.")
      } else {
        n1a <- ceiling(sel$n1 * fac); n2a <- ceiling(sel$n2 * fac)
        paste0(lead, "plan to enroll ", format(n1a + n2a, big.mark = ","),
               " in total (treatment = ", format(n2a, big.mark = ","),
               ", control = ", format(n1a, big.mark = ","), ").")
      }
    }
  }

  # --- Copyable parameters memo ----------------------------------------------
  memo_text <- function(sel, r, adj) {
    k <- max(1, round(orv(input$n_outcomes, 1)))
    attr <- orv(input$attrition, 0); tk <- orv(input$takeup, 1)
    method_lab <- if (r$design == "ind") {
      if (r$method == "mde") "vary MDES (solve for sample size)" else "vary sample size (solve for MDE)"
    } else {
      if (r$method == "mde") "vary MDES and cluster size (cluster number fixed)"
      else "vary MDES and cluster number (cluster size fixed)"
    }
    lines <- c(
      "PARAMETERS USED",
      paste0("Generated: ", format(Sys.Date(), "%d %b %Y")),
      paste0("Test: two-sample ", if (r$test == "means") "means" else "proportions",
             "  |  Design: ", if (r$design == "ind") "individual-level" else "cluster"),
      paste0("Method: ", method_lab),
      paste0("Power target: ", input$power, "  |  Alpha: ", input$alpha,
             if (k > 1) paste0(" (Bonferroni-adjusted to ", signif(alpha_eff(), 3),
                               " for ", k, " primary outcomes)") else ""),
      if (r$design == "ind") paste0("nratio (treatment/control): ", input$nratio),
      if (r$test == "means") paste0("Control mean: ", input$control_mean, "  |  SD: ", input$sd)
        else paste0("Control proportion: ", input$control_prop),
      if (r$design == "cl") paste0("ICC: ", input$icc,
        if (r$method == "mde") paste0("  |  Clusters: ", input$K1, " control, ", input$K2, " treatment")
        else paste0("  |  Cluster sizes: ", input$M1, " control, ", input$M2, " treatment")),
      paste0("Selected ", if (r$vary_n) "sample size: " else "MDES: ", input$pick),
      paste0("Attrition assumed: ", round(attr * 100, 1), "%  |  Take-up assumed: ", round(tk * 100, 1), "%"),
      "",
      "RESULT",
      sel$text,
      if (!is.null(adj)) adj
    )
    paste(lines, collapse = "\n")
  }

  output$summary_ui <- renderUI({
    sel <- pick_result()
    req(sel)
    r <- results()
    adj <- adjust_text(sel, r)
    tagList(
      h4("Step 8: Result"),
      div(class = "summary-box",
          sel$text,
          if (!is.null(adj)) tags$p(class = "summary-adjust", adj)),
      if (!is.na(sel$N)) tags$pre(class = "memo", memo_text(sel, r, adj)),
      br()
    )
  })

  # --- Power curve at the selected MDES --------------------------------------
  achieved_power <- function(Ntot, sel, r) {
    a <- alpha_eff(); za <- qnorm(1 - a / 2)
    if (r$test == "means") {
      d <- abs(sel$delta); sd <- input$sd
      if (r$design == "ind") {
        n1 <- Ntot / (1 + input$nratio); n2 <- Ntot - n1
        pnorm(d / (sd * sqrt(1 / n1 + 1 / n2)) - za)
      } else if (r$method == "mde") {
        M <- Ntot / (input$K1 + input$K2); D <- 1 + (M - 1) * input$icc
        pnorm(d / (sd * sqrt(D / (input$K1 * M) + D / (input$K2 * M))) - za)
      } else {
        K <- Ntot / (input$M1 + input$M2)
        D1 <- 1 + (input$M1 - 1) * input$icc; D2 <- 1 + (input$M2 - 1) * input$icc
        pnorm(d / (sd * sqrt(D1 / (K * input$M1) + D2 / (K * input$M2))) - za)
      }
    } else {
      p1 <- input$control_prop
      p2 <- min(max(p1 + sel$delta, 1e-6), 1 - 1e-6)
      if (r$design == "ind") {
        rr <- input$nratio; n1 <- Ntot / (1 + rr)
        dd <- abs(p2 - p1)
        pbar <- (p1 + rr * p2) / (1 + rr); qbar <- 1 - pbar
        se0 <- sqrt((1 + 1 / rr) * pbar * qbar)
        se1 <- sqrt(p1 * (1 - p1) + p2 * (1 - p2) / rr)
        pnorm((dd * sqrt(n1) - za * se0) / se1)
      } else if (r$method == "mde") {
        M <- Ntot / (input$K1 + input$K2)
        props_cluster_power(p1, p2, a, input$K1, input$K2, M, M, input$icc)
      } else {
        K <- Ntot / (input$M1 + input$M2)
        props_cluster_power(p1, p2, a, K, K, input$M1, input$M2, input$icc)
      }
    }
  }

  output$power_curve_block <- renderUI({
    sel <- pick_result()
    req(sel, !is.na(sel$N), !is.na(sel$delta), sel$delta != 0)
    tagList(
      h4("Power curve at your selected MDES"),
      p(class = "stata-green",
        "Achieved power as the sample grows, holding the MDES fixed (before any attrition/take-up adjustment). ",
        "Dashed lines mark your target power and the required sample size."),
      plotlyOutput("power_curve", height = "300px"),
      br()
    )
  })

  output$power_curve <- renderPlotly({
    sel <- pick_result()
    req(sel, !is.na(sel$N), !is.na(sel$delta), sel$delta != 0)
    r <- results()
    Ns <- seq(max(4, ceiling(0.15 * sel$N)), ceiling(2 * sel$N), length.out = 60)
    pw <- vapply(Ns, achieved_power, numeric(1), sel = sel, r = r)
    dfp <- data.frame(
      N = Ns, power = pw,
      hover = sprintf("N: %s<br>Power: %s%%",
                      formatC(round(Ns), format = "d", big.mark = ","),
                      formatC(pw * 100, format = "f", digits = 1)))
    g <- ggplot(dfp, aes(x = N, y = power, group = 1, text = hover)) +
      geom_line(color = "#2C5F8A", linewidth = 1.1) +
      geom_hline(yintercept = input$power, linetype = "dashed", color = "#5C6B63") +
      geom_vline(xintercept = sel$N, linetype = "dashed", color = "#1F6E43") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Total sample size (N)", y = "Achieved power") +
      theme_minimal(base_size = 13)
    ggplotly(g, tooltip = "text") |> config(displayModeBar = FALSE)
  })
}

# Self-launching: sourcing this file (or clicking Source / Run App in RStudio)
# starts the app, whatever the file is named.
shiny::runApp(shinyApp(ui, server))
