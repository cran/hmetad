## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 150, fig.width = 8, fig.height = 6,
  fig.align = "center", out.width = "75%"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(tidyverse)
library(brms)
library(tidybayes)
library(hmetad)

## ----gumbel_min_r-------------------------------------------------------------
gumbel_min_lcdf <- function(x, g) {
  log1p(-exp(-exp(x - g)))
}
gumbel_min_lccdf <- function(x, g) {
  -exp(x - g)
}

## ----gumbel_min_stan, results=FALSE-------------------------------------------
gumbel_min <- stanvar(
  scode = "
real gumbel_min_lcdf(real x, real g) {
  return log1m_exp(-exp(x - g));
}
real gumbel_min_lccdf(real x, real g) {
  return -exp(x - g);
}",
  block = "functions"
)

## ----data---------------------------------------------------------------------
d <- sim_metad(
  N_trials = 10000, dprime = 1.5, c = .1, log_M = -.5,
  c2_0_diff = c(.25, .5, .25), c2_1_diff = c(.1, .5, .25),
  lcdf = gumbel_min_lcdf, lccdf = gumbel_min_lccdf
)

## ----model, results=FALSE, message=FALSE, warning=FALSE-----------------------
m <- fit_metad(N ~ 1,
  data = d,
  prior = prior(normal(0, 1), class = Intercept) +
    prior(normal(0, 1), class = dprime) +
    prior(normal(0, 1), class = c) +
    prior(lognormal(-1, 1), class = metac2zero1diff) +
    prior(lognormal(-1, 1), class = metac2zero2diff) +
    prior(lognormal(-1, 1), class = metac2zero3diff) +
    prior(lognormal(-1, 1), class = metac2one1diff) +
    prior(lognormal(-1, 1), class = metac2one2diff) +
    prior(lognormal(-1, 1), class = metac2one3diff),
  distribution = "gumbel_min", stanvars = gumbel_min,
)

## ----echo=FALSE---------------------------------------------------------------
summary(m)

## ----roc1---------------------------------------------------------------------
# psuedo type-1 ROC
tibble(.row = 1) |>
  add_roc1_draws(m, bounds = TRUE) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(
    x = p_fa, xmin = p_fa.lower, xmax = p_fa.upper,
    y = p_hit, ymin = p_hit.lower, ymax = p_hit.upper
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(orientation = "y", width = .01) +
  geom_errorbar(orientation = "x", width = .01) +
  geom_point() +
  geom_line() +
  coord_fixed(xlim = 0:1, ylim = 0:1, expand = FALSE) +
  xlab("P(False Alarm)") +
  ylab("P(Hit)") +
  theme_bw(18)

## ----roc2---------------------------------------------------------------------
# type 2 ROC
roc2_draws(m, tibble(.row = 1), bounds = TRUE) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response = factor(response)) |>
  ggplot(aes(
    x = p_fa2, xmin = p_fa2.lower, xmax = p_fa2.upper,
    y = p_hit2, ymin = p_hit2.lower, ymax = p_hit2.upper,
    color = response
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(orientation = "y", width = .01) +
  geom_errorbar(orientation = "x", width = .01) +
  geom_point() +
  geom_line() +
  coord_fixed(xlim = 0:1, ylim = 0:1, expand = FALSE) +
  xlab("P(Type 2 False Alarm)") +
  ylab("P(Type 2 Hit)") +
  theme_bw(18)

