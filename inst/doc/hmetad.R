## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 150, fig.width = 8, fig.height = 6,
  fig.align = "center", out.width = "75%"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(tidyverse)
library(tidybayes)
library(hmetad)

d <- sim_metad(
  N_trials = 1000, dprime = .75, c = -.5, log_M = -1,
  c2_0 = c(.25, .75, 1), c2_1 = c(.5, 1, 1.25)
)

## ----data, echo=FALSE---------------------------------------------------------
d <- d |> select(trial, stimulus, response, confidence)
d

## ----joint_response, echo=FALSE-----------------------------------------------
d.joint_response <- d |>
  ungroup() |>
  mutate(joint_response = joint_response(response, confidence, max(confidence))) |>
  select(trial, joint_response)
d.joint_response

## ----convert_joint_response---------------------------------------------------
d.joint_response |>
  mutate(
    response = type1_response(joint_response, K = 4),
    confidence = type2_response(joint_response, K = 4)
  )

## ----convert_separate_responses-----------------------------------------------
d |>
  mutate(joint_response = joint_response(response, confidence, K = 4))

## ----to_unsigned--------------------------------------------------------------
to_unsigned(c(-1, 1))

## ----to_signed----------------------------------------------------------------
to_signed(c(0, 1))

## ----aggregate_metad----------------------------------------------------------
d.summary <- aggregate_metad(d)

## ----echo=FALSE---------------------------------------------------------------
d.summary

## ----aggregate_metad_response-------------------------------------------------
aggregate_metad(d, .name = "y")

## ----aggregate_condition, eval=FALSE------------------------------------------
# aggregate_metad(d, participant, condition)

## ----model_fitting, results=FALSE, message=FALSE, warning=FALSE---------------
m <- fit_metad(N ~ 1,
  data = d,
  prior = prior(normal(0, 1), class = Intercept) +
    prior(normal(0, 1), class = dprime) +
    prior(normal(0, 1), class = c) +
    prior(lognormal(0, 1), class = metac2zero1diff) +
    prior(lognormal(0, 1), class = metac2zero2diff) +
    prior(lognormal(0, 1), class = metac2one1diff) +
    prior(lognormal(0, 1), class = metac2one2diff)
)

## ----echo=FALSE---------------------------------------------------------------
summary(m)

## ----eval=FALSE---------------------------------------------------------------
# K <- n_distinct(d$confidence)
# 
# m <- brm(bf(...),
#   data = aggregate_metad(d, ...),
#   family = metad(K = K, ...),
#   stanvars = stanvars_metad(K = K, ...),
#   ...
# )

## ----parameters---------------------------------------------------------------
draws.metad <- tibble(.row = 1) |>
  add_linpred_draws_metad(m)

## ----echo=FALSE---------------------------------------------------------------
print(draws.metad)

## -----------------------------------------------------------------------------
draws.metad <- tibble(.row = 1) |>
  add_linpred_draws_metad(m, pivot_longer = TRUE)

## ----echo=FALSE---------------------------------------------------------------
print(draws.metad)

## -----------------------------------------------------------------------------
draws.metad |>
  median_qi()

## ----post_pred1---------------------------------------------------------------
draws.predicted <- predicted_draws_metad(m, d.summary)

## ----echo=FALSE---------------------------------------------------------------
draws.predicted

## -----------------------------------------------------------------------------
draws.predicted |>
  group_by(.row, stimulus, joint_response, response, confidence) |>
  median_qi(.prediction) |>
  group_by(.row) |>
  mutate(N = t(d.summary$N[.row, ])) |>
  ggplot(aes(x = joint_response)) +
  geom_col(aes(y = N), fill = "grey80") +
  geom_pointrange(aes(y = .prediction, ymin = .lower, ymax = .upper)) +
  facet_wrap(~stimulus, labeller = label_both) +
  theme_classic(18)

## ----epred_draws--------------------------------------------------------------
draws.epred <- epred_draws_metad(m, newdata = tibble(.row = 1))

## ----echo=FALSE---------------------------------------------------------------
draws.epred

## ----epred--------------------------------------------------------------------
draws.epred |>
  group_by(.row, stimulus, joint_response, response, confidence) |>
  median_qi(.epred) |>
  group_by(.row) |>
  mutate(.true = t(response_probabilities(d.summary$N[.row, ]))) |>
  ggplot(aes(x = joint_response)) +
  geom_col(aes(y = .true), fill = "grey80") +
  geom_pointrange(aes(y = .epred, ymin = .lower, ymax = .upper)) +
  scale_alpha_discrete(range = c(.25, 1)) +
  facet_wrap(~stimulus, labeller = label_both) +
  theme_classic(18)

## ----mean_confidence----------------------------------------------------------
tibble(.row = 1) |>
  add_mean_confidence_draws(m) |>
  median_qi(.epred) |>
  left_join(d |>
    group_by(stimulus, response) |>
    summarize(.true = mean(confidence)))

## ----mean_confidence2---------------------------------------------------------
tibble(.row = 1) |>
  add_mean_confidence_draws(m, by_stimulus = FALSE) |>
  median_qi(.epred) |>
  left_join(d |>
    group_by(response) |>
    summarize(.true = mean(confidence)))

## ----mean_confidence3---------------------------------------------------------
tibble(.row = 1) |>
  add_mean_confidence_draws(m, by_response = FALSE) |>
  median_qi(.epred) |>
  left_join(d |>
    group_by(stimulus) |>
    summarize(.true = mean(confidence)))

## ----mean_confidence4---------------------------------------------------------
tibble(.row = 1) |>
  add_mean_confidence_draws(m, by_stimulus = FALSE, by_response = FALSE) |>
  median_qi(.epred) |>
  bind_cols(d |>
    ungroup() |>
    summarize(.true = mean(confidence)))

## ----metacognitive_bias_draws-------------------------------------------------
tibble(.row = 1) |>
  add_metacognitive_bias_draws(m) |>
  median_qi()

## ----roc1_draws---------------------------------------------------------------
draws.roc1 <- tibble(.row = 1) |>
  add_roc1_draws(m)

## ----echo=FALSE---------------------------------------------------------------
draws.roc1

## ----roc1---------------------------------------------------------------------
draws.roc1 |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(
    x = p_fa, xmin = p_fa.lower, xmax = p_fa.upper,
    y = p_hit, ymin = p_hit.lower, ymax = p_hit.upper
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(orientation = "y", width = .01) +
  geom_errorbar(orientation = "x", width = .01) +
  geom_line() +
  coord_fixed(xlim = 0:1, ylim = 0:1, expand = FALSE) +
  xlab("P(False Alarm)") +
  ylab("P(Hit)") +
  theme_bw(18)

## ----roc2---------------------------------------------------------------------
draws.roc2 <- tibble(.row = 1) |>
  add_roc2_draws(m)

## ----echo=FALSE---------------------------------------------------------------
draws.roc2

## -----------------------------------------------------------------------------
draws.roc2 |>
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
  geom_line() +
  coord_fixed(xlim = 0:1, ylim = 0:1, expand = FALSE) +
  xlab("P(Type 2 False Alarm)") +
  ylab("P(Type 2 Hit)") +
  theme_bw(18)

