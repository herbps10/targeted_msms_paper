library(tidyverse)

results_path <- "simulation_study_1/results"

results_summarized <- read_rds(glue::glue("{results_path}/results_summarized.rds")) |>
  filter(!is.na(scenario)) |>
  mutate(
    setting = case_when(
      scenario == 1 ~ "(a) μ correct, π correct",
      scenario == 2 ~ "(b) μ incorrect, π correct",
      scenario == 3 ~ "(c) μ correct, π incorrect",
      scenario == 4 ~ "(d) μ incorrect, π incorrect"
    ),
    Estimator = case_when(
      estimator == "bayes" ~ "Bayes TMLE",
      estimator == "tmle" ~ "TMLE",
      estimator == "onestep" ~ "One-step"
    ),
    Parameter = ifelse(term == "(Intercept)", "\u03B2\u2081", "\u03B2\u2082")
  )

# Coverage
results_summarized |>
  select(setting, Estimator, Parameter, N, coverage) |>
  ggplot(aes(x = factor(N), y = coverage, color = Parameter, shape = Estimator)) +
  geom_hline(yintercept = 0.95, lty = 2, color = "black") +
  geom_point(size = 2) +
  geom_line(aes(group = paste0(Parameter, Estimator)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "Empirical 95% Coverage") +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  cowplot::background_grid()
ggsave("~/targeted_msms_paper/simulation_study_1/plots/simulation_1_coverage.pdf", width = 7, height = 5, device = cairo_pdf)

# MAE
results_summarized |>
  select(setting, Estimator, Parameter, N, mae) |>
  ggplot(aes(x = factor(N), y = mae, color = Parameter, shape = Estimator)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste0(Parameter, Estimator)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "Mean Absolute Error") +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  cowplot::background_grid()
ggsave("~/targeted_msms_paper/simulation_study_1/plots/simulation_1_mae.pdf", width = 7, height = 5, device = cairo_pdf)
