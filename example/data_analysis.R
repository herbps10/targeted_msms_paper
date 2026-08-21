library(tidyverse)
library(torch)
library(automsm)

dat <- read_csv("example/data/pnas_data.csv")

X <- c(
  "curr_use_inj_adj_bl",
  "area_bl",
  "age_group_bl",
  "sex_age_bl",
  "total_alive_bl",
  "ever_use_bl",
  "edu_primary_bl",
  "work_bl",
  "religion_r_bl",
  "ethnicity_r_bl"
)
A <- "treatment"
Y <- "curr_use_inj_adj"

learners_trt <- c("SL.glm", "SL.glm.interaction", "SL.ranger", "SL.mean")
learners_outcome <- c("SL.glm", "SL.glm.interaction", "SL.ranger", "SL.mean")

set.seed(10016)
fit <- automsm::cate(
  dat[, c(X, A, Y)], X, A, Y,
  formula = ~1 + total_alive_bl,
  outcome_type = "binomial",
  tmle = TRUE,
  bayes = bayes_control(
    chains = 4,
    warmup = 5e2,
    draws = 4e3,
    prior = \(beta) sum(dnorm(as.numeric(beta), mean = 0, sd = 5, log = TRUE)),
    seed = 10001
  ),
  nuisance = nuisance_control(
    learners_trt = learners_trt,
    learners_outcome = learners_outcome
  ),
)

fit$bayes_tmle$diagnostics
tidy(fit)

range1 <- range(post$`beta[1]`)
range2 <- range(post$`beta[2]`)

post |>
  pivot_longer(everything()) |>
  filter(name %in% c("beta[1]", "beta[2]")) |>
  mutate(name = ifelse(name == "beta[1]", "\u03B2\u2081", "\u03B2\u2082")) |>
  ggplot(aes(x = value)) +
  geom_line(aes(color = "TMLE", x = x, y = y), data = tibble(name = "\u03B2\u2081", x = seq(range1[1], range1[2], length.out = 100), y = dnorm(x, mean = fit$tmle$est[1], sd = fit$tmle$se[1]))) +
  geom_line(aes(color = "TMLE", x = x, y = y), data = tibble(name = "\u03B2\u2082", x = seq(range2[1], range2[2], length.out = 100), y = dnorm(x, mean = fit$tmle$est[2], sd = fit$tmle$se[2]))) +
  geom_density(aes(color = "Bayes TMLE")) +
  guides(color = "none") +
  facet_wrap(~name, scales = "free") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  cowplot::panel_border() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
ggsave("example/plots/pnas_posterior.pdf", width = 7, height = 4, device = cairo_pdf)

preds <- tibble(total_alive_bl = 0:9)
preds$post <- predict(fit, newdata = data.frame(total_alive_bl = 0:9), estimator = "bayes", type = "rvar")

ui <- preds |>
  unnest("post") |>
  group_by(total_alive_bl) |>
  ggdist::median_qi(post, .width = c(0.5, 0.8, 0.95))

ggplot(ui, aes(x = total_alive_bl, y = post)) +
  ggdist::geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  geom_line() +
  scale_fill_brewer(name = "") +
  scale_x_continuous(breaks = c(0:9)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  labs(x = "Number of children at baseline", y = "Intervention effect size", title = "Intervention effect size increases with number of children")
ggsave("example/plots/pnas_posterior_cate2.pdf", width = 7, height = 5, device = cairo_pdf)
