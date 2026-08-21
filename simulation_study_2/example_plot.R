library(tidyverse)
library(patchwork)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_2")

# Load simulation files
source(glue::glue("{basepath}/simulate.R"))

treatments <- 25
n_knots <- 5
data <- simulate_data(seed = 1, N = 500, treatments = treatments, sigma = 0.25, linear = TRUE)
learners_trt <- c("SL.ranger", "SL.glm", "SL.glm.interaction", "SL.earth")
learners_outcome <- c("SL.ranger", "SL.glm", "SL.glm.interaction", "SL.earth")

testpoints <- tibble(a = seq(1, 25, length.out = 100))
f <- ~-1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments))
mat <- model.matrix(f, data = testpoints)
colnames(mat) <- 1:7
splines <- as.data.frame(mat) |>
  mutate(a = testpoints$a) |>
  pivot_longer(cols = `1`:`7`)

fit <- automsm::dose_response(
  data,
  c("x1", "x2", "x3"), "a", "y",
  formula = ~-1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments)),
  outcome_type = "continuous",
  tmle = TRUE,
  bayes = FALSE,
  nuisance = nuisance_control(
    learners_outcome = learners_outcome,
    learners_trt = learners_trt,
    outer_folds = 5,
    epsilon = 1e-3
  )
)

newdata <- tibble(a = 1:25)
newdata$point <- predict(fit, newdata = newdata, estimator = "tmle")
newdata$ui <- predict(fit, newdata = newdata, estimator = "tmle", type = "rvar")
newdata$lower <- map_dbl(newdata$ui, quantile, 0.025)
newdata$upper <- map_dbl(newdata$ui, quantile, 0.975)

p1 <- ggplot(newdata, aes(a, point)) +
  geom_line(aes(color = "B-spline NP-MSM")) +
  geom_line(aes(y = lower, color = "B-spline NP-MSM"), lty = 2) +
  geom_line(aes(y = upper, color = "B-spline NP-MSM"), lty = 2) +
  geom_function(fun = \(x) 2 / x + 1 / (26 - x), aes(color = "True dose-response curve")) +
  scale_color_discrete(name = "") +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  theme(legend.position = "top") +
  labs(x = "", y = "Counterfactual mean", title = "Example B-spline NP-MSM\nfor dose-response curve")

spline_labels <- splines |>
  group_by(name) |>
  filter(value == max(value)) |>
  filter(a == min(a))

p2 <- ggplot(splines, aes(a, value, color = name)) +
  geom_line(lty = 1, alpha = 0.75) +
  geom_text(data = spline_labels, aes(label = name), vjust = 1.75, alpha = 0.75, size = 4) +
  cowplot::theme_cowplot() +
  labs(x = "A", subtitle = "B-spline basis functions") +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p1 / p2 + plot_layout(heights = c(5, 2))

ggsave("simulation_study_2/plots/bspline-msm-example.pdf", width = 7, height = 4, device = cairo_pdf)
