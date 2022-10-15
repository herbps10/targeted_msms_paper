library(tidyverse)
library(tidybayes)
library(ggdist)

frequentist <- read_csv("~/TargetedMSM/results/pnas_frequentist.csv")

frequentist %>% mutate_at(vars(beta, lower, upper), signif, 2)

z <- frequentist$beta / (frequentist$se / sqrt(frequentist$n))
2 * pnorm(-abs(z))

posterior <- read_csv("~/TargetedMSM/results/pnas_posterior.csv")

range1 <- range(posterior$beta1)
range2 <- range(posterior$beta2)

cairo_pdf("~/TargetedMSM/plots/pnas_posterior.pdf", width = 9, height = 5)
posterior %>%
  pivot_longer(everything()) %>%
  filter(name %in% c("beta1", "beta2")) %>%
  mutate(name = ifelse(name == "beta1", "\u03B2\u2081", "\u03B2\u2082")) %>%
  ggplot(aes(x = value)) +
  geom_line(aes(color = "Frequentist", x = x, y = y), data = tibble(name = "\u03B2\u2081", x = seq(range1[1], range1[2], length.out = 100), y = dnorm(x, mean = frequentist$beta[[1]], sd = frequentist$se[[1]] / sqrt(frequentist$n[[1]])))) +
  geom_line(aes(color = "Frequentist", x = x, y = y), data = tibble(name = "\u03B2\u2082", x = seq(range2[1], range2[2], length.out = 100), y = dnorm(x, mean = frequentist$beta[[2]], sd = frequentist$se[[2]] / sqrt(frequentist$n[[1]])))) +
  geom_density(aes(color = "Bayesian")) +
  guides(color = "none") +
  facet_wrap(~name, scales = "free") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  cowplot::panel_border()
dev.off()


cairo_pdf("~/TargetedMSM/plots/pnas_posterior_cate.pdf", width = 9, height = 7)
tibble(
  children = 0:9
) %>%
  full_join(posterior, by = character()) %>%
  mutate(y = beta1 + beta2 * children) %>%
  group_by(children) %>%
  ggplot(aes(x = factor(children), y = y)) +
  #geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
  stat_eye() +
  #geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0) +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  labs(x = "Number of children at baseline", y = "Intervention effect size", title = "Intervention effect size increases with number of children")
dev.off()

cairo_pdf("~/TargetedMSM/plots/pnas_posterior_example.pdf", width = 7, height = 5)
tibble(
  children = 0:9,
  y = -0.01 + 0.03 * children
) %>%
  ggplot(aes(x = children, y = y)) +
  geom_line() +
  scale_x_continuous(expand = c(0, 0), breaks = 0:9) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  annotate("text", 6, 0.15, label = expression(beta[1]==-0.01), size = 7, hjust = 0) +
  annotate("text", 6, 0.13, label = expression(beta[2]==0.03),  size = 7, hjust = 0) +
  labs(x = "Number of children at baseline", y = "Intervention effect size")
dev.off()

cairo_pdf("~/TargetedMSM/plots/pnas_posterior_cate2.pdf", width = 9, height = 7)
tibble(
  children = 0:9
) %>%
  full_join(posterior, by = character()) %>%
  mutate(y = beta1 + beta2 * children) %>%
  group_by(children) %>%
  median_qi(y, .width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(x = children, y = y)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
  #geom_abline(intercept = frequentist$beta[1], slope = frequentist$beta[2], lty = 2, color = "red") +
  scale_fill_brewer() +
  scale_x_continuous(expand = c(0, 0), breaks = 0:9) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  labs(x = "Number of children at baseline", y = "Intervention effect size", title = "Intervention effect size increases with number of children")
dev.off()

tibble(
  children = 0:9
) %>%
  full_join(sample_n(mutate(posterior, index = 1:n()), 50), by = character()) %>%
  mutate(y = beta1 + beta2 * children) %>%
  ggplot(aes(x = children, y = y, group = index)) +
  geom_line() +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  labs(x = "Number of children at baseline", y = "Estimated CATE")
