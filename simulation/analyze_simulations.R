library(tidyverse)

beta0 <- c(
  0.29346633095464136,
  0.11530484173880995
)

dat <- read_csv("~/TargetedMSM/results/results.csv") %>%
  mutate(index = 1:n())
indices <- dat %>% filter(N == 50) %>%
  group_by(g_correct, Q_correct) %>%
  sample_n(500)

dat <- dat %>%
  filter((N == 50 & (index %in% indices$index)) | (N > 50))
  
dat <- dat %>% 
  group_by(N, g_correct, Q_correct) %>%
  summarize(
    n = n(),
    bias_beta1 = mean(abs(beta1 - beta0[1])),
    bias_beta2 = mean(abs(beta2 - beta0[2])),
    
    bias_beta1_bayes = mean(abs(bayes_beta1 - beta0[1])),
    bias_beta2_bayes = mean(abs(bayes_beta2 - beta0[2])),
    
    mse_beta1 = mean((beta1 - beta0[1])^2),
    mse_beta2 = mean((beta2 - beta0[2])^2),
    
    mse_beta1_bayes = mean((bayes_beta1 - beta0[1])^2),
    mse_beta2_bayes = mean((bayes_beta2 - beta0[2])^2),
    
    ci_beta1 = mean(beta1_upper - beta1_lower),
    ci_beta2 = mean(beta2_upper - beta2_lower),
    
    ci_beta1_bayes = mean(bayes_beta1_upper - bayes_beta1_lower),
    ci_beta2_bayes = mean(bayes_beta2_upper - bayes_beta2_lower),
    
    coverage_beta1 = mean(beta1_lower <= beta0[1] & beta1_upper >= beta0[1]),
    coverage_beta2 = mean(beta2_lower <= beta0[2] & beta2_upper >= beta0[2]),
    
    coverage_beta1 = mean(beta1_lower <= beta0[1] & beta1_upper >= beta0[1]),
    coverage_beta2 = mean(beta2_lower <= beta0[2] & beta2_upper >= beta0[2]),
    
    coverage_beta1_bayes = mean(bayes_beta1_lower <= beta0[1] & bayes_beta1_upper >= beta0[1]),
    coverage_beta2_bayes = mean(bayes_beta2_lower <= beta0[2] & bayes_beta2_upper >= beta0[2])
  ) %>%
  mutate(setting = case_when(
    g_correct == TRUE & Q_correct == TRUE ~ "(a) Q\u305 correct, g correct",
    g_correct == TRUE & Q_correct == FALSE ~ "(b) Q\u305 incorrect, g correct",
    g_correct == FALSE & Q_correct == TRUE ~ "(c) Q\u305 correct, g incorrect",
    g_correct == FALSE & Q_correct == FALSE ~ "(d) Q\u305 incorrect, g incorrect",
  ))

# Table
tab <- dat %>%
  arrange(-g_correct, -Q_correct, N) %>%
  select(g_correct, Q_correct, N, coverage_beta1, coverage_beta1_bayes, bias_beta1, bias_beta1_bayes, coverage_beta2, coverage_beta2_bayes, bias_beta2, bias_beta2_bayes) %>%
  pivot_longer(coverage_beta1:bias_beta2_bayes) %>%
  mutate(estimator = ifelse(str_detect(name, "_bayes"), "Bayesian", "Frequentist"),
         name = str_replace(name, "_bayes", "")) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  mutate_at(vars(contains("bias"), contains("coverage")), signif, 2) %>%
  ungroup() %>%
  select(-g_correct, -Q_correct) %>%
  mutate(N = as.character(N))

tab[seq(2, nrow(tab), 2), 1] <- ""
  
knitr::kable(tab, format = "latex")

# Coverage
cairo_pdf("~/TargetedMSM/coverage.pdf")
dat %>%
  pivot_longer(coverage_beta1:coverage_beta2_bayes) %>%
  mutate(Estimator = ifelse(str_detect(name, "bayes"), "Bayesian", "Frequentist"),
         Parameter = str_replace(name, "coverage_", ""),
         Parameter = str_replace(Parameter, "_bayes", ""),
         Parameter = ifelse(Parameter == "beta1", "\u03B2\u2081", "\u03B2\u2082")) %>%
  ggplot(aes(x = factor(N), y = value, color = Parameter, shape = Estimator)) +
  geom_hline(yintercept = 0.95, lty = 2, color = "black") +
  geom_point(size = 2) +
  geom_line(aes(group = paste0(name, Estimator)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "Empirical 95% Coverage") +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  cowplot::background_grid()
dev.off()

ggsave("~/TargetedMSM/coverage.pdf", width = 7, height = 5)

# Bias
cairo_pdf("~/TargetedMSM/bias.pdf")
dat %>%
  pivot_longer(c(bias_beta1, bias_beta1_bayes, bias_beta2, bias_beta2_bayes)) %>%
  mutate(Estimator = ifelse(str_detect(name, "bayes"), "Bayesian", "Frequentist"),
         Parameter = str_replace(name, "bias_", ""),
         Parameter = str_replace(Parameter, "_bayes", ""),
         Parameter = ifelse(Parameter == "beta1", "\u03B2\u2081", "\u03B2\u2082")) %>%
  ggplot(aes(x = factor(N), y = value, color = Parameter, shape = Estimator)) +
  geom_point() +
  geom_line(aes(group = paste0(name, Estimator)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "Absolute Bias") +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  cowplot::background_grid()
dev.off()

ggsave("~/TargetedMSM/bias.pdf", width = 7, height = 5)

# CI width
dat %>%
  pivot_longer(ci_beta1:ci_beta2_bayes) %>%
  mutate(bayes = str_detect(name, "bayes"),
         name = str_replace(name, "ci_", ""),
         name = str_replace(name, "_bayes", "")) %>%
  ggplot(aes(x = factor(N), y = value, color = name, shape = bayes)) +
  geom_point() +
  geom_line(aes(group = paste0(name, bayes)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "CI")

# MSE
dat %>%
  pivot_longer(mse_beta1:mse_beta2_bayes) %>%
  mutate(bayes = str_detect(name, "bayes"),
         name = str_replace(name, "ci_", ""),
         name = str_replace(name, "_bayes", "")) %>%
  ggplot(aes(x = factor(N), y = value, color = name, shape = bayes)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste0(name, bayes)), alpha = 0.25) +
  facet_wrap(~setting) +
  labs(x = "N", y = "MSE")
ggsave("~/TargetedMSM/mse.pdf", width = 7, height = 5)

