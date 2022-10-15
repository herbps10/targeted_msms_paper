library(tidyverse)
library(tmle)
library(superleaner)
library(haven)

covariates <- c("curr_use_inj_adj_bl",
  "area_bl",
  "age_group_bl",
  "sex_age_bl",
  "total_alive_bl",
  "ever_use_bl",
  "edu_primary_bl",
  "work_bl",
  "religion_r_bl",
  "ethnicity_r_bl")

dat <- read_dta("~/TargetedMSM/pnas_data/MFPS_PNAS_MAIN_ITT_DATA.dta") %>%
  mutate_at(vars(treatment, curr_use_inj_adj_bl, area_bl, age_group_bl, ever_use_bl, edu_primary_bl, work_bl, religion_r_bl, ethnicity_r_bl), as.factor) %>%
  filter(year == 2018) %>%
  select(curr_use_inj_adj, treatment, !!covariates) %>%
  drop_na()

cairo_pdf("~/TargetedMSM/plots/pnas_effect_modifier.pdf", width = 9, height = 4)
count(dat, total_alive_bl) %>%
  ggplot(aes(x = factor(total_alive_bl), y = n)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +
  geom_text(aes(label = n), vjust = 0, nudge_y = 10) +
  cowplot::theme_cowplot() +
  cowplot::background_grid(major = "y") +
  cowplot::panel_border() +
  labs(x = "Number of children at baseline")
dev.off()
 
mean(as.numeric(dat$treatment) - 1)

# Replicate original analysis
# Unadjusted
lm(curr_use_inj_adj ~ 1 + treatment, data = filter(dat, year == 2018))

# Adjusted
f <- as.formula(glue::glue("curr_use_inj_adj ~ 1 + treatment + {str_c(covariates, collapse = ' + ')}"))
lm(f, data = dat)

SL.library <- c(
  "SL.glm", 
  "SL.xgboost", 
  "SL.glm.interaction",
  "SL.rpart", 
  "SL.mean",
  "SL.ranger", 
  "SL.glmnet",
  "SL.hal9001"
)

fit <- tmle(
  Y = as.numeric(dat$curr_use_inj_adj),
  A = as.numeric(dat$treatment) - 1,
  W = as.data.frame(dat[, covariates]),
  Q.SL.library = SL.library,
  g.SL.library = SL.library,
  family = "binomial",
  verbose = TRUE
)

results <- tibble(
  age = 13:28,
  fit = map(age, function(age) lm(curr_use_inj_adj ~ 1 + treatment, data = filter(dat, sex_age_bl == age))),
  coef = map_dbl(fit, function(fit) coef(fit)[[2]])
)
ggplot(results, aes(x = age, y = coef)) +
  geom_point()

results <- tibble(
  total_alive = 1:6,
  fit = map(total_alive, function(total_alive) lm(curr_use_inj_adj ~ 1 + treatment, data = filter(dat, total_alive_bl == total_alive))),
  coef = map_dbl(fit, function(fit) coef(fit)[[2]])
)
ggplot(results, aes(x = total_alive, y = coef)) +
  geom_point()

plot(dat$sex_age_bl, fit$Qstar[, "Q1W"] - fit$Qstar[, "Q0W"])

augmented_data <- dat %>%
  mutate(Qbar0 = fit$Qinit$Q[, "Q0W"],
         Qbar1 = fit$Qinit$Q[, "Q1W"],
         Qbar  = ifelse(treatment == 1, Qbar1, Qbar0),
         g     = fit$g$g1W,
         treatment = as.numeric(treatment) - 1)

write_csv(augmented_data, "~/TargetedMSM/pnas_data.csv")

     
