library(tidyverse)
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

dat <- read_dta("example/data/MFPS_PNAS_MAIN_ITT_DATA.dta") %>%
  mutate_at(vars(treatment, curr_use_inj_adj_bl, area_bl, age_group_bl, ever_use_bl, edu_primary_bl, work_bl, religion_r_bl, ethnicity_r_bl), as.factor) %>%
  filter(year == 2018) %>%
  select(curr_use_inj_adj, treatment, !!covariates) %>%
  drop_na()

cairo_pdf("~plots/pnas_effect_modifier.pdf", width = 9, height = 4)
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
lm(curr_use_inj_adj ~ 1 + treatment, data = dat)

# Adjusted
f <- as.formula(glue::glue("curr_use_inj_adj ~ 1 + treatment + {str_c(covariates, collapse = ' + ')}"))
lm(f, data = dat)

write_csv(dat, "~example/data/pnas_data.csv")


