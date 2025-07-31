library(dplyr)
library(purrr)
library(survival)
install.packages("survminer")
library(survminer)

# 1. Specify the diseases you want to plot
#    Replace these with the actual column names from `dis2`
indicated_diseases <- c("130892-0.0", "130708-0.0", "131306-0.0")

# 2. Function to build the dataset and return a survminer plot
make_ci_plot <- function(disease_field) {
  df_plot <- base_cohort %>%
    left_join(dis2 %>% select(eid, !!sym(disease_field))
              %>% rename(diag_date = !!sym(disease_field)),
              by = "eid") %>%
    left_join(bio_covs %>% select(eid, abs_gap, ares_q), by = "eid") %>%
    filter(is.na(death_date) | death_date > date_bsampling,
           is.na(diag_date)  | diag_date  > date_bsampling) %>%
    mutate(
      end_date = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event    = as.integer(!is.na(diag_date) &
                              diag_date <= pmin(death_date, censor_date, na.rm = TRUE)),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25,
      ares_q   = factor(ares_q, levels = 1:5)
    ) %>%
    filter(time_yrs > 1, ares_q %in% c("1","3","5"))

  # Fit survival curves
  fit <- survfit(Surv(time_yrs, event) ~ ares_q, data = df_plot)

  # Plot cumulative incidence = 1 – survival
  ggsurvplot(
    fit,
    data        = df_plot,
    fun         = "event",
    palette     = c("steelblue", "forestgreen", "firebrick"),
    legend.title= "res_abs quintile",
    legend.labs = c("Q1","Q3","Q5"),
    xlab        = "Years since baseline",
    ylab        = "Cumulative incidence",
    title       = paste("Disease:", disease_field),
    risk.table  = TRUE,
    risk.table.height = 0.25
  )
}

# 3. Generate one plot per disease
ci_plots <- map(indicated_diseases, make_ci_plot)

# 4. Display them
# If you’re in an interactive session:
walk(ci_plots, print)

# Or, to arrange them in a grid:
library(gridExtra)
grid.arrange(grobs = map(ci_plots, `[[`, "plot"), ncol = 1)
