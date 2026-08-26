source("config/paths.R")


library(tidyverse)

resdir <- bmrc_path("results/ldsc/")
files_2 <- list.files(resdir, pattern = "\\.log$", full.names = TRUE)

results_trait <- map_dfr(files_2, extract_ldsc) %>%
  mutate(fdr = p.adjust(p),
         sig2 = ifelse(fdr < 0.05, "•", ifelse(p < 0.05, "·", ""))) 

x_levels <- tribble(~trait,    ~name, ~type,
                    "Circadian Acceleration (CA)","res_GWAS",                "Blood",
                    #"Chronotype (binary)", "Jones_2018_morning_person", "Self-report",
                    #"Chronotype (binary) 2", "chrono_b", "Self-report",
                    "Chronotype (quantitative)", "chrono_v2_GWAS",      "Self-report",
                    "Relative amplitude", "RAgroup_lin", "Accelerometer",
                    "L5 time","ACC_L5_TIME", "Accelerometer",
                    "M10 time","ACC_M10_TIME", "Accelerometer",
                    "Sleep midpoint", "ACC_SLEEP_MIDP", "Accelerometer",
                    "Diurnal inactivity","ACC_DIURNAL_INACT", "Accelerometer",
                    "Sleep duration","ACC_SLEEP_DUR", "Accelerometer",
                    "Sleep duration sd","ACC_SLEEP_DUR_SD", "Accelerometer",
                    "Sleep efficiency","ACC_SLEEP_EFF", "Accelerometer",
                    "Num sleep episodes","ACC_N_SLEEP_EPISODES", "Accelerometer"
)

trait_map <- tribble(
  ~trait,          ~label,                          ~file, ~class,
  "AD",            "Alzheimer’s disease",            "AD_sumstats_Jansenetal_2019sept_tidy", "Psych",
  "T2D",           "Type 2 diabetes",                "EUR_Metal_LDSC-CORR_Neff", "Cardiometabolic",
  "RA_autoimmune", "Rheumatoid arthritis",           "GCST90132223_buildGRCh37", "Immune",
  "DBP",           "DBP",       "BP-ICE_EUR_DBP_transformed_15-04-2020", "Cardiometabolic",
  "HTN",           "Hypertension",                   "BP-ICE_EUR_HTN_15-04-2020", "Cardiometabolic",
#  "PP",            "Pulse pressure",                 "BP-ICE_EUR_PP_transformed_15-04-2020", "Cardiometabolic",
#  "SBP",           "Systolic blood pressure",        "BP-ICE_EUR_SBP_transformed_15-04-2020", "Cardiometabolic",

  "CAD",           "Coronary artery disease",        "GCST90132314_buildGRCh37", "Cardiometabolic",
  "mtDNA_CN",      "mtDNA copy number",              "mtDNA_CN", "Other",
  "BMI",           "Body mass index",                "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED", "Cardiometabolic",
  "FG",            "Fasting glucose",                "MAGIC_FastingGlucose", "Cardiometabolic",
  "HipAdjBMI",     "BMI-adjusted hip circumference", "GCST012229", "Cardiometabolic",
  "pgc-mdd2024_Clin_eur",      "MDD Clin.",    "pgc-mdd2024_Clin_eur",                                     "Psych",
  "pgc-mdd2024_EHR_eur",       "MDD EHR",     "pgc-mdd2024_EHR_eur",                                      "Psych",
  "pgc-mdd2024_Quest_eur",     "MDD Quest",   "pgc-mdd2024_Quest_eur",                                    "Psych",
  "PGC3_SCZ_wave3",            "Schizophrenia",                      "PGC3_SCZ_wave3",                                           "Psych",
  "bip2024",                  "Bipolar disorder",              "bip2024",                                                  "Psych",
  "Insomnia", "Insomnia", "Jansen_2018_insomnia_broad", "Sleep"
)

out<- results_trait %>%
  mutate(p1 = factor(p1_name, levels = x_levels$name, labels = x_levels$trait),
         p2 = factor(p2_name, levels = trait_map$file, labels = trait_map$label)) %>%
  filter(!is.na(p1)) %>% filter(!is.na(p2)) %>%
  filter(p1 == "Circadian Acceleration (CA)") %>%
  select(p1, p2, rg, se, z, p, fdr)

pcomp <- results_trait %>%
  #filter(p < 0.05) %>%
  mutate(p1 = factor(p1_name, levels = x_levels$name, labels = x_levels$trait),
         p2 = factor(p2_name, levels = trait_map$file, labels = trait_map$label)) %>%
  left_join(x_levels, by = c("p1_name" = "name")) %>%
  filter(!is.na(p1)) %>%
  filter(!is.na(p2)) %>%
  left_join(trait_map, by = c("p2" = "label")) %>%
  mutate(class = factor(class, levels = c( "Sleep","Psych","Cardiometabolic", "Immune", "Other")),
         type = factor(type, levels = c("Blood", "Self-report", "Accelerometer"))) %>%
  ggplot(aes(x = p1, y = p2, fill = rg)) +
  geom_tile(aes(fill = rg), color = "black", linewidth = 0.3) +
  geom_text(aes(label = sig2), size = 4) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-1, 1),
    name = "rg"
  ) +
  scale_x_discrete(position = "top") +
  #scale_y_discrete(expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  facet_grid(rows = vars(class),cols = vars(type), scales = "free", space = "free", switch = "x") +
  #facet_grid(cols = vars(type), scales = "free", space = "free",switch = "x") +
  theme_classic(base_size = 10) +
  theme(#panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.placement  = "inside",
        legend.position = "none",
        #panel.spacing.y = unit(0.2, "lines"),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        #strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1),
        plot.margin            = margin(1, r = 20, 1, 1))
pcomp
ggsave(paste0("OUTPUTS/rg_traits.png"), pcomp, width = 6.5, height = 8.2)
