
split_sex_age_r2 <- function(data,
                             main1  = "sex",
                             main2  = "age_recruitment",
                             r2_col = "pr2",          # change to "t_r2" if that’s your column name
                             drop_interaction = TRUE) # keep it? FALSE ⇒ keep with pr2 = 0
{
  # tidy-eval handle for the R² column
  r2 <- sym(r2_col)

  data %>%
    group_by(phen, title) %>%
    mutate(
      is_inter = term == "sex:age_recruitment",

      # one interaction R² per group (0 if none)
      inter_r2 = max(if_else(is_inter, !!r2, 0), na.rm = TRUE),

      # allocate half-half, zero the interaction itself
      !!r2_col := case_when(
        term == str_to_lower(main1) ~ coalesce(!!r2, 0) + inter_r2 / 2,
        term == str_to_lower(main2) ~ coalesce(!!r2, 0) + inter_r2 / 2,
        is_inter                          ~ 0,
        TRUE                              ~ !!r2
      )
    ) %>%
    { if (drop_interaction) filter(., !is_inter) else . } %>%      # drop or keep
    select(-inter_r2, -is_inter) %>%
    ungroup()
}
