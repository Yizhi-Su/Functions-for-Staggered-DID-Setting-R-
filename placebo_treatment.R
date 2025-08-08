# placebo treatment in staggered DID setting
rm(list = ls())


library(haven)
library(dplyr)
library(fixest)
library(DIDmultiplegt)
library(dplyr)

## --------------------------TWFE version-------------------------- ##
placebo_twfe_generator <- function(dat, year_var, id_var, treat_var, outcome_var, n = 500,
                                   treat_first_years_exclude = NULL,
                                   fe_formula = NULL,     # optional: specify a custom fixed effects structure, such as "id + year"
                                   cluster_var = NULL     # optional: cluster variable
) {
  
  df <- dat
  
  # Step 1: identify never-treated individuals
  never_treated_ids <- df %>%
    group_by(.data[[id_var]]) %>%
    summarize(ever_treated = sum(.data[[treat_var]], na.rm = TRUE) > 0, .groups = "drop") %>%
    filter(!ever_treated) %>%
    pull(.data[[id_var]])
  
  # Step 2: identify the first year of treatment
  treat_first_years <- df %>%
    group_by(.data[[id_var]]) %>%
    filter(.data[[treat_var]] == 1) %>%
    summarize(first_year = min(.data[[year_var]]), .groups = "drop") %>%
    pull(first_year) %>%
    unique()
  
  # Step 3: optionally exclude actual treatment years that should not be included in the placebo simulation
  if (!is.null(treat_first_years_exclude)) {
    treat_first_years <- treat_first_years[!treat_first_years %in% treat_first_years_exclude]
  }
  
  # Step 4: placebo TWFE
  single_placebo_run <- function() {
    df_placebo <- df
    df_placebo$placebo_treat <- 0
    
    for (y in treat_first_years) {
      # the number of individuals treated in a specific year
      n_treated <- df %>%
        group_by(.data[[id_var]]) %>%
        filter(.data[[treat_var]] == 1) %>%
        summarize(first_year = min(.data[[year_var]]), .groups = "drop") %>%
        filter(first_year == y) %>%
        nrow()
      
      # get placbo treated sample from never-treated ones
      sampled_ids <- sample(never_treated_ids, n_treated)
      
      # assign placebo treatment
      df_placebo <- df_placebo %>%
        mutate(placebo_treat = ifelse(.data[[id_var]] %in% sampled_ids & .data[[year_var]] >= y, 1, placebo_treat))
    }
    
    # Step 5: use TWFE
    # The default fixed effects structure is: id + year.
    if (is.null(fe_formula)) {
      fe_formula <- paste0(id_var, " + ", year_var)
    }
    
    formula_str <- paste0(outcome_var, " ~ placebo_treat | ", fe_formula)
    fe_model <- feols(as.formula(formula_str), data = df_placebo,
                      cluster = cluster_var)
    
    # get the estimated coefficient
    return(coef(fe_model)["placebo_treat"])
  }
  
  # Step 6: repeat for n times
  placebo_estimates <- replicate(n, single_placebo_run())
  
  return(placebo_estimates)
}





## --------------------------DID-multiplegt version-------------------------- ##
placebo_did_generator <- function(dat, year_var, id_var, treat_var, outcome_var, n = 500,
                                  treat_first_years_exclude = NULL,   # optional: exclude actual treatment years from placebo simulation
                                  effect_length = 5,                  # optional: number of dynamic effect periods
                                  ...) {                              # other optional parameters passed to did_multiplegt_dyn
  
  df <- dat
  
  # Step 1: identify never-treated individuals
  never_treated_ids <- df %>%
    group_by(.data[[id_var]]) %>%
    summarize(ever_treated = sum(.data[[treat_var]], na.rm = TRUE) > 0, .groups = "drop") %>%
    filter(!ever_treated) %>%
    pull(.data[[id_var]])
  
  # Step 2: identify the first year of treatment for treated units
  treat_first_years <- df %>%
    group_by(.data[[id_var]]) %>%
    filter(.data[[treat_var]] == 1) %>%
    summarize(first_year = min(.data[[year_var]]), .groups = "drop") %>%
    pull(first_year) %>%
    unique()
  
  # Step 3: optionally exclude actual treatment years that should not be included in the placebo simulation
  if (!is.null(treat_first_years_exclude)) {
    treat_first_years <- treat_first_years[!treat_first_years %in% treat_first_years_exclude]
  }
  
  # Step 4: placebo DID using did_multiplegt_dyn
  single_placebo_run <- function() {
    df_placebo <- df
    df_placebo$placebo_treat <- 0
    
    for (y in treat_first_years) {
      # the number of individuals treated in a specific year
      n_treated <- df %>%
        group_by(.data[[id_var]]) %>%
        summarize(first_year = min(.data[[year_var]][.data[[treat_var]] == 1], na.rm = TRUE),
                  .groups = "drop") %>%
        filter(first_year == y) %>%
        nrow()
      
      # get placebo treated sample from never-treated ones
      sampled_ids <- sample(never_treated_ids, n_treated)
      
      # assign placebo treatment from year y onward
      df_placebo <- df_placebo %>%
        mutate(placebo_treat = ifelse(.data[[id_var]] %in% sampled_ids & .data[[year_var]] >= y, 1, placebo_treat))
    }
    
    # Step 5: estimate placebo effect using did_multiplegt_dyn
    model <- did_multiplegt_dyn(
      df_placebo,
      outcome_var,
      id_var,
      year_var,
      placebo_treat,
      effects = effect_length,
      placebo = 1,
      only_never_switchers = TRUE,
      graph_off = TRUE
    )
    
    # get the estimated placebo ATE
    return(model$results$ATE[1])
  }
  
  # Step 6: repeat for n times
  placebo_estimates <- replicate(n, single_placebo_run())
  
  return(placebo_estimates)
}