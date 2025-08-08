match_single_variable <- function(dat, year_var, id_var, treat_var, outcome_var,
                                  match_years = 5, eval_years = 2,
                                  tolerance = 0.2) {
  df <- dat
  
  # Step 1: Identify treatment years
  treat_years <- df %>%
    filter(.data[[treat_var]] > 0) %>%
    pull(.data[[year_var]]) %>%
    unique() %>%
    sort()
  
  all_matched <- list()
  
  for (T in treat_years) {
    # Step 2: Define matching and evaluation periods
    match_period <- (T - match_years):(T - eval_years - 1)
    eval_period <- (T - eval_years):(T - 1)
    
    # Step 3: Identify treated and never-treated units
    treated_ids <- df %>%
      filter(.data[[treat_var]] == T) %>%
      pull(.data[[id_var]]) %>%
      unique()
    
    control_ids <- df %>%
      filter(.data[[treat_var]] == 0) %>%
      pull(.data[[id_var]]) %>%
      unique()
    
    # Step 4: Build matching matrices for outcome variable
    treated_matrix <- df %>%
      filter(.data[[id_var]] %in% treated_ids, .data[[year_var]] %in% match_period) %>%
      select(id = .data[[id_var]], year = .data[[year_var]], outcome = .data[[outcome_var]]) %>%
      pivot_wider(names_from = year, values_from = outcome)
    
    control_matrix <- df %>%
      filter(.data[[id_var]] %in% control_ids, .data[[year_var]] %in% match_period) %>%
      select(id = .data[[id_var]], year = .data[[year_var]], outcome = .data[[outcome_var]]) %>%
      pivot_wider(names_from = year, values_from = outcome)
    
    # Step 5: Match each treated unit to the best control
    matched_pairs <- treated_matrix %>%
      rowwise() %>%
      mutate(
        best_match = {
          this_vec <- c_across(matches("^\\d{4}"))
          distances <- control_matrix %>%
            select(-id) %>%
            as.matrix() %>%
            apply(1, function(x) sum((x - this_vec)^2, na.rm = TRUE))
          control_matrix$id[which.min(distances)]
        }
      ) %>%
      select(treated = id, control = best_match)
    
    # Step 6: Evaluate match quality in eval period
    matched_pairs <- matched_pairs %>%
      rowwise() %>%
      mutate(
        diff_check = {
          t_vals <- df %>%
            filter(.data[[id_var]] == treated, .data[[year_var]] %in% eval_period) %>%
            arrange(.data[[year_var]]) %>%
            pull(.data[[outcome_var]])
          c_vals <- df %>%
            filter(.data[[id_var]] == control, .data[[year_var]] %in% eval_period) %>%
            arrange(.data[[year_var]]) %>%
            pull(.data[[outcome_var]])
          if (length(t_vals) < eval_years || length(c_vals) < eval_years) Inf else mean(abs(t_vals - c_vals))
        }
      ) %>%
      filter(diff_check <= tolerance)
    
    # Step 7: Build matched sample
    keep_ids <- unique(c(matched_pairs$treated, matched_pairs$control))
    matched_df <- df %>%
      filter(.data[[id_var]] %in% keep_ids) %>%
      mutate(matched_year = T,
             matched_group = case_when(
               .data[[id_var]] %in% matched_pairs$treated ~ "treated",
               .data[[id_var]] %in% matched_pairs$control ~ "control",
               TRUE ~ NA_character_
             ))
    
    all_matched[[as.character(T)]] <- matched_df
  }
  
  # Step 8: Combine all years
  final_matched <- bind_rows(all_matched)
  return(final_matched)
}
