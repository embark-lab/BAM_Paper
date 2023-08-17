library(dplyr)
library(haven)
load('BAM_redcap_long_merged.2023-08-11.RData')
BAM_redcap_long_merged$group <- ifelse(BAM_redcap_long_merged$group == "FP", "BAM", BAM_redcap_long_merged$group)

# Create data with pre-post completers 
# Find IDs of individuals who have at least baseline and post timepoints
post_completers <- unique(with(BAM_redcap_long_merged, id[ave(timepoint %in% c("baseline", "post"), id, FUN = function(x) sum(x) >= 2)]))

# Subset the dataset using these IDs
post_subset <- BAM_redcap_long_merged[BAM_redcap_long_merged$id %in% post_completers, ]

fu_completers <- unique(with(BAM_redcap_long_merged, id[ave(timepoint %in% c("baseline", "8wk"), id, FUN = function(x) sum(x) >= 2)]))

fu_subset <- BAM_redcap_long_merged[BAM_redcap_long_merged$id %in% fu_completers, ]

target_vars <- c('gffs_sum_25', 'umb_total_sum_25', 'epsi_neg_obese_25', 'sataq_average_25', 'ibssr_mean_25')



# Function to calculate Cohen's d effect size
calculate_cohens_d <- function(mean1, mean2, sd_pooled) {
  mean_diff <- mean1 - mean2
  return(mean_diff / sd_pooled)
}


# Function to calculate pooled standard deviation
pooled_sd <- function(n1, sd1, n2, sd2) {
  return(sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))
}

effect_sizes_list <- list()

for (i in target_vars) {
  # Filter the data for the current target variable and only completers
  filtered_data <- BAM_redcap_long_merged %>%
    filter(id %in% post_completers) %>%
    select(group, timepoint, !!sym(i))
  

  
  # Calculate the means and standard deviations for each timepoint within each group
  means_sd <- filtered_data %>%
    group_by(group, timepoint) %>%
    summarise(
      Mean = mean(!!sym(i), na.rm = TRUE),
      SD = sd(!!sym(i), na.rm = TRUE)
    )
  
  # Extract the mean and standard deviation for each timepoint
  baseline_mean <- means_sd %>%
    filter(timepoint == "baseline") %>%
    pull(Mean)
  post_mean <- means_sd %>%
    filter(timepoint == "post") %>%
    pull(Mean)
  week8_mean <- means_sd %>%
    filter(timepoint == "8wk") %>%
    pull(Mean)
  
  baseline_sd <- means_sd %>%
    filter(timepoint == "baseline") %>%
    pull(SD)
  post_sd <- means_sd %>%
    filter(timepoint == "post") %>%
    pull(SD)
  week8_sd <- means_sd %>%
    filter(timepoint == "8wk") %>%
    pull(SD)

  # Calculate Cohen's d effect size for baseline-post and baseline-8wk comparisons within each group
  cohen_d_baseline_post <- calculate_cohens_d(post_mean, baseline_mean, pooled_sd(nrow(filtered_data %>% filter(timepoint == "post")), post_sd, nrow(filtered_data %>% filter(timepoint == "baseline")), baseline_sd))
  
  
  cohen_d_baseline_8_week <- calculate_cohens_d(week8_mean, baseline_mean, pooled_sd(nrow(filtered_data %>% filter(timepoint == "8wk")), week8_sd, nrow(filtered_data %>% filter(timepoint == "baseline")), baseline_sd))
  
  
  # Create a data frame to store the effect sizes for the current target variable
  result <- data.frame(
    group = c("BAM","BP"),
    Variable = i,
    Cohens_d_baseline_post = cohen_d_baseline_post,
    Cohens_d_baseline_8_week = cohen_d_baseline_8_week
  )
  
  effect_sizes_list[[i]] <- result
}

# Combine all effect sizes for different target variables into a single data frame
effect_sizes_df_postsub <- do.call(rbind, effect_sizes_list)
effect_sizes_df_postsub <- effect_sizes_df_postsub |> 
  select(-Cohens_d_baseline_8_week)

n <- unique(post_subset$id)


effect_sizes_list_8wk <- list()

for (i in target_vars) {
  # Filter the data for the current target variable and only completers
  filtered_data <- BAM_redcap_long_merged %>%
    filter(id %in% fu_completers) %>%
    select(group, timepoint, !!sym(i))
  
  
  
  # Calculate the means and standard deviations for each timepoint within each group
  means_sd <- filtered_data %>%
    group_by(group, timepoint) %>%
    summarise(
      Mean = mean(!!sym(i), na.rm = TRUE),
      SD = sd(!!sym(i), na.rm = TRUE)
    )
  
  # Extract the mean and standard deviation for each timepoint
  baseline_mean <- means_sd %>%
    filter(timepoint == "baseline") %>%
    pull(Mean)
  post_mean <- means_sd %>%
    filter(timepoint == "post") %>%
    pull(Mean)
  week8_mean <- means_sd %>%
    filter(timepoint == "8wk") %>%
    pull(Mean)
  
  baseline_sd <- means_sd %>%
    filter(timepoint == "baseline") %>%
    pull(SD)
  post_sd <- means_sd %>%
    filter(timepoint == "post") %>%
    pull(SD)
  week8_sd <- means_sd %>%
    filter(timepoint == "8wk") %>%
    pull(SD)
  
  # Calculate Cohen's d effect size for baseline-post and baseline-8wk comparisons within each group
  cohen_d_baseline_post <- calculate_cohens_d(post_mean, baseline_mean, pooled_sd(nrow(filtered_data %>% filter(timepoint == "post")), post_sd, nrow(filtered_data %>% filter(timepoint == "baseline")), baseline_sd))
  
  
  cohen_d_baseline_8_week <- calculate_cohens_d(week8_mean, baseline_mean, pooled_sd(nrow(filtered_data %>% filter(timepoint == "8wk")), week8_sd, nrow(filtered_data %>% filter(timepoint == "baseline")), baseline_sd))
  
  
  # Create a data frame to store the effect sizes for the current target variable
  result <- data.frame(
    group = c("BAM","BP"),
    Variable = i,
    Cohens_d_baseline_post = cohen_d_baseline_post,
    Cohens_d_baseline_8_week = cohen_d_baseline_8_week
  )
  
  effect_sizes_list_8wk[[i]] <- result
}


# Combine all effect sizes for different target variables into a single data frame
effect_sizes_df_fusub <- do.call(rbind, effect_sizes_list)
effect_sizes_df_fusub <- effect_sizes_df_fusub |> 
  select(-Cohens_d_baseline_post)

subset_effects <- full_join(effect_sizes_df_postsub, effect_sizes_df_fusub)


n_fu <- unique(fu_subset$id)

