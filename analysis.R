library(tidyverse)
library(dragon)
library(cowplot)

# Settings and functions -----------------------------------
theme_set(theme_light() + theme(legend.position = "none"))
pauling_threshold <- 0.11 # for wMEEcv low vs high

add_time_groups <- function(df, column)
{
  #The time intervals should be changed to:
  #  Group 1: 4.34 < time <= 1.8
  #  Group 2: 1.8 < time <= 0.6
  #  Group 3: < 0.6
  df %>%
    # case_when() works
    dplyr::mutate(time_group = case_when(
      {{column}} >= 1.8 & {{column}} < 4.34 ~ "Group 1",
      {{column}} < 1.8 & {{column}} >= 0.6 ~ "Group 2",
      {{column}} < 0.6 ~ "Group 3")) %>%
    # remove >=4.34 (the NAs)
    drop_na(time_group) 
}

anova_and_plot <- function(df, title)
{
  aov(w_cov_pauling ~ time_group, data = df) %>% 
    TukeyHSD() %>%
    broom::tidy() %>%
    ggplot(aes(y = contrast,
               x = estimate, 
               xmin = conf.low,
               xmax = conf.high)) +
    geom_linerange(color = "red") + 
    geom_point() +
    geom_vline(xintercept=0) + 
    xlab("Mean difference in wMEEcv (with 95% CI)") +
    ggtitle(title)
}



wmeecv_boxplot <- function(df, title)
{
  ggplot(df, aes(x = time_group, 
                 y = w_cov_pauling, 
                 fill = time_group)) + 
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_brewer(palette = "Pastel1") +
    ylab("wMEEcv") + 
    xlab("Time group") +
    ggtitle(title) 
}


high_to_low_barplots <- function(df, title)
{
  df %>%
    dplyr::count(time_group, wmeecv_discrete) %>%
    tidyr::pivot_wider(names_from = wmeecv_discrete, values_from = n) %>%
    dplyr::mutate(high_to_low_ratio = high/low,
                  total = high + low) %>%
    ggplot(aes(x = time_group, fill = time_group, y = high_to_low_ratio)) +
    geom_col(color = "black") + 
    scale_fill_brewer(palette = "Pastel1") +
    geom_text(aes(y = high_to_low_ratio + 1, label = glue::glue("Total minerals: {total}"))) +
    labs(x = "Time (Gy)", y = "Ratio of number of high/low wMEEcv minerals formed") + 
    ggtitle(title)
}



# Prepare data ---------------------------------


# Build the full mineral network
dragon:::element_info %>%
  dplyr::pull(element) -> all_elements
full_network <- dragon:::initialize_network(elements_of_interest = all_elements)

# Grab minerals for analysis and add time groups - MAX ONLY
full_network$nodes %>%
  dplyr::filter(group == "mineral") %>%
  dplyr::select(name = id, w_cov_pauling, max_age) %>%
  add_time_groups(max_age) %>%
  dplyr::mutate(wmeecv_discrete = ifelse(w_cov_pauling >= pauling_threshold, "high", "low")) %>%
  tidyr::drop_na() -> max_minerals_time

# Now, ALL minerals regardless of max age
full_network$locality_info %>%
  dplyr::select(name = mineral_name, max_age_locality) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(dummycount = 1:dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    max_minerals_time %>% select(name,w_cov_pauling,wmeecv_discrete)
  ) %>%
  dplyr::distinct() %>%
  add_time_groups(max_age_locality) -> all_minerals_time


# Boxplot of wMEEcv values for each time group ---------------------------------
max_boxes <- wmeecv_boxplot(max_minerals_time, "Maximum age per mineral only")
all_boxes <- wmeecv_boxplot(all_minerals_time, "All minerals formed")
box_grid <- plot_grid(max_boxes, all_boxes, nrow=1, labels="auto")
save_plot("meecv_boxplots.png", box_grid, base_width = 8, base_height = 3)

# ANOVA coefficients among time groups ------------------------------------------
max_aov_plot <- anova_and_plot(max_minerals_time, "Maximum age per mineral only")
all_aov_plot <- anova_and_plot(all_minerals_time, "All minerals formed")
aov_grid <- plot_grid(max_aov_plot, all_aov_plot, nrow=1, labels="auto", scale=0.95)
save_plot("anova_results.png", aov_grid, base_width = 8, base_height = 3)


# Ratio of high/low minerals formed among time groups ------------------------------------------
max_ratio_barplot <- high_to_low_barplots(max_minerals_time, "Maximum age per mineral only")
all_ratio_barplot <- high_to_low_barplots(all_minerals_time, "All minerals formed")
barplot_grid <- plot_grid(max_ratio_barplot, all_ratio_barplot, nrow=1, labels="auto", scale=0.95)
save_plot("barplot_ratios.png", barplot_grid, base_width = 12, base_height = 5)

# Chi-squared test ----------------------------------------------------
max_table <- table(max_minerals_time$time_group, max_minerals_time$wmeecv_discrete)
chisq.test(max_table) # p-value = 7.422e-12

all_table <- table(all_minerals_time$time_group, all_minerals_time$wmeecv_discrete)
chisq.test(all_table) #  p-value < 2.2e-16
