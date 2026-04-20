rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(readr)
library(readxl)
library(ggsignif)

# load aging data
basedir <- "C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS"

# Load all datasets and store in lists
maple_datasets <- list(
  cov1 = read.csv(file.path(basedir, "maple_output_for_normalized_cov1.csv")),
  cov5 = read.csv(file.path(basedir, "maple_output_for_normalized_cov5.csv")),
  cov10 = read.csv(file.path(basedir, "maple_output_for_normalized_cov10.csv"))
)

ic_datasets <- list(
  cov1 = read.csv(file.path(basedir, "predictions-2026-04-10_cov1.csv")),
  cov5 = read.csv(file.path(basedir, "predictions-2026-04-10_cov5.csv")),
  cov10 = read.csv(file.path(basedir, "predictions-2026-04-10_cov10.csv"))
)


# Make Group column from Sample or sample_id column
extract_group_maple <- function(data) {
  data %>% 
    dplyr::mutate(Group = sample_id) %>%
    dplyr::mutate(Group = gsub("-.*", "", Group)) %>%
    dplyr::mutate(Group = factor(Group, levels = c("NIR", "IR2Gy24h", "IR2Gy6d", "IR10Gy24h", "IR10Gy6d")))
}
extract_group_ic <- function(data) {
  data %>% 
    dplyr::mutate(Group = Sample) %>%
    dplyr::mutate(Group = gsub("\\..*", "", Group)) %>%
    dplyr::mutate(Group = factor(Group, levels = c("NIR", "IR2Gy24h", "IR2Gy6d", "IR10Gy24h", "IR10Gy6d")))
}

# Extract Dose and Time variables from Group factor for analysis
# Group levels: NIR, IR2Gy24h, IR2Gy6d, IR10Gy24h, IR10Gy6d
extract_dose_time <- function(data) {
  data %>%
    dplyr::mutate(
      Dose = dplyr::case_when(
        Group == "NIR" ~ 0,
        grepl("2Gy", Group) ~ 2,
        grepl("10Gy", Group) ~ 10,
        TRUE ~ NA_real_
      ),
      Time = dplyr::case_when(
        Group == "NIR" ~ 0,
        grepl("24h", Group) ~ 1,
        grepl("6d", Group) ~ 6,
        TRUE ~ NA_real_
      )
    )
}

# Apply group extraction and dose/time variables to all datasets
maple_datasets <- lapply(maple_datasets, function(data) {
  extract_group_maple(data) %>% extract_dose_time()
})

ic_datasets <- lapply(ic_datasets, function(data) {
  extract_group_ic(data) %>% extract_dose_time()
})

# Redirect all output to text file
output_file <- file.path(basedir, "statistical_test_results_aging.txt")
sink(output_file, append = FALSE, split = TRUE)
cat("Statistical Analysis Results\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Linear models for MAPLE data: test effects of Time, Dose, and Time*Dose on Pred_Age
for (cov_name in names(maple_datasets)) {
  cat("\n=== Linear Model for maple_", cov_name, ": Pred_Age ~ Dose * Time ===\n", sep = "")
  lm_model <- lm(Pred_Age ~ Dose * Time, data = maple_datasets[[cov_name]])
  print(summary(lm_model))
  print(anova(lm_model))
}

# Linear models for IC data: test effects of Time, Dose, and Time*Dose on DNAmIC
for (cov_name in names(ic_datasets)) {
  cat("\n=== Linear Model for ic_", cov_name, ": DNAmIC ~ Dose * Time ===\n", sep = "")
  lm_model <- lm(DNAmIC ~ Dose * Time, data = ic_datasets[[cov_name]])
  print(summary(lm_model))
  print(anova(lm_model))
}

# ANOVA testing Group effects with post hoc tests using NIR as reference
library(multcomp)

# ANOVA for MAPLE data
for (cov_name in names(maple_datasets)) {
  cat("\n=== ANOVA for maple_", cov_name, ": Pred_Age ~ Group ===\n", sep = "")
  aov_model <- aov(Pred_Age ~ Group, data = maple_datasets[[cov_name]])
  print(summary(aov_model))
  
  cat("\n--- Post hoc test (Dunnett): All groups vs NIR ---\n")
  posthoc_model <- glht(aov_model, linfct = mcp(Group = "Dunnett"))
  print(summary(posthoc_model))
}

# ANOVA for IC data
for (cov_name in names(ic_datasets)) {
  cat("\n=== ANOVA for ic_", cov_name, ": DNAmIC ~ Group ===\n", sep = "")
  aov_model <- aov(DNAmIC ~ Group, data = ic_datasets[[cov_name]])
  print(summary(aov_model))
  
  cat("\n--- Post hoc test (Dunnett): All groups vs NIR ---\n")
  posthoc_model <- glht(aov_model, linfct = mcp(Group = "Dunnett"))
  print(summary(posthoc_model))
}

# Pairwise comparison: IR10Gy6d vs NIR only
cat("\n\n========================================\n")
cat("PAIRWISE COMPARISON: IR10Gy6d vs NIR\n")
cat("========================================\n")

# MAPLE datasets
for (cov_name in names(maple_datasets)) {
  subset_data <- maple_datasets[[cov_name]] %>% dplyr::filter(Group %in% c("NIR", "IR10Gy6d"))
  
  cat("\n=== maple_", cov_name, ": Pred_Age (IR10Gy6d vs NIR) ===\n", sep = "")
  cat("\n--- Two-sample t-test (equal variances) ---\n")
  print(t.test(Pred_Age ~ Group, data = subset_data, var.equal = TRUE))
  
  cat("\n--- Welch's t-test (unequal variances) ---\n")
  print(t.test(Pred_Age ~ Group, data = subset_data, var.equal = FALSE))
  
  cat("\n--- Mann-Whitney U test (non-parametric) ---\n")
  print(wilcox.test(Pred_Age ~ Group, data = subset_data))
}

# IC datasets
for (cov_name in names(ic_datasets)) {
  subset_data <- ic_datasets[[cov_name]] %>% dplyr::filter(Group %in% c("NIR", "IR10Gy6d"))
  
  cat("\n=== ic_", cov_name, ": DNAmIC (IR10Gy6d vs NIR) ===\n", sep = "")
  cat("\n--- Two-sample t-test (equal variances) ---\n")
  print(t.test(DNAmIC ~ Group, data = subset_data, var.equal = TRUE))
  
  cat("\n--- Welch's t-test (unequal variances) ---\n")
  print(t.test(DNAmIC ~ Group, data = subset_data, var.equal = FALSE))
  
  cat("\n--- Mann-Whitney U test (non-parametric) ---\n")
  print(wilcox.test(DNAmIC ~ Group, data = subset_data))
}

# Close output file for statistical tests
sink()
cat("\nResults written to:", output_file, "\n")

# Calculate summary statistics (mean and SE) for plotting
calc_summary <- function(data, response_var) {
  data %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      mean = mean(!!rlang::sym(response_var), na.rm = TRUE),
      se = sd(!!rlang::sym(response_var), na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
}

# Create plotting function
create_plot <- function(data, response_var, dataset_name, pval) {
  summary_data <- calc_summary(data, response_var)
  
  y_min <- min(summary_data$mean - summary_data$se) * 0.95
  y_max <- max(summary_data$mean + summary_data$se) * 1.05
  
  y_label <- paste(response_var, "(mean ± SE)")
  
  ggplot(summary_data, aes(x = Group, y = mean, fill = Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 0.8) +
    geom_signif(comparisons = list(c("NIR", "IR10Gy6d")),
                annotations = sprintf("p = %.3g", pval),
                y_position = y_max * 0.98,
                tip_length = 0.02,
                textsize = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(title = dataset_name,
         x = "Group",
         y = y_label) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

# Generate plots for all datasets
maple_plots <- list()
ic_plots <- list()

for (cov_name in names(maple_datasets)) {
  subset_data <- maple_datasets[[cov_name]] %>% dplyr::filter(Group %in% c("NIR", "IR10Gy6d"))
  pval <- t.test(Pred_Age ~ Group, data = subset_data, var.equal = TRUE)$p.value
  
  plot_name <- paste0("MAPLE ", toupper(cov_name))
  maple_plots[[cov_name]] <- create_plot(maple_datasets[[cov_name]], "Pred_Age", plot_name, pval)
}

for (cov_name in names(ic_datasets)) {
  subset_data <- ic_datasets[[cov_name]] %>% dplyr::filter(Group %in% c("NIR", "IR10Gy6d"))
  pval <- t.test(DNAmIC ~ Group, data = subset_data, var.equal = TRUE)$p.value
  
  plot_name <- paste0("IC ", toupper(cov_name))
  ic_plots[[cov_name]] <- create_plot(ic_datasets[[cov_name]], "DNAmIC", plot_name, pval)
}

# Combine plots in 2x3 grid: MAPLE row (top), IC row (bottom)
combined_plot <- plot_grid(
  maple_plots$cov1, maple_plots$cov5, maple_plots$cov10,
  ic_plots$cov1, ic_plots$cov5, ic_plots$cov10,
  ncol = 3, nrow = 2,
  labels = c("A", "B", "C", "D", "E", "F")
)

# Save plots
plot_file <- file.path(basedir, "aging_group_means_plots.pdf")
ggsave(plot_file, combined_plot, width = 15, height = 10)
cat("\nPlots saved to:", plot_file, "\n")



