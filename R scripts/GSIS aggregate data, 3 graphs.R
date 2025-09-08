# ---- Load Required Libraries ----
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
# ---- Optional: Install writexl if not already installed ----
# install.packages("writexl")
library(writexl)


# ---- Load the Excel File ----
df <- read_excel("C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//NormalizedStims.xlsx")

# ---- Filter Out Omitted Samples ----
df <- df %>%
  filter(tolower(trimws(Omit)) != "y" | is.na(Omit))

# ---- Define Dose Order and Custom Alpha Values ----
dose_order <- c("10uM", "1uM", "100nM")
alpha_map <- c("10uM" = 1,
               "1uM" = 1,
               "100nM" = 1,
               "NA" = 1.0)

# ---- Add Grouping and Alpha Columns ----
df <- df %>%
  mutate(Dose = factor(Dose, levels = dose_order),
         Group = paste0(Sample, "_", Dose),
         AlphaVal = alpha_map[as.character(Dose)])

# ---- Order Group by Sample and Dose ----
group_order_df <- df %>%
  distinct(Sample, Dose, Group) %>%
  arrange(
    case_when(
      grepl("DMSO", Sample) ~ 0,
      TRUE ~ 1
    ),
    Sample,
    Dose
  )

df$Group <- factor(df$Group, levels = group_order_df$Group)

# ---- Separate DMSO background ----
dmso_vals <- df %>%
  filter(Sample == "DMSO") %>%
  pull(Log2)

# ---- One-sided empirical p-value: test for INCREASE vs DMSO ----
# H0: mean(treatment) <= mean(DMSO); HA: mean(treatment) > mean(DMSO)
# We estimate the null by resampling means from dmso_vals with the same n as the treatment group.
emp_pval_increase <- function(x, ref, B = 100000) {
  x <- x[is.finite(x)]
  ref <- ref[is.finite(ref)]
  if (length(x) < 1 || length(ref) < 2) return(NA_real_)
  obs <- mean(x)
  # sample with replacement if DMSO pool is smaller than group size
  replace_flag <- length(ref) < length(x)
  sim_means <- replicate(B, mean(sample(ref, length(x), replace = replace_flag)))
  # conservative one-sided p-value (include ties)
  (sum(sim_means >= obs) + 1) / (B + 1)
}

# ---- Calculate empirical p-values per group (keep Dose) ----
treatment_stats <- df %>%
  filter(Sample != "DMSO") %>%
  group_by(Group) %>%
  summarise(
    Sample   = first(Sample),
    Dose     = first(Dose),             # <— keep Dose for easy filtering later
    n        = n(),
    mean_log2= mean(Log2, na.rm = TRUE),
    pval     = emp_pval(Log2, dmso_vals),
    .groups  = "drop"
  ) %>%
  mutate(
    pval_adj = p.adjust(pval, method = "BH"),
    sig_label = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

# ---- Save p-value summary to Excel ----
write_xlsx(treatment_stats, path = "C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//empirical_p_values.xlsx")

make_plot_for_dose <- function(dose_label) {
  df_sub <- df %>%
    filter(grepl("DMSO|Untreated", Sample) | Dose == dose_label) %>%
    droplevels()
  
  y_pos_sub <- df_sub %>%
    group_by(Group) %>%
    summarise(y = max(Log2, na.rm = TRUE) + 0.1, .groups = "drop")
  
  label_df_sub <- treatment_stats %>%
    filter(Dose == dose_label) %>%
    left_join(y_pos_sub, by = "Group")
  
  df_sub$Group <- droplevels(df_sub$Group)
  
  # Create a mapping from Group -> Sample for axis labels
  x_labels <- df_sub %>%
    distinct(Group, Sample) %>%
    arrange(Group)
  
  p <- ggplot(df_sub, aes(x = Group, y = Log2,
                          color = Sample,
                          shape = `Cell Source`,
                          alpha = AlphaVal)) +
    geom_jitter(width = 0.3, size = 2, show.legend = FALSE) +   # hide legend here
    
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 color = "black", width = 0.5, inherit.aes = FALSE,
                 aes(x = Group, y = Log2), show.legend = FALSE) +
    stat_summary(fun = mean, geom = "crossbar",
                 fatten = 0, width = 0.8, color = "black", inherit.aes = FALSE,
                 aes(x = Group, y = Log2), show.legend = FALSE) +
    
    geom_text(data = label_df_sub, aes(x = Group, y = y, label = sig_label),
              inherit.aes = FALSE, size = 5, vjust = 0) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    
    theme_minimal() +
    labs(title = paste0("Log2 Index by Treatment — Controls + ", dose_label),
         y = "Log2 Normalized STIM Index",
         x = "Treatment",
         shape = "Cell Source") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    
    scale_x_discrete(labels = setNames(x_labels$Sample, x_labels$Group)) +
    
    # Hide legends for color + alpha, keep only Cell Source (shape)
    guides(color = "none", alpha = "none")
  
  return(p)
}


# ---- Build and save the three separate plots ----
doses_to_plot <- c("100nM", "1uM", "10uM")
plots <- lapply(doses_to_plot, make_plot_for_dose)

# Print them to the active device (RStudio will show all three)
plots[[1]]
plots[[2]]
plots[[3]]

