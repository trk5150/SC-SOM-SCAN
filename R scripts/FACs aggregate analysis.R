# ---- Load Required Libraries ----
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
# ---- Optional: Install writexl if not already installed ----
# install.packages("writexl")
library(writexl)


# ---- Load the Excel File ----
df <- read_excel("C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//FACs aggregate.xlsx")

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
      grepl("Untreated", Sample) ~ 1,
      grepl("Human Islet", Sample) ~ 2,
      TRUE ~ 3
    ),
    Sample,
    Dose
  )

df$Group <- factor(df$Group, levels = group_order_df$Group)

# ---- Separate DMSO background ----
dmso_vals <- df %>%
  filter(Sample == "DMSO") %>%
  pull(Log2)

# ---- Empirical p-value function ----
emp_pval <- function(x, ref, B = 10000) {
  obs <- mean(x)
  sim_means <- replicate(B, mean(sample(ref, length(x), replace = FALSE)))
  mean(abs(sim_means - mean(ref)) >= abs(obs - mean(ref)))
}

# ---- Calculate empirical p-values per group ----
treatment_stats <- df %>%
  filter(Sample != "DMSO") %>%
  group_by(Group) %>%
  summarise(
    Sample = first(Sample),
    n = n(),
    mean_log2 = mean(Log2, na.rm = TRUE),
    pval = emp_pval(Log2, dmso_vals)
  ) %>%
  ungroup() %>%
  mutate(
    pval_adj = p.adjust(pval, method = "BH"),
    sig_label = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) #%>%
#filter(sig_label != "")

# ---- Determine Y-Axis Position for Labels ----
y_pos <- df %>%
  group_by(Group) %>%
  summarise(y = max(Log2, na.rm = TRUE) + 0.1)

label_df <- left_join(treatment_stats, y_pos, by = "Group")

# ---- Save p-value summary to Excel ----
write_xlsx(treatment_stats, path = "C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//FACs empirical_p_values.xlsx")


# ---- Final Plot ----
ggplot(df, aes(x = Group, y = Log2,
               color = Sample,
               shape = `Cell Source`,
               alpha = AlphaVal)) +
  geom_jitter(width = 0.3, size = 2) +
  
  # Error bars (mean ± SEM)
  stat_summary(fun.data = mean_se, geom = "errorbar",
               color = "black", width = 0.8, inherit.aes = FALSE,
               aes(x = Group, y = Log2)) +
  
  # Mean point
  stat_summary(fun = mean, geom = "crossbar",
               fatten = 0, width = 1.2, color = "black", inherit.aes = FALSE,
               aes(x = Group, y = Log2)) +
  
  # Significance labels
  geom_text(data = label_df, aes(x = Group, y = y, label = sig_label),
            inherit.aes = FALSE, size = 5, vjust = 0) +
  
  # Baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  theme_minimal() +
  labs(title = "Log2 beta abundance",
       y = "Log2 Normalized beat cell abundance",
       x = "Treatment Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
