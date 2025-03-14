---
  title: "The Power of Synergy: Cognitive-Motor Training Boosts BDNF More Than Isolated Training in Older Adults"
author: "EG"
date: "14-03-2025"

# 1) Data set-up and management

# Load Packages
install.packages("readxl")
install.packages("plyr")
install.packages("lme4")
install.packages("nlme")
install.packages("emmeans")
install.packages("robustlmm")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("MASS")
install.packages("sjPlot")
install.packages("effectsize")

library(readxl)
library(lmerTest)
library(plyr)
library(lme4)
library(nlme)
library(emmeans)
library(robustlmm)
library(dplyr)
library(ggplot2)
library(MASS)
library(sjPlot)
library(effectsize)

# Read the dataset
Data <- read_excel("DATA-Germen.xlsx", sheet = "Sheet1")

# Check data structure
str(Data)
summary(Data)

# Calculate missing data percentage
missing_data <- colSums(is.na(Data)) / nrow(Data) * 100
print(missing_data)

# 2) Baseline Participant Characteristics by group

# Ensure Group variable is categorical
Data$Group <- as.factor(Data$Group)

# Count the number of participants by Group and Sex
sex_distribution <- Data %>%
  group_by(Group, Sex) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Display results
print(sex_distribution)


# Function to calculate Mean and SE
mean_se <- function(x) {
  n <- sum(!is.na(x))
  mean_value <- mean(x, na.rm = TRUE)
  se_value <- sd(x, na.rm = TRUE) / sqrt(n)
  return(c(Mean = mean_value, SE = se_value))
}

# Compute mean and SE for each numeric variable and show variable names
baseline_stats <- Data %>%
  group_by(Group) %>%
  summarise(
    across(
      .cols = c(Age, Education, STACC_con_trials, STACC_incon_trials, ConflictEffect, 
                STmeanRT_correct_trials, NBmeanRT_correct_target_trials, NBmeanRT_correct_nontarget_trials, 
                NBmedianRT_correct_target_trials, NBmedianRT_correct_nontarget_trials, 
                NBsdRT_correct_target_trials, NBsdRT_correct_nontarget_trials, 
                DSST_correct, DSST_error, CTSB, BDNF), 
      .fns = mean_se, 
      .names = "{.col}_mean_se"
    )
  )

# Display results
print(baseline_stats)

# 3) Reshape the data to bring pre and post values together
Data_wide <- Data %>%
  pivot_wider(names_from = Time, values_from = c(BDNF, CTSB, Age, Education, 
                                                 STACC_con_trials, STACC_incon_trials, 
                                                 ConflictEffect, STmeanRT_correct_trials, 
                                                 NBmeanRT_correct_target_trials, 
                                                 NBmeanRT_correct_nontarget_trials, 
                                                 NBmedianRT_correct_target_trials, 
                                                 NBmedianRT_correct_nontarget_trials, 
                                                 NBsdRT_correct_target_trials, 
                                                 NBsdRT_correct_nontarget_trials, 
                                                 DSST_correct, DSST_error))

# Compute delta values (Post - Pre)
Data_wide <- Data_wide %>%
  mutate(
    DeltaBDNF = BDNF_2 - BDNF_1,
    DeltaCTSB = CTSB_2 - CTSB_1,
    DeltaSTACC_con = STACC_con_trials_2 - STACC_con_trials_1,
    DeltaSTACC_incon = STACC_incon_trials_2 - STACC_incon_trials_1,
    DeltaConflictEffect = ConflictEffect_2 - ConflictEffect_1,
    DeltaSTmeanRT = STmeanRT_correct_trials_2 - STmeanRT_correct_trials_1,
    DeltaNBmeanRT_target = NBmeanRT_correct_target_trials_2 - NBmeanRT_correct_target_trials_1,
    DeltaNBmeanRT_nontarget = NBmeanRT_correct_nontarget_trials_2 - NBmeanRT_correct_nontarget_trials_1,
    DeltaNBmedianRT_target = NBmedianRT_correct_target_trials_2 - NBmedianRT_correct_target_trials_1,
    DeltaNBmedianRT_nontarget = NBmedianRT_correct_nontarget_trials_2 - NBmedianRT_correct_nontarget_trials_1,
    DeltaNBsdRT_target = NBsdRT_correct_target_trials_2 - NBsdRT_correct_target_trials_1,
    DeltaNBsdRT_nontarget = NBsdRT_correct_nontarget_trials_2 - NBsdRT_correct_nontarget_trials_1,
    DeltaDSST_correct = DSST_correct_2 - DSST_correct_1,
    DeltaDSST_error = DSST_error_2 - DSST_error_1
  )

# Print the updated dataset with delta variables
print(Data_wide)

# 4) Linear Mixed Model for Myokines and Cognitive Outcomes

# 4.1 Function to run LMEM for Myokines (without Education)
run_lme_myokine <- function(outcome) {
  model_formula <- as.formula(paste(outcome, "~ Group * Time + Age + Sex + (1 | ID)"))
  model <- lmer(model_formula, data = Data, REML = FALSE)
  print(summary(model))
  print(anova(model, type = "III"))
  print(emmeans(model, pairwise ~ Time | Group, adjust = "bonferroni"))
  return(model)
}

# 4.2. Function to run LMEM for Cognitive Measures (with Education)
run_lme_cognitive <- function(outcome) {
  model_formula <- as.formula(paste(outcome, "~ Group * Time + Age + Sex + Education + (1 | ID)"))
  model <- lmer(model_formula, data = Data, REML = FALSE)
  print(summary(model))
  print(anova(model, type = "III"))
  print(emmeans(model, pairwise ~ Time | Group, adjust = "bonferroni"))
  return(model)
}

# Run LMEM for Myokines 
model_BDNF <- run_lme_myokine("BDNF")
model_CTSB <- run_lme_myokine("CTSB")

# Run LMEM for Cognitive Measures 
model_DSST <- run_lme_cognitive("DSST_correct")
model_NBACC <- run_lme_cognitive("NBACC_trials")
model_NBmeanRT <- run_lme_cognitive("NBmeanRT_correct_trials")
model_SimonEffect <- run_lme_cognitive("SimonEffect")
model_ConflictEffect <- run_lme_cognitive("ConflictEffect")

# 4.3. Function to check assumptions after LMEM
check_lme_assumptions <- function(model) {
  # 1) Linearity & Homoscedasticity: Residual vs. Fitted plot
  plot(model, resid(.) ~ fitted(.), main = "Residuals vs. Fitted") 
  
  # 2) Normality of Residuals: QQ-plot
  qqnorm(resid(model))
  qqline(resid(model), col = "red", lwd = 2)
  
  # 3) Homoscedasticity: Scale-location plot
  plot(model, sqrt(abs(resid(.))) ~ fitted(.), main = "Scale-Location")
  
  # 4) Histogram of Residuals
  hist(resid(model), main = "Histogram of Residuals", col = "lightblue", breaks = 20)
  
  # 5) Random Effects Normality Check
  ranef_data <- ranef(model)$ID
  qqnorm(ranef_data[,1], main = "QQ Plot for Random Effects")
  qqline(ranef_data[,1], col = "red", lwd = 2)
  
  # 6) Multicollinearity Check (VIF)
  library(car)
  vif_values <- vif(model)
  print("Variance Inflation Factors (VIF):")
  print(vif_values)
  
}
# 4.4. Apply assumption checks after each LMEM
check_lme_assumptions(model_BDNF)
check_lme_assumptions(model_CTSB)
check_lme_assumptions(model_DSST)
check_lme_assumptions(model_NBACC)
check_lme_assumptions(model_NBmeanRT)
check_lme_assumptions(model_SimonEffect)
check_lme_assumptions(model_ConflictEffect)

# 5) Visualization 
# 5.1. Boxplots for BDNF and CTSB

# Function to create boxplots
plot_boxplot <- function(data, outcome, y_label, y_limit) {
  cognitive_motor_x <- which(levels(data$Group) == "Cognitive + Motor")
  
  ggplot(data, aes(x = Group, y = !!sym(outcome), fill = Time)) +
    geom_boxplot(width = 0.5, position = position_dodge(width = 0.5), outlier.shape = NA, coef = 0.6, lwd = 0.05) +
    scale_fill_manual(values = c("Pre" = "gray", "Post" = "pink")) +  
    scale_y_continuous(limits = c(0, y_limit), breaks = seq(0, y_limit, by = 2500)) +
    geom_segment(aes(x = cognitive_motor_x - 0.2, xend = cognitive_motor_x + 0.2, 
                     y = y_limit - 1000, yend = y_limit - 1000), size = 0.5, color = "black") +
    annotate("text", x = cognitive_motor_x, y = y_limit - 500, label = "* p = 0.02", size = 5) +
    labs(x = NULL, y = y_label, fill = "Time") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3, color = "black"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.ticks.y = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.3, "cm"),
      plot.margin = margin(1, 1, 1, 1),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "cm"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14)
    )
}

# Generate boxplots
plot_boxplot(Data_Germen, "BDNF", "BDNF (pg/ml)", 13500)
plot_boxplot(Data_Germen, "CTSB", "CTSB (pg/ml)", 5000)

# 5.2. Line plots for cognitive outcomes


# List of cognitive test variables
cognitive_tests <- c("SimonEffect", "DSST_correct", "NBmeanRT_correct_trials", "ConflictEffect", "NBmeanRT_correct_trials")

# Function to create a line plot for each cognitive test
plot_cognitive_test <- function(test_variable, y_label) {
  Data_Cognitive <- Data %>%
    filter(!is.na(!!sym(test_variable))) %>%
    mutate(Group = factor(Group, levels = c("Cognitive", "Motor", "Cognitive + Motor")),
           Time = factor(Time, levels = c("Pre", "Post")))
  
  ggplot(Data_Cognitive, aes(x = Time, y = !!sym(test_variable), group = Group, color = Group)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Cognitive" = "#FFB6C1", "Motor" = "gray", "Cognitive + Motor" = "lightblue")) +
    labs(x = "Pre-Intervention                Post-Intervention", y = y_label, color = "Group") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.8, color = "black"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold")
    )
}

# Generate plots for all cognitive tests
plots <- lapply(cognitive_tests, function(test) {
  plot_cognitive_test(test, paste(test, "(ms)"))  # Customize y-label if needed
})

# Display all plots in a grid (Optional)
wrap_plots(plots)


# 6) Correlation Analysis

# Function to test normality for each cognitive outcome
normality_test <- function(variable) {
  print(paste("Shapiro-Wilk normality test for", variable))
  for (group in 0:2) {
    cat("\nGroup:", group, "\n")
    print(shapiro.test(subset(Data, Group == group)[[variable]]))
  }
}

# Run normality test for all cognitive variables
variables <- c("DeltaBDNF", "DeltaDSST", "DeltaCTSB", "DeltaNBAcc", "DeltaNBRT", "DeltaSE", "DeltaCE")
for (var in variables) {
  normality_test(var)
}

# Function to compute correlation
correlation_test <- function(x, y, method) {
  if (sum(!is.na(x) & !is.na(y)) > 2) { # At least 3 valid data points required
    test <- cor.test(x, y, use = "pairwise.complete.obs", method = method)
    return(c(correlation = test$estimate, p_value = test$p.value))
  } else {
    return(c(correlation = NA, p_value = NA))
  }
}

# Compute Pearson and Spearman correlations by group
cor_results <- Data %>%
  group_by(Group) %>%
  summarise(
    BDNF_DSST_pearson_cor = correlation_test(DeltaBDNF, DeltaDSST, "pearson")["correlation"],
    BDNF_DSST_pearson_p = correlation_test(DeltaBDNF, DeltaDSST, "pearson")["p_value"],
    BDNF_NBAcc_pearson_cor = correlation_test(DeltaBDNF, DeltaNBAcc, "spearman")["correlation"],
    BDNF_NBAcc_pearson_p = correlation_test(DeltaBDNF, DeltaNBAcc, "spearman")["p_value"],
    BDNF_NBRT_pearson_cor = correlation_test(DeltaBDNF, DeltaNBRT, "spearman")["correlation"],
    BDNF_NBRT_pearson_p = correlation_test(DeltaBDNF, DeltaNBRT, "spearman")["p_value"],
    BDNF_SE_pearson_cor = correlation_test(DeltaBDNF, DeltaSE, "spearman")["correlation"],
    BDNF_SE_pearson_p = correlation_test(DeltaBDNF, DeltaSE, "spearman")["p_value"],
    BDNF_CE_spearman_cor = correlation_test(DeltaBDNF, DeltaCE, "spearman")["correlation"],
    BDNF_CE_spearman_p = correlation_test(DeltaBDNF, DeltaCE, "spearman")["p_value"],
    CTSB_DSST_spearman_cor = correlation_test(DeltaCTSB, DeltaDSST, "spearman")["correlation"],
    CTSB_DSST_spearman_p = correlation_test(DeltaCTSB, DeltaDSST, "spearman")["p_value"],
    CTSB_CE_spearman_cor = correlation_test(DeltaCTSB, DeltaCE, "spearman")["correlation"],
    CTSB_CE_spearman_p = correlation_test(DeltaCTSB, DeltaCE, "spearman")["p_value"],
    CTSB_NBAcc_spearman_cor = correlation_test(DeltaCTSB, DeltaNBAcc, "spearman")["correlation"],
    CTSB_NBAcc_spearman_p = correlation_test(DeltaCTSB, DeltaNBAcc, "spearman")["p_value"],
    CTSB_NBRT_spearman_cor = correlation_test(DeltaCTSB, DeltaNBRT, "spearman")["correlation"],
    CTSB_NBRT_spearman_p = correlation_test(DeltaCTSB, DeltaNBRT, "spearman")["p_value"],
    CTSB_SE_spearman_cor = correlation_test(DeltaCTSB, DeltaSE, "spearman")["correlation"],
    CTSB_SE_spearman_p = correlation_test(DeltaCTSB, DeltaSE, "spearman")["p_value"],
    .groups = "drop"
  )

# Display correlation results
print(cor_results)

# 7) Moderation Analysis

# Define function for moderation analysis
run_moderation_analysis <- function(outcome_var) {
  formula <- as.formula(paste(outcome_var, "~ DeltaBDNF_centered * Group"))
  model <- rlm(formula, data = DATA_Germen)
  
  cat("\nModeration Analysis for:", outcome_var, "\n")
  print(summary(model))
  
  return(model)
}

# List of cognitive outcomes
cognitive_outcomes <- c("DeltaCE", "DeltaSE", 
                        "DeltaNBAcc", "DeltaNBRT", "DeltaDSST")

# Run the models for all outcomes
moderation_models <- lapply(cognitive_outcomes, run_moderation_analysis)
