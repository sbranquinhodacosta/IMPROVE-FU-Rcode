# This script contains the code used to analyse the data for the manuscript 
# " INSERT TITLE " 

data <- SR1_eventCIN2p_UE_excl
data <- data %>%
  mutate(CIN2_total = ifelse(worst_Histo_Total == "CIN2", 1, 0))
data <- data %>%
  mutate(CIN3_total = ifelse(worst_Histo_Total == "CIN3", 1, 0))

############################################################################################################################################################
# Table 1. 
library(dplyr)
library(gtsummary)
data
data2 <- data %>%
  mutate(
    ResHPV1 = factor(ResHPV1, levels = c("HPV negative", "HPV positive")),
    CIN2_total = factor(CIN2_total, levels = c(0,1)),
    CIN3_total = factor(CIN3_total, levels = c(0,1)),
    Cancer_total = factor(Cancer_total, levels = c(0,1)),
    SR2_eligible = factor(SR2_eligible, levels = c(0,1)),
    screen_round2 = factor(screen_round2, levels = c(0,1)),
    screen_round2_HPV = factor(screen_round2_HPV, levels = c(0,1))
    )

table1 <- data2 %>%
  select(
    randomiz_baseline,
    age_randomiz,
    ResHPV1,
    CIN2_total,
    CIN3_total,
    Cancer_total,
    SR2_eligible,
    screen_round2,
    screen_round2_HPV
  ) %>%
  tbl_summary(
    by = randomiz_baseline,
    
    statistic = list(
      age_randomiz ~ "{median} ({min}, {max})"),
    
    digits = age_randomiz ~ 1,
    
    label = list(
      age_randomiz ~ "Age at randomization",
      ResHPV1 ~ "HPV result",
      CIN2_total ~ "CIN2+",
      CIN3_total ~ "CIN3+",
      Cancer_total ~ "Cancer",
      SR2_eligible ~ "Eligible for screening round 2",
      screen_round2 ~ "Attended screening round 2",
      screen_round2_HPV ~ "HPV test at round 2")
  ) %>%
  add_overall() %>%       # adds Total column
  modify_header(label ~ "**Variable**") %>%
  bold_labels(); table1

############################################################################################################################################################
#Table 3. Two-round Cumulative risks and risk difference according to baseline randomisation with subgroups
# The example code concerns outcome "CIN3+" and can be used for "CIN2+" by changing the outcome-related variables in the code.
# The K-M Plots correspond to what is presented in Figure 2.

# Create the new variable 'timetoCIN3p'
data <- data %>%
  mutate(timetoCIN3p = case_when(
    CIN3plus_total == 1 | Hysterectomy == 1 ~ CIN3plus_time,  # First condition
    CIN3plus_total == 0 & CIN3plus_time == 4.0 ~ 4.0,        # Second condition
    CIN3plus_total == 0 & CIN3plus_time > 4.0 ~ 8.0,         # Third condition
    TRUE ~ CIN3plus_time                                      # Default case if none of the above conditions match
  )) %>%
  mutate(timetoCIN3p = round(timetoCIN3p, 1))  # Convert to character with 1 decimal place

##############################################################################
###                       Analysis - CIN3 - All women                      ###
##############################################################################
Data_analysis <- data

Data_analysis$randomiz_baseline <- as.factor(Data_analysis$randomiz_baseline)
Data_analysis$CIN3plus_time <- as.numeric(Data_analysis$CIN3plus_time)
Data_analysis$timetoCIN3p <- as.numeric(Data_analysis$timetoCIN3p)
Data_analysis$CIN3plus_total <- as.numeric(Data_analysis$CIN3plus_total)

#####
# Cumulative risk / baseline Randomization - CIN3+
library(survival)
library(dplyr)

# Survival object
surv_obj_CIN3p <- Surv(Data_analysis$timetoCIN3p,
                       Data_analysis$CIN3plus_total == 1)

# Log-rank test
survdiff_cin3 <- survdiff(
  Surv(timetoCIN3p, CIN3plus_total) ~ randomiz_baseline,
  data = Data_analysis,
  rho = 0
)

survdiff_cin3
survdiff_cin3$pvalue

# Kaplan-Meier fit
fit1_CIN3p <- survfit(
  surv_obj_CIN3p ~ randomiz_baseline,
  data = Data_analysis
)

summary(fit1_CIN3p)
summary(fit1_CIN3p, times = c(2,4,6,8))

# Cumulative incidence dataframe
cuminc_fit1_df <- data.frame(
  time = fit1_CIN3p$time,
  n.event = fit1_CIN3p$n.event,
  cuminc = (1 - fit1_CIN3p$surv) * 100,
  std_error_cuminc = fit1_CIN3p$std.err * 100,
  ci = paste0(
    round((1 - fit1_CIN3p$upper) * 100, 4),
    "-",
    round((1 - fit1_CIN3p$lower) * 100, 4)
  )
) %>%
  mutate(across(where(is.numeric), round, 4)) %>%
  mutate(
    Group = rep(names(fit1_CIN3p$strata), fit1_CIN3p$strata),
    Group = gsub("randomiz_baseline=", "", Group)
  ) %>%
  arrange(Group, time)

# Extract cumulative incidence at timepoint==8
calculate_cuminc <- function(data, time_point) {
  data %>%
    filter(time == time_point) %>%
    select(
      group = Group,
      time,
      cuminc,
      se_cuminc = std_error_cuminc,
      ci_cuminc = ci
    )
}

cuminc_2SR <- calculate_cuminc(cuminc_fit1_df, 8); cuminc_2SR  # Cumulative incidence at 8 years

#######
# Risk Difference at timepoint==8
library(dplyr)

time_point <- 8

summary_fit <- summary(fit1_CIN3p, times = time_point)

risk_df <- data.frame(
  group = gsub("randomiz_baseline=", "", summary_fit$strata),
  surv = summary_fit$surv,
  se_surv = summary_fit$std.err
) %>%
  mutate(
    cum_risk = 1 - surv
  )

risk_SS <- risk_df %>% filter(group == "SS")
risk_CS <- risk_df %>% filter(group == "CS")

# Risk difference
risk_diff <- risk_SS$cum_risk - risk_CS$cum_risk

# SE risk difference
se_diff <- sqrt(risk_SS$se_surv^2 + risk_CS$se_surv^2)

# 95% CI
lower_ci <- risk_diff - 1.96 * se_diff
upper_ci <- risk_diff + 1.96 * se_diff

risk_difference_results <- data.frame(
  time = time_point,
  risk_difference = risk_diff * 100,
  lower_ci = lower_ci * 100,
  upper_ci = upper_ci * 100
)

risk_difference_results

#######
# K-M Plots - Figure 2.

#Set colors for the lines 
custom_colors <- c("dark blue", "tomato1")

plot <- ggsurvplot(fit1_CIN3p,
                   conf.int = TRUE,
                   conf.int.alpha = 0.05,
                   conf.int.style = c("step"),
                   surv.scale="percent", 
                   risk.table = FALSE, 
                   palette = custom_colors,
                   censor = FALSE,
                   risk.table.col = "strata", # Change risk table color by groups
                   ggtheme = theme_minimal(), # Change ggplot2 theme
                   risk.table.theme = theme(box.color = "white"),
                   xlim=c(0, 8),
                   ylim=c(0.0, 0.05),
                   xlab = "Time since baseline HPV test (years)",
                   ylab = "Cumulative risk - CIN3+",
                   fun = "event",
                   risk.table.y.text = TRUE,
                   risk.table.fontsize = 3,
                   tables.y.text = FALSE,
                   legend = c(0.2,0.95),
                   legend.title="",
                   legend.labs = c("Clinician-sample", "Self-sample"),
                   legend.background = element_rect(fill='transparent'),
                   legend.box.background = element_rect(fill='transparent'),
                   font.legend = list(size = 9, color = "black"))
plot$plot + theme(plot.title = element_text(hjust=0.5))

plot$plot +
  annotate("text",
           x = -0.35, y = 0.05,
           label = "a)",
           fontface = "bold",
           size = 5,
           hjust = 0, vjust = 1) +
  theme(plot.title = element_text(hjust = 0.5))

##############################################################################
###              Analysis - CIN3 - Baseline HPV positive women             ###
##############################################################################
# Run code to subset the data
Data_analysis <- subset(data, data$ResHPV1=="HPV positive") #HPV positive women only   

# Run code lines [79 to 183] to obtain the values for Cumulative Risk and Risk Difference

#######
# K-M Plots - Figure 2.

#Set colors for the lines 
custom_colors <- c("dark blue", "tomato1")

plot <- ggsurvplot(fit1_CIN3p,
                   conf.int = TRUE,
                   conf.int.alpha = 0.05,
                   conf.int.style = c("step"),
                   surv.scale="percent", 
                   risk.table = FALSE, 
                   palette = custom_colors,
                   censor = FALSE,
                   risk.table.col = "strata", # Change risk table color by groups
                   ggtheme = theme_minimal(), # Change ggplot2 theme
                   risk.table.theme = theme(box.color = "white"),
                   xlim=c(0, 8),
                   ylim=c(0.0, 0.35),
                   xlab = "Time since baseline HPV-positive result (years)",
                   ylab = "Cumulative risk - CIN3+",
                   fun = "event",
                   risk.table.y.text = TRUE,
                   risk.table.fontsize = 3,
                   tables.y.text = FALSE,
                   legend = c(0.2,0.95),
                   legend.title="",
                   legend.labs = c("Clinician-sample", "Self-sample"),
                   legend.background = element_rect(fill='transparent'),
                   legend.box.background = element_rect(fill='transparent'),
                   font.legend = list(size = 9, color = "black"))
plot$plot + theme(plot.title = element_text(hjust=0.5))

plot$plot +
  annotate("text",
           x = -0.35, y = 0.35,
           label = "c)",
           fontface = "bold",
           size = 5,
           hjust = 0, vjust = 1) +
  theme(plot.title = element_text(hjust = 0.5))

##############################################################################
###              Analysis - CIN3 - Baseline HPV positive women             ###
##############################################################################
# Run code to subset the data
Data_analysis <- subset(data, data$ResHPV1=="HPV negative") #HPV negative women only

# Run code lines [79 to 183] to obtain the values for Cumulative Risk and Risk Difference

#######
# K-M Plots - Figure 2.

#Set colors for the lines 
custom_colors <- c("dark blue", "tomato1")

plot <- ggsurvplot(fit1_CIN3p,
                   conf.int = TRUE,
                   conf.int.alpha = 0.05,
                   conf.int.style = c("step"),
                   surv.scale="percent", 
                   risk.table = FALSE, 
                   palette = custom_colors,
                   censor = FALSE,
                   risk.table.col = "strata", # Change risk table color by groups
                   ggtheme = theme_minimal(), # Change ggplot2 theme
                   risk.table.theme = theme(box.color = "white"),
                   xlim=c(0, 8),
                   ylim=c(0.0, 0.05),
                   xlab = "Time since baseline HPV-negative result (years)",
                   ylab = "Cumulative risk - CIN3+",
                   fun = "event",
                   risk.table.y.text = TRUE,
                   risk.table.fontsize = 3,
                   tables.y.text = FALSE,
                   legend = c(0.2,0.95),
                   legend.title="",
                   legend.labs = c("Clinician-sample", "Self-sample"),
                   legend.background = element_rect(fill='transparent'),
                   legend.box.background = element_rect(fill='transparent'),
                   font.legend = list(size = 9, color = "black"))
plot$plot + theme(plot.title = element_text(hjust=0.5))

plot$plot +
  annotate("text",
           x = -0.35, y = 0.05,
           label = "e)",
           fontface = "bold",
           size = 5,
           hjust = 0, vjust = 1) +
  theme(plot.title = element_text(hjust = 0.5))

##############################################################################
###          Analysis - CIN3 - Age at baseline: Under 40 years-old         ###
##############################################################################
data <- data %>%
  mutate(age_cat = ifelse(age_randomiz<40, 1, 0))

# Run code to subset the data
Data_analysis <- subset(data, data$age_cat==1) # Women U40 yo only   

# Run code lines [79 to 183] to obtain the values for Cumulative Risk and Risk Difference

##############################################################################
###          Analysis - CIN3 - Age at baseline: 40 years & older           ###
##############################################################################

# Run code to subset the data
Data_analysis <- subset(data, data$age_cat==0) # Women O40 yo only   

# Run code lines [79 to 183] to obtain the values for Cumulative Risk and Risk Difference


############################################################################################################################################################
#Table 4. Two-round Cumulative risks and risk difference according to baseline randomisation and 
# Estimates obtained using the same code lines [79 to 183] but with data subset according to data$ResHPV1 & data$age_cat


