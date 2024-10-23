
# Figure 2 R

library(survival)
library(survminer)
library(ggprism)
######### Survival analysis #####################
cli.1<-read.csv("source data 1.csv")

survfit(data =cli.1, Surv(PFS, PFS_status) ~ treatment) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c( "#2c7fb8","#e6550d", "#fdc086","#A0CEE8"),
             legend = "right",
             legend.title = "",
             legend.labs = c("chemo", "ICI-chemo"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))




survfit(data =cli.1, Surv(OS, OS_status) ~ treatment) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c( "#2c7fb8","#e6550d", "#fdc086","#A0CEE8"),
             legend = "right",
             legend.title = "",
             legend.labs = c("chemo", "ICI-chemo"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))


###### Group calculation of ctDNA positive vs ctDNA negative under the condition of treatment HR value (PFS)


cli.positive<- cli.1 %>%
  filter(ctDNA_positive_C1D1=="yes")
cox_fit <- coxph(Surv(PFS, PFS_status) ~ treatment, data = cli.positive)
summary_cox<- summary(cox_fit)

# Extract HR, 95% CI, and p-values for entire dataset "cli.jco"
hr_cli <- signif(summary_cox$coefficients[1, 2], 3)
lower_cli <- signif(summary_cox$conf.int[1, 3], 2)
upper_cli <- round(summary_cox$conf.int[1, 4], 2)
p_value_cli <- signif(summary_cox$coefficients[1, 5], 2)

# Step 4: Merge all the results into a single DataFrame
df.positive <- data.frame(
  Group = c("ctDNA postive (N = 299)"),
  HR = hr_cli,
  Lower_CI = lower_cli,
  Upper_CI = upper_cli,
  P_Value = p_value_cli)


cli.negative<- cli.1 %>%
  filter(ctDNA_positive_C1D1=="no")
cox_fit <- coxph(Surv(PFS, PFS_status) ~ treatment, data = cli.negative)
summary_cox<- summary(cox_fit)

# Extract HR, 95% CI, and p-values for entire dataset "cli.jco"
hr_cli <- signif(summary_cox$coefficients[1, 2], 3)
lower_cli <- signif(summary_cox$conf.int[1, 3], 2)
upper_cli <- round(summary_cox$conf.int[1, 4], 2)
p_value_cli <- signif(summary_cox$coefficients[1, 5], 2)

# Step 4: Merge all the results into a single DataFrame
df.negative <- data.frame(
  Group = c("ctDNA negative (N = 94)"),
  HR = hr_cli,
  Lower_CI = lower_cli,
  Upper_CI = upper_cli,
  P_Value = p_value_cli)

df<-rbind(df.positive,df.negative)
df

###### Group calculation of ctDNA positive vs ctDNA negative under the condition of treatment HR value (OS)

cli.positive<- cli.1 %>%
  filter(ctDNA_positive_C1D1=="yes")
cox_fit <- coxph(Surv(OS, OS_status) ~ treatment, data = cli.positive)
summary_cox<- summary(cox_fit)

# Extract HR, 95% CI, and p-values for entire dataset "cli.jco"
hr_cli <- signif(summary_cox$coefficients[1, 2], 3)
lower_cli <- signif(summary_cox$conf.int[1, 3], 2)
upper_cli <- round(summary_cox$conf.int[1, 4], 2)
p_value_cli <- signif(summary_cox$coefficients[1, 5], 2)

#Merge all the results into a single DataFrame
df.positive <- data.frame(
  Group = c("ctDNA postive (N = 299)"),
  HR = hr_cli,
  Lower_CI = lower_cli,
  Upper_CI = upper_cli,
  P_Value = p_value_cli)


cli.negative<- cli.1 %>%
  filter(ctDNA_positive_C1D1=="no")
cox_fit <- coxph(Surv(OS, OS_status) ~ treatment, data = cli.negative)
summary_cox<- summary(cox_fit)

# Extract HR, 95% CI, and p-values for entire dataset 
hr_cli <- signif(summary_cox$coefficients[1, 2], 3)
lower_cli <- signif(summary_cox$conf.int[1, 3], 2)
upper_cli <- round(summary_cox$conf.int[1, 4], 2)
p_value_cli <- signif(summary_cox$coefficients[1, 5], 2)

# Merge all the results into a single DataFrame
df.negative <- data.frame(
  Group = c("ctDNA negative (N = 94)"),
  HR = hr_cli,
  Lower_CI = lower_cli,
  Upper_CI = upper_cli,
  P_Value = p_value_cli)

df<-rbind(df.positive,df.negative)
df
########### Create a forest plot for HR values of treatment under different ctDNA states #############


library(ggplot2)

# Define a custom theme
prism_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(family = "", size = 18))

# Plot the graph with the desired order of Variables
ggplot(df, aes(x = Group, y = HR, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_linerange(colour = "#2c7fb8", size = 3.0) +
  geom_point(colour = "#2c7fb8", size = 6.0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  xlab("") +
  ylab("") +
  prism_theme +
  scale_y_continuous(limits = c(0, 3)) +
  scale_x_discrete(limits = df$Group) +
  theme_minimal() + # Using the default font from theme_minimal()
  geom_text(aes(label = paste0("HR (95% CI)\n", round(HR, 2), " (", round(Lower_CI, 2), "-", round(Upper_CI, 2), ") p =", signif(P_Value, 2))),
            size = 5, y = 1.5, hjust = -0.1) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16),  # Customize axis labels font size
        axis.text.y = element_text(size = 16),  # Customize y-axis labels font size
        axis.ticks.x = element_line(colour = "black"),
        axis.title = element_text(size = 16),  # Customize axis titles font size
        axis.title.y = element_blank(), # Remove y-axis title
        panel.grid = element_blank()) # Remove background grid lines


###############Calculate p interaction


cox_result <- function(time, event) {
  data <- cli.1 
  cox_fit <- coxph(as.formula(paste0("Surv(", time, ", ", event, ") ~ treatment + ctDNA_positive_C1D1 + treatment*ctDNA_positive_C1D1")), data = data)
  cox_summary <- summary(cox_fit)
  
  hr_treatment <- signif(cox_summary$coefficients[1, 2], 3)
  lower_treatment <- signif(cox_summary$conf.int[1, 3], 3)
  upper_treatment <- signif(cox_summary$conf.int[1, 4], 3)
  p.value_treatment <- signif(cox_summary$coefficients[1, 5], 3)
  
  
  hr_ctDNA <- signif(cox_summary$coefficients[2, 2], 3)
  lower_ctDNA <- signif(cox_summary$conf.int[2, 3], 3)
  upper_ctDNA <- signif(cox_summary$conf.int[2, 4], 3)
  p.value_ctDNA <- signif(cox_summary$coefficients[2, 5], 3)
  
  hr_interact <- signif(cox_summary$coefficients[3, 2], 3)
  lower_interact <- signif(cox_summary$conf.int[3, 3], 3)
  upper_interact <- signif(cox_summary$conf.int[3, 4], 3)
  p.value_interact <- signif(cox_summary$coefficients[3,5],3)
  
  
  return(data.frame(hr_treatment, lower_treatment, upper_treatment,
                    p.value_treatment,
                    hr_ctDNA,
                    lower_ctDNA,
                    upper_ctDNA,
                    p.value_ctDNA,
                    hr_interact,
                    lower_interact,
                    upper_interact,
                    p.value_interact))
}

results <- cox_result(time = "PFS", event = "PFS_status")
results <- cox_result(time = "OS", event = "OS_status")
