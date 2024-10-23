###########Figure 3 and Extended Data 5.
rm(list = ls())
gc()

cli.1<-read.csv("source.3.csv")

# Select the columns of interest and the target variables
data_subset <- cli.1 %>%
  filter(ctDNA_positive_C1D1=="no")%>%
  select(treatment, PFS, PFS_status, OS, OS_status, PDL1.bina, KEAP1, KRAS, STK11, TP53, TMB.bina, TNB.bina, wGII.bina)

str(data_subset)

data_subset<- data_subset %>%
  mutate_at(vars(tail(names(data_subset), 8)), as.numeric)
# Create an empty dataframe to store regression results
results <- data.frame(
  Variable = character(),
  Value = numeric(),
  HR_PFS = numeric(),
  CI_lower_PFS = numeric(),
  CI_upper_PFS = numeric(),
  p_value_PFS = numeric(),
  HR_OS = numeric(),
  CI_lower_OS = numeric(),
  CI_upper_OS = numeric(),
  p_value_OS = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each column of interest
for (col_name in c("PDL1.bina", "KEAP1", "KRAS", "STK11", "TP53", "TMB.bina", "TNB.bina", "wGII.bina")) {
  
  
  # Loop through each value (0 and 1)
  for (value in c(0, 1)) {
    # Create a subset containing only the data for the current column and value, converting the column to numeric
    
    subset_data <- data_subset %>%
      filter(.[,col_name] == value)
    
    
    # Perform Cox regression for PFS
    cox_model_PFS <- coxph(Surv(PFS_new, PFS_status_new) ~ treatment, data = subset_data)
    cox_summary_PFS <- summary(cox_model_PFS)
    cox_summary_PFS[["conf.int"]]
    # Extract HR, CI, and p-value for PFS
    hr_PFS <- cox_summary_PFS$coefficients["treatmentICI", "exp(coef)"]
    ci_lower_PFS <- cox_summary_PFS$conf.int["treatmentICI", "lower .95"] 
    ci_upper_PFS <- cox_summary_PFS$conf.int["treatmentICI", "upper .95"] 
    p_value_PFS <- cox_summary_PFS$coefficients["treatmentICI", "Pr(>|z|)"]
    
    # Perform Cox regression for OS
    cox_model_OS <- coxph(Surv(OS_new, OS_status_new) ~ treatment, data = subset_data)
    cox_summary_OS <- summary(cox_model_OS)
    
    # Extract HR, CI, and p-value for OS
    hr_OS <- cox_summary_OS$coefficients["treatmentICI", "exp(coef)"]
    ci_lower_OS <- cox_summary_OS$conf.int["treatmentICI", "lower .95"] 
    ci_upper_OS <- cox_summary_OS$conf.int["treatmentICI", "upper .95"] 
    p_value_OS <- cox_summary_OS$coefficients["treatmentICI", "Pr(>|z|)"]
    
    ## Add results to the dataframe
    results <- rbind(results, data.frame(
      Variable = col_name,
      Value = value,
      HR_PFS = hr_PFS,
      CI_lower_PFS = ci_lower_PFS,
      CI_upper_PFS = ci_upper_PFS,
      p_value_PFS = p_value_PFS,
      HR_OS = hr_OS,
      CI_lower_OS = ci_lower_OS,
      CI_upper_OS = ci_upper_OS,
      p_value_OS = p_value_OS
    ))
  }
  
  
}

# Display the results
print(results)
