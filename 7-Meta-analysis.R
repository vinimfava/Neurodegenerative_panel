library(readr)
library(tidyverse)

# Import the results
Genewise_results <- read_delim("Genewise_results.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       col_types = cols(T1R_affected = col_number(), 
                                                        T1R_free = col_number(), P = col_number(), 
                                                        BETA = col_number(), SE = col_number()), 
                                       trim_ws = TRUE)

# FDR correction per disease
Genewise_results <- Genewise_results %>%
  group_by(Disease, MAF) %>%
  mutate(P_adjusted = p.adjust(P, method = "BH")) %>%
  ungroup()

# calculate the OR
Genewise_results <- Genewise_results %>% 
  mutate(OR = exp(BETA)) %>%
  mutate(L95 = exp(BETA - 1.96 * SE),
         U95 = exp(BETA + 1.96 * SE))

# Inverse of the variance meta-analysis code
META_P <- function(genes, MAF = 0.01) {
# Input the beta values and standard errors
  data <- Genewise_results %>% filter(Gene %in% genes, MAF == !!MAF) %>% as.data.frame()

  # Calculate weights (inverse variance)
  data$weight <- 1 / (data$SE^2)
  
  # Meta-analysis beta (weighted mean)
  meta_beta <- sum(data$BETA * data$weight) / sum(data$weight)
  
  # Meta-analysis standard error
  meta_se <- sqrt(1 / sum(data$weight))
  
  # Z-score
  z_score <- meta_beta / meta_se
  
  # P-value
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  # Odds ratio
  OR <- exp(meta_beta)
  L95 <- exp(meta_beta - 1.96 * meta_se)
  U95 = exp(meta_beta + 1.96 * meta_se)
  
  # Results
  cat("Meta-analysis Results:\n")
  cat("Beta:", meta_beta, "\n")
  cat("Standard Error:", meta_se, "\n")
  cat("Odds ratio", OR,"(",L95," - ",U95,")", "\n")
  cat("Z-score:", z_score, "\n")
  cat("P-value:", p_value, "\n")
}

META_P(genes = c("PRKN","PINK1"))
META_P(genes = c("PRKN","PINK1","ADORA1"))
META_P(genes = c("PRKN","PINK1","TBK1"))
META_P(genes = c("PRKN","PINK1","TBK1","ADORA1"))

META_P(genes = c("LRRK2","GAK"), MAF = 0.05)
