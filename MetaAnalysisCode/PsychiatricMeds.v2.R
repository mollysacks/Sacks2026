###############################
### Meds Meta-analysis v2   ###
###.    November 2025       ###
###############################

## Setup
setwd('/path/to/your/MetaAnalysis/')
library(dplyr)
library(metafor)

cohorts_all <- c('UKBB.EUR', 'UKBB.AFR', 'UKBB.ASN', 
                 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 
                 'MyCode.EUR', 'MyCode.AFR', 'MyCode.LATNAT', 
                 'EstBB', 
                 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT')
cohorts <- c('UKBB.EUR', 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 'MyCode.EUR', 'EstBB', 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT')
cohorts_EUR <- c('UKBB.EUR', 'AoU.EUR', 'MyCode.EUR', 'EstBB', 'MVP.EUR')

dir.create('PsychiatricMedication.v2')
dir.create('PsychiatricMedication.v2/fixed_effects')
dir.create('PsychiatricMedication.v2/random_effects')
dir.create('PsychiatricMedication.v2/raw_quantiles')
dir.create('PsychiatricMedication.v2/scaled_quantiles')
dir.create('PsychiatricMedication.v2/tables')


######################
## define functions ##
######################

fixed_effect_ma <- function(input){
  fixed <- rma.uni(yi=input[which(!is.na(input$se)),'estimate'], sei=input[which(!is.na(input$se)),'se'],method="FE")
  output <- c(sum(input$count), fixed$b, fixed$se, fixed$zval, fixed$pval, fixed$ci.lb, fixed$ci.ub)
  return(output)
}

random_effect_ma <- function(input){
  random <- rma.uni(yi=input[which(!is.na(input$se)),'estimate'], sei=input[which(!is.na(input$se)),'se'],method="REML")
  output <- c(sum(input$count), random$b, random$se, random$zval, random$pval, random$ci.lb, random$ci.ub)
  return(output)
}

fixed_effect_ma_int <- function(input){
  fixed <- rma.uni(yi=input[which(!is.na(input$intercept_se)),'intercept_estimate'], sei=input[which(!is.na(input$intercept_se)),'intercept_se'],method="FE")
  output <- c(sum(input$count), fixed$b, fixed$se, fixed$zval, fixed$pval, fixed$ci.lb, fixed$ci.ub)
  return(output)
}

random_effect_ma_int <- function(input){
  random <- rma.uni(yi=input[which(!is.na(input$intercept_se)),'intercept_estimate'], sei=input[which(!is.na(input$intercept_se)),'intercept_se'],method="REML")
  output <- c(sum(input$count), random$b, random$se, random$zval, random$pval, random$ci.lb, random$ci.ub)
  return(output)
}


meta_linear_quantiles <- function(dat) {
  # initialize reconstructed distribution
  all_samples <- c()
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  # for each cohort, simulate log-normal distribution based on quantiles
  for (i in 1:nrow(dat)) {
    row <- dat[i, ]
    qvals <- as.numeric(row[c("pct025_raw","pct25_raw","pct50_raw","pct75_raw","pct975_raw")])
    # Skip invalid rows
    if (any(is.na(qvals)) || any(qvals <= 0)) {
      next
    }
    
    # Create linear approximation function
    qfun <- approxfun(probs, qvals, method = "linear", rule = 2)
    
    # Sample from the approximated CDF
    samples <- qfun(runif(row$count))
    all_samples <- c(all_samples, samples)
  }
  
  # Compute pooled quantiles
  qs <- quantile(all_samples, probs = probs)
  out <- c(length(all_samples), qs[1], qs[2], qs[3], qs[4], qs[5])
  return(out)
}

meta_normal_quantiles <- function(dat) {
  # initialize reconstructed distribution
  all_samples <- c()
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  # for each cohort, simulate log-normal distribution based on quantiles
  for (i in 1:nrow(dat)) {
    row <- dat[i, ]
    qvals <- as.numeric(row[c("pct025","pct25","pct50","pct75","pct975")])
    # Skip invalid rows
    if (any(is.na(qvals))) {
      next
    }
    
    # Estimate mean and SD for this study
    mu <- qvals[3]  # median ~ mean
    sigma <- (qvals[4] - qvals[2]) / (qnorm(0.75) - qnorm(0.25))  # IQR method
    
    # simulate n samples
    samples <- rnorm(row$count, mean = mu, sd = sigma)
    all_samples <- c(all_samples, samples)
  }
  
  # Compute pooled quantiles
  qs <- quantile(all_samples, probs = probs)
  out <- c(length(all_samples), qs[1], qs[2], qs[3], qs[4], qs[5])
  return(out)
}

meta_r2 <- function(r2, n, model) {
  # r2: vector of R² values from each study
  # n:  vector of corresponding sample sizes
  
  # check input
  if(length(r2) != length(n)) stop("r2 and n must be same length")
  
  # if R2 is negative (could be for adj R2), set to 0
  r2[r2 < 0] <- 0
  
  # transform R² to Fisher's Z via correlation r = sqrt(R²)
  r <- sqrt(r2)
  fisher_z <- atanh(r)
  
  # compute standard errors (approx for Fisher's Z)
  se <- 1 / sqrt(n - 3)
  
  # meta-analyze using fixed or random effects (here: random)
  m <- rma.uni(yi = fisher_z, sei = se, method = model)
  
  # back-transform pooled effect
  pooled_r <- tanh(m$b)
  pooled_r2 <- pooled_r^2
  
  # confidence intervals
  ci_lb_r <- tanh(m$ci.lb)
  ci_ub_r <- tanh(m$ci.ub)
  
  ci_lb_r2 <- ci_lb_r^2
  ci_ub_r2 <- ci_ub_r^2
  
  # results
  return(c(pooled_r2, ci_lb_r2, ci_ub_r2))
}


#####################################
### 1. Medication main effect.    ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "1_med_main_effect.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_raw) <- c('medication', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_scaled) <- c('medication', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this medication
  medication <- unique(table_all$medication)[i]
  input <- table_all[which(table_all$medication == medication),]
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(medication), random)
  # Raw quantiles
  raw <- meta_linear_quantiles(input)
  quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication), raw)
  # Scaled quantiles
  scaled <- meta_normal_quantiles(input)
  quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication), scaled)
  
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 2. CNV x Med                  ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "2_CNVxmed.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$locus))){
    locus <- unique(table_all$locus)[j]
    input <- table_all[which(table_all$medication == medication & table_all$locus == locus),]
    
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, locus), fixed)
    # Random effects
    # Random effects
    random <- tryCatch(
      random_effect_ma(input),
      error = function(e) {
        # fallback: fixed effects == random with tau^2 = 0
        fixed_effect_ma(input)
      }
    )
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, locus), random)
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))

######################################
### 2. CNV vs. Med Chi2            ###
######################################

noCNV_Nmapping <- c("UKBB.EUR"=305436,
                    "AoU.EUR"=212192,
                    "AoU.AAM"=78896,
                    "AoU.LATNAT"=45880,
                    "MyCode.EUR"=101064,
                    "EstBB"=58075,
                    "MVP.EUR"=418539,
                    "MVP.AFR"=112820,
                    "MVP.LATNAT"=59863)


#####################################
### 3. Medication by CNV          ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "3_med_by_CNV.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('medication', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('medication', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$locus))){
    locus <- unique(table_all$locus)[j]
    input <- table_all[which(table_all$medication == medication & table_all$locus == locus),]
    
    if(nrow(input[which(!is.na(input$se)),]) == 0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, locus), fixed)
    # Random effects
    random <- tryCatch(
      random_effect_ma(input),
      error = function(e) {
        # fallback: fixed effects == random with tau^2 = 0
        fixed_effect_ma(input)
      }
    )
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, locus), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, locus), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, locus), scaled)
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))



#####################################
### 4. CNV by med                 ###
#####################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "4_CNV_by_med.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('medication', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('medication', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$locus))){
    locus <- unique(table_all$locus)[j]
    input <- table_all[which(table_all$medication == medication & table_all$locus == locus),]
    
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, locus), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, locus), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, locus), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, locus), scaled)
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))

############################
### 5. Med x PRS         ###
############################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "5_MedXPRS.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  input <- table_all[which(table_all$medication == medication),]
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(medication), random)
  
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
          
#####################################
### 6. Medication by PRS.         ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "6_Med_by_PRS.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('medication', 'quartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('medication', 'quartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('medication', 'quartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('medication', 'quartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$quartile))){
    quartile <- unique(table_all$quartile)[j]
    input <- table_all[which(table_all$medication == medication & table_all$quartile == quartile),]
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }

    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, quartile), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, quartile), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, quartile), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, quartile), scaled)
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 7. PRS by Medication          ###
#####################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "7_PRS_by_Med.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=15))
colnames(fixed_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                             'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=15))
colnames(random_effects) <- c('medication', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper',
                              'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  input <- table_all[which(table_all$medication == medication),]
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_int <- fixed_effect_ma_int(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication), fixed, fixed_int)
  # Random effects
  random <- random_effect_ma(input)
  random_int <- random_effect_ma_int(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(medication), random, random_int)
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))

#####################################
### 8. Medication x sex           ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "8_MedXSex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('medication', 'count_males', 'count_females',  'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('medication', 'count_males', 'count_females', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  input <- table_all[which(table_all$medication == medication),]
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  count_males <- sum(input$count_males)
  count_females <- sum(input$count_females)
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, count_males, count_females), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(medication, count_males, count_females), random)
  
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))


#####################################
### 9. Medication by sex          ###
#####################################



# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "9_Med_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c("X", "medication", "sex", "count", "estimate", "se", "p", "mean_raw", "sd_raw",   
                         "pct025_raw", "pct25_raw", "pct50_raw", "pct75_raw", "pct975_raw", "mean", "sd",
                         "pct025", "pct25", "pct50", "pct75", "pct975", "cohort" )

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$sex))){
    sex <- unique(table_all$sex)[j]
    input <- table_all[which(table_all$medication == medication & table_all$sex == sex),]
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, sex), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, sex), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, sex), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, sex), scaled)
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 10. Medication by CNV by sex  ###
#####################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "10_Med_by_CNV_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('locus', 'medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('locus', 'medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('locus', 'medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('locus', 'medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$sex))){
    sex <- unique(table_all$sex)[j]
    for(k in seq_along(unique(table_all$locus))){
      locus <- unique(table_all$locus)[k]
      input <- table_all[which(table_all$medication == medication & table_all$sex == sex),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, medication, sex), fixed)
      # Random effects
      random <- random_effect_ma(input)
      random_effects[nrow(random_effects) + 1,] <- c(c(locus, medication, sex), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, medication, sex), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, medication, sex), scaled)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 11. CNV by med by sex.        ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "11_CNV_by_Med_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('locus', 'medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('locus', 'medication', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('locus', 'medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('locus', 'medication', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$sex))){
    sex <- unique(table_all$sex)[j]
    for(k in seq_along(unique(table_all$locus))){
      locus <- unique(table_all$locus)[k]
      input <- table_all[which(table_all$medication == medication & table_all$sex == sex & table_all$locus == locus),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, medication, sex), fixed)
      # Random effects
      random <- tryCatch(
        random_effect_ma(input),
        error = function(e) {
          # fallback: fixed effects == random with tau^2 = 0
          fixed_effect_ma(input)
        }
      )
      random_effects[nrow(random_effects) + 1,] <- c(c(locus, medication, sex), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, medication, sex), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, medication, sex), scaled)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 12. PRS by CNV by med         ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "12_PRS_by_CNV_by_Med.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=16))
colnames(fixed_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                             'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=16))
colnames(random_effects) <- c('medication', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper',
                              'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(k in seq_along(unique(table_all$locus))){
    locus <- unique(table_all$locus)[k]
    input <- table_all[which(table_all$medication == medication & table_all$locus == locus),]
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_int <- fixed_effect_ma_int(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, locus), fixed, fixed_int)
    # Random effects
    random <- random_effect_ma(input)
    random_int <- random_effect_ma_int(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, locus), random, random_int)
    
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))


#####################################
### 13. Medication by PRS by CNV  ###
#####################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "13_Meds_by_PRS_by_CNV.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('medication', 'quartile', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('medication', 'quartile', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('medication', 'quartile', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('medication', 'quartile', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$quartile))){
    quartile <- unique(table_all$quartile)[j]
    for(k in seq_along(unique(table_all$locus))){
      locus <- unique(table_all$locus)[k]
      input <- table_all[which(table_all$medication == medication & table_all$quartile == quartile & table_all$locus == locus),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, quartile, locus), fixed)
      # Random effects
      random <- random_effect_ma(input)
      random_effects[nrow(random_effects) + 1,] <- c(c(medication, quartile, locus), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, quartile, locus), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, quartile, locus), scaled)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 14. CNV by PRS by meds        ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "14_CNV_by_PRS_by_meds.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('medication', 'quartile', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('medication', 'quartile', 'locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('medication', 'quartile', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('medication', 'quartile', 'locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$quartile))){
    quartile <- unique(table_all$quartile)[j]
    for(k in seq_along(unique(table_all$locus))){
      locus <- unique(table_all$locus)[k]
      input <- table_all[which(table_all$medication == medication & table_all$quartile == quartile & table_all$locus == locus),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, quartile, locus), fixed)
      # Random effects
      random <- tryCatch(
        random_effect_ma(input),
        error = function(e) {
          # fallback: fixed effects == random with tau^2 = 0
          fixed_effect_ma(input)
        }
      )
      random_effects[nrow(random_effects) + 1,] <- c(c(medication, quartile, locus), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(medication, quartile, locus), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(medication, quartile, locus), scaled)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('PsychiatricMedication.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('PsychiatricMedication.v2/scaled_quantiles/', out_string))


#####################################
### 15. PRS by CNV by med by sex  ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "15_PRS_by_CNV_by_med_by_sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=17))
colnames(fixed_effects) <- c('medication', 'locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                             'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=17))
colnames(random_effects) <- c('medication', 'locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper',
                              'count2', 'intercept_estimate', 'intercept_se', 'intercept_z', 'intercept_p', 'intercept_ci_lower', 'intercept_ci_upper')

for(i in seq_along(unique(table_all$medication))){
  # limit to this locus
  medication <- unique(table_all$medication)[i]
  for(j in seq_along(unique(table_all$sex))){
    # limit to this locus
    sex <- unique(table_all$sex)[j]
    for(k in seq_along(unique(table_all$locus))){
      locus <- unique(table_all$locus)[k]
      input <- table_all[which(table_all$medication == medication & table_all$locus == locus & table_all$sex == sex),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_int <- fixed_effect_ma_int(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, locus, sex), fixed, fixed_int)
      # Random effects
      random <- tryCatch(
        random_effect_ma(input),
        error = function(e) {
          # fallback: fixed effects == random with tau^2 = 0
          fixed_effect_ma(input)
        }
      )
      random_int <- tryCatch(
        random_effect_ma_int(input),
        error = function(e) {
          # fallback: fixed effects == random with tau^2 = 0
          fixed_effect_ma_int(input)
        }
      )
      random_effects[nrow(random_effects) + 1,] <- c(c(medication, locus, sex), random, random_int)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))


#########################################
### 16. CNV main effect meds factor   ###
#########################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "16_CNV_main_effect_meds_factor.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(j in seq_along(unique(table_all$locus))){
  locus <- unique(table_all$locus)[j]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))

#####################################
### 17.  Variance explained       ###
#####################################


anova_all <- data.frame(matrix(nrow=0, ncol=6))
colnames(anova_all) <- c("model", "Rsq", "Rsq_adj", "n", "m", "cohort")

anova_fixed <- data.frame(matrix(nrow=0, ncol=7))
colnames(anova_fixed) <- c("model", "Rsq", "ci_lower", "ci_upper", "Rsq_adj", "ci_lower_adj", "ci_upper_adj")
anova_random <- data.frame(matrix(nrow=0, ncol=7))
colnames(anova_random) <- c("model", "Rsq", "ci_lower", "ci_upper", "Rsq_adj", "ci_lower_adj", "ci_upper_adj")


models <- c("CNV_meds_ANOVA.csv",
            "CNV_PRS_meds_additive_ANOVA.csv", 
            "CNV_PRS_meds_interaction_ANOVA.csv",
            "covariates_meds_ANOVA.csv",
            "PRS_meds_ANOVA.csv")
models_EUR <- c("CNV_PRS_meds_additive_ANOVA.csv", 
                "CNV_PRS_meds_interaction_ANOVA.csv",
                "PRS_meds_ANOVA.csv")


for(j in seq_along(models)){
  model <- models[j]
  model_name <- substr(model, 1, nchar(model) - 10)
  input <- data.frame(matrix(nrow=0, ncol=6))
  colnames(input) <- c("model", "Rsq", "Rsq_adj", "n", "m", "cohort")
  
  for(i in seq_along(cohorts)){
    cohort <- cohorts[i]
    if(!(cohort %in% cohorts_EUR) & (model %in% models_EUR)){
      next
    }
    anova_df <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/17_variance_explained/', model))
    rownames(anova_df) <- anova_df$X
    Rsq <- 1 - (anova_df["Residuals", "Sum.Sq"]/sum(anova_df$Sum.Sq))
    Rsq_adj <- 1 - ((1 - Rsq) * (sum(anova_df$Df)/anova_df["Residuals", "Df"]))
    n <-  1 + sum(anova_df$Df) # sample size
    m <- sum(anova_df$Df) - anova_df["Residuals", "Df"] # number of predictors
    anova_all[nrow(anova_all) + 1,] <- c(model_name, Rsq, Rsq_adj, n, m, cohort)
    input[nrow(input) + 1,] <- c(model_name, Rsq, Rsq_adj, n, m, cohort)
  }
  
  input$Rsq <- as.numeric(input$Rsq)
  input$Rsq_adj <- as.numeric(input$Rsq_adj)
  input$n <- as.numeric(input$n)
  input$m <- as.numeric(input$m)
  
  # apply Fisher's transformation
  input$ZR2 <- 0.5 * log((1 + sqrt(input$Rsq)) / (1 - sqrt(input$Rsq)))
  input$ZR2_adj <- 0.5 * log((1 + sqrt(input$Rsq_adj)) / (1 - sqrt(input$Rsq_adj)))
  
  fixed <- rma.uni(yi = input$ZR2, sei = 1/sqrt(input$n - 3), measure = "ZCOR", method="FE")
  ma_R2 <- ((exp(2 * fixed$b) - 1) / (exp(2 * fixed$b) + 1)) ** 2
  ma_R2_ci_lower <- ((exp(2 * fixed$ci.lb) - 1) / (exp(2 * fixed$ci.lb) + 1)) ** 2
  ma_R2_ci_upper <- ((exp(2 * fixed$ci.ub) - 1) / (exp(2 * fixed$ci.ub) + 1)) ** 2
  
  fixed_adj <- rma.uni(yi = input$ZR2_adj, sei = 1/sqrt(input$n - 3), measure = "ZCOR", method="FE")
  ma_R2_adj <- ((exp(2 * fixed_adj$b) - 1) / (exp(2 * fixed_adj$b) + 1)) ** 2
  ma_R2_ci_lower_adj <- ((exp(2 * fixed_adj$ci.lb) - 1) / (exp(2 * fixed_adj$ci.lb) + 1)) ** 2
  ma_R2_ci_upper_adj <- ((exp(2 * fixed_adj$ci.ub) - 1) / (exp(2 * fixed_adj$ci.ub) + 1)) ** 2
  
  anova_fixed[nrow(anova_fixed) + 1,] <- c(model_name, ma_R2, ma_R2_ci_lower, ma_R2_ci_upper, ma_R2_adj, ma_R2_ci_lower_adj, ma_R2_ci_upper_adj)
  
  random <- rma.uni(yi = input$ZR2, sei = 1/sqrt(input$n - 3), measure = "ZCOR", method="REML")
  ma_R2 <- ((exp(2 * random$b) - 1) / (exp(2 * random$b) + 1)) ** 2
  ma_R2_ci_lower <- ((exp(2 * random$ci.lb) - 1) / (exp(2 * random$ci.lb) + 1)) ** 2
  ma_R2_ci_upper <- ((exp(2 * random$ci.ub) - 1) / (exp(2 * random$ci.ub) + 1)) ** 2
  
  random_adj <- rma.uni(yi = input$ZR2_adj, sei = 1/sqrt(input$n - 3), measure = "ZCOR", method="REML")
  ma_R2_adj <- ((exp(2 * random_adj$b) - 1) / (exp(2 * random_adj$b) + 1)) ** 2
  ma_R2_ci_lower_adj <- ((exp(2 * random_adj$ci.lb) - 1) / (exp(2 * random_adj$ci.lb) + 1)) ** 2
  ma_R2_ci_upper_adj <- ((exp(2 * random_adj$ci.ub) - 1) / (exp(2 * random_adj$ci.ub) + 1)) ** 2
  
  anova_random[nrow(anova_random) + 1,] <- c(model_name, ma_R2, ma_R2_ci_lower, ma_R2_ci_upper, ma_R2_adj, ma_R2_ci_lower_adj, ma_R2_ci_upper_adj)
}

write.csv(anova_all, 'PsychiatricMedication.v2/tables/17_variance_explained.csv')
write.csv(anova_fixed, 'PsychiatricMedication.v2/fixed_effects/17_variance_explained.csv')
write.csv(anova_random, 'PsychiatricMedication.v2/random_effects/17_variance_explained.csv')

#####################################
### 18. CNV causal mediation.     ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "18_CNV_causal_mediation.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

table_all$se <- (table_all$ci_upper - table_all$ci_lower) / (1.96*2)

fixed_effects <- data.frame(matrix(nrow=0, ncol=14))
colnames(fixed_effects) <- c('CNV', 'medication', 'Effect',
                             "nCNV_med", "nCNV_nomed", "n_noCNV_med", "n_noCNV_nomed",
                             'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=14))
colnames(random_effects) <- c('CNV', 'medication', 'Effect',
                              "nCNV_med", "nCNV_nomed", "n_noCNV_med", "n_noCNV_nomed",
                              'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$CNV))){
  CNV <- unique(table_all$CNV)[i]
  for(j in seq_along(unique(table_all$medication))){
    medication <- unique(table_all$medication)[j]
    for(k in seq_along(unique(table_all$Effect))){
      Effect <- unique(table_all$Effect)[k]
      input <- table_all[which((table_all$CNV == CNV) & (table_all$medication == medication) & (table_all$Effect == Effect)),]
      input <- input[which(input$se != 0),]
      if(nrow(input[which(!is.na(input$se)),]) ==0){
        next
      }
      
      nCNV_med <- sum(input$nCNV_med)
      nCNV_nomed <- sum(input$nCNV_nomed)
      n_noCNV_med <- sum(input$n_noCNV_med)
      n_noCNV_nomed <- sum(input$n_noCNV_nomed)
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(CNV, medication, Effect, 
                                                     nCNV_med, nCNV_nomed, n_noCNV_med, n_noCNV_nomed), fixed)
      # Random effects
      random <- random_effect_ma(input)
      random_effects[nrow(random_effects) + 1,] <- c(c(CNV, medication, Effect, 
                                                       nCNV_med, nCNV_nomed, n_noCNV_med, n_noCNV_nomed), random)
    }
  }
}

write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))

#####################################
### 19. PRS causal mediation.     ###
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "19_PRS_causal_mediation.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/BMI/PsychiatricMedicationAnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

table_all$se <- (table_all$ci_upper - table_all$ci_lower) / (1.96*2)

fixed_effects <- data.frame(matrix(nrow=0, ncol=13))
colnames(fixed_effects) <- c('medication', 'Effect',
                             "nCNV_med", "nCNV_nomed", "n_noCNV_med", "n_noCNV_nomed",
                             'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=13))
colnames(random_effects) <- c('medication', 'Effect',
                              "nCNV_med", "nCNV_nomed", "n_noCNV_med", "n_noCNV_nomed",
                              'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')


for(j in seq_along(unique(table_all$medication))){
  medication <- unique(table_all$medication)[j]
  for(k in seq_along(unique(table_all$Effect))){
    Effect <- unique(table_all$Effect)[k]
    input <- table_all[which((table_all$CNV == CNV) & (table_all$medication == medication) & (table_all$Effect == Effect)),]
    input <- input[which(input$se != 0),]
    if(nrow(input[which(!is.na(input$se)),]) ==0){
      next
    }
    
    nCNV_med <- sum(input$nCNV_med)
    nCNV_nomed <- sum(input$nCNV_nomed)
    n_noCNV_med <- sum(input$n_noCNV_med)
    n_noCNV_nomed <- sum(input$n_noCNV_nomed)
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(medication, Effect, 
                                                   nCNV_med, nCNV_nomed, n_noCNV_med, n_noCNV_nomed), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(medication, Effect, 
                                                     nCNV_med, nCNV_nomed, n_noCNV_med, n_noCNV_nomed), random)
  }
}


write.csv(table_all, paste0('PsychiatricMedication.v2/tables/', out_string))
write.csv(fixed_effects, paste0('PsychiatricMedication.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('PsychiatricMedication.v2/random_effects/', out_string))

################################################
### 20. CNV causal mediation w/ interaction  ###
################################################

################################################
### 21. PRS causal mediation w/ interaction  ###
################################################

