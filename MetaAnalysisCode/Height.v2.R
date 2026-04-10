###############################
### Height Meta-analysis v2 ###
###.    October 2025        ###
###############################

## Setup
setwd('/path/to/your/MetaAnalysis/')
library(dplyr)
library(metafor)

cohorts_all <- c('UKBB.EUR', 'UKBB.AFR', 'UKBB.ASN', 
                 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 
                 'MyCode.EUR', 'MyCode.AFR', 'MyCode.LATNAT', 
                 'EstBB', 
                 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT', 'ClinicalAdults')
cohorts <- c('UKBB.EUR', 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 'MyCode.EUR', 'EstBB', 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT', 'ClinicalAdults')
cohorts_EUR <- c('UKBB.EUR', 'AoU.EUR', 'MyCode.EUR', 'EstBB', 'MVP.EUR', 'ClinicalAdults')

# cohorts <- c('UKBB', 'EstBB', 'Geisinger', 'MVP', 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT')
# cohorts_EUR <- c('UKBB', 'EstBB', 'Geisinger', 'MVP', 'AoU.EUR')

dir.create('Height.v2')
dir.create('Height.v2/fixed_effects')
dir.create('Height.v2/random_effects')
dir.create('Height.v2/raw_quantiles')
dir.create('Height.v2/scaled_quantiles')
dir.create('Height.v2/tables')


######################
## define functions ##
######################

fixed_effect_ma <- function(input){
  fixed <- rma.uni(yi=input[which(!is.na(input$se)),'estimate'], sei=input[which(!is.na(input$se)),'se'],method="FE")
  output <- c(sum(input[which(!is.na(input$se)),'count']), fixed$b, fixed$se, fixed$zval, fixed$pval, fixed$ci.lb, fixed$ci.ub)
  return(output)
}

random_effect_ma <- function(input) {
  tryCatch({
    random <- metafor::rma.uni(
      yi = input$estimate,
      sei = input$se,
      method = "REML"
    )
    c(sum(input$count, na.rm = TRUE),
      random$b, random$se, random$zval,
      random$pval, random$ci.lb, random$ci.ub)
  },
  error = function(e) rep(NA, 7))
}

meta_linear_quantiles <- function(dat) {
  # initialize reconstructed distribution
  all_samples <- c()
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  # for each cohort, simulate log-normal distribution based on quantiles
  for (i in 1:nrow(dat)) {
    row <- dat[i, ]
    qvals <- as.numeric(row[c("pct025_raw","pct25_raw","pct50_raw","pct75_raw","pct975_raw")])
    if(row$cohort %in% c('MyCode.EUR', 'MyCode.LATNAT', 'MyCode.AFR')){
      qvals <- qvals * 2.54
    }
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
  keep <- r2 < 1
  r2 <- r2[keep]
  n  <- n[keep]
  if(length(r2) != length(n)) stop("r2 and n must be same length")
  
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




########################
## 2. CNV main effect ##
########################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "2_CNV_main_effect.csv"

for(i in seq_along(cohorts_all)){
  cohort <- cohorts_all[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_raw) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")
  
quantiles_scaled <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_scaled) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
  # Raw quantiles
  raw <- meta_linear_quantiles(input)
  quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus), raw)
  # Scaled quantiles
  scaled <- meta_normal_quantiles(input)
  quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus), scaled)
  
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))

#############################################
## 2b. CNV main effect sex stratified meta ##
#############################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "9_CNV_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_raw) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_scaled) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
  
}

write.csv(table_all, paste0('Height.v2/tables/', '2b_CNV_main_effect.csv'))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', '2b_CNV_main_effect.csv'))
write.csv(random_effects, paste0('Height.v2/random_effects/', '2b_CNV_main_effect.csv'))


########################
## 3. PRS main effect ##
########################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "3_PRS_main_effect.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table <- table[which(table$X %in% c('(Intercept)', 'PRS')),]
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c('coeff', 'estimate', 'se', 'z', 'p', 'cohort')

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in c('(Intercept)', 'PRS')){
  input <- table_all[which(table_all$coeff == i),]
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(i, fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(i, random)
  
}


write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))

##############################
## 4. CNV x PRS interaction ##
##############################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "4_CNVxPRS_interaction.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))

############################
## 5. CNV by PRS quartile ##
############################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "5_CNV_by_PRSquartile.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('locus', 'PRSquartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('locus', 'PRSquartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'PRSquartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'PRSquartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  for(j in (1:4)){
    input <- table_all[which(table_all$locus == locus & table_all$PRSquartile == j),]
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, j), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, j), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, j), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, j), scaled)
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))


############################
## 6. CNV effect on PRS   ##
############################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "6_CNV_effect_on_PRS.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


############################
## 7. PRS by CNV          ##
############################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "7_PRS_by_CNV.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c("","locus","count","estimate","se","p","r2","r2_adj","intercept","intercept_se","intercept_p", "cohort")

fixed_effects <- data.frame(matrix(nrow=0, ncol=21))
colnames(fixed_effects) <-  c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

random_effects <- data.frame(matrix(nrow=0, ncol=21))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

for(i in seq_along(unique(table_all$locus))){
  locus <- unique(table_all$locus)[i]
  # limit to this locus
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  input <- input[which(!is.na(input$se)),]
  if(nrow(input) == 1){
    next
  }
  
  input_intercept <- input[,c("locus", "count", "intercept","intercept_se")]
  colnames(input_intercept) <- c("locus", "count", "estimate", "se")
  
  # Fixed effects
  fixed_slope <- fixed_effect_ma(input)
  fixed_r2 <- meta_r2(input$r2, input$count, 'FE')
  fixed_r2_adj <- meta_r2(input$r2_adj, input$count, 'FE')
  fixed_intercept <- fixed_effect_ma(input_intercept)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(locus, fixed_slope, fixed_r2, fixed_r2_adj, fixed_intercept)
  
  # Random effects
  random_slope <- random_effect_ma(input)
  random_r2 <- meta_r2(input$r2, input$count, 'REML')
  random_r2_adj <- meta_r2(input$r2_adj, input$count, 'REML')
  random_intercept <- random_effect_ma(input_intercept)
  random_effects[nrow(random_effects) + 1,] <- c(locus, random_slope, random_r2, random_r2_adj, random_intercept)
  
}


write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


##############################
## 8. CNV x Sex interaction ##
##############################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "8_CNVxSex_interaction.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


#######################
## 9. CNV by Sex     ##
#######################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "9_CNV_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  for(sex in c('male', 'female')){
    # limit to this locus
    locus <- unique(table_all$locus)[i]
    input <- table_all[which(table_all$locus == locus & table_all$sex == sex),]
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, sex), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, sex), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, sex), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, sex), scaled)
    
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))

#######################
## 10. PRSxSex       ##
#######################



# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "10_PRSxSex.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table <- table[which(table$X %in% c('(Intercept)', 'PRS:sex2', 'PRS')),]
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c('coeff', 'estimate', 'se', 'z', 'p', 'cohort')

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in c('(Intercept)', 'PRS:sex2', 'PRS')){
  input <- table_all[which(table_all$coeff == i),]
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(i, fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(i, random)
  
}


write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


###########################
## 11. PRS by Sex        ##
###########################



# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "11_PRS_by_Sex.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c('coeff', 'estimate', 'se', 'z', 'p', 'sex', 'cohort')

table_all$coeff <- ifelse(table_all$sex == 'male', table_all$coeff, substr(table_all$coeff, 1, nchar(table_all$coeff) - 1))
table_all <- table_all[which(table_all$coeff %in% c('(Intercept)', 'PRS')),]

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in c('(Intercept)', 'PRS')){
  for(sex in c('male', 'female')){
    input <- table_all[which(table_all$coeff == i & table_all$sex == sex),]
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(i, sex), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(i, sex), random)
    
  }
}


write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


###########################
## 12. PRS by CNV by Sex ##
###########################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "12_PRS_by_CNV_by_Sex.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=22))
colnames(fixed_effects) <-  c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

random_effects <- data.frame(matrix(nrow=0, ncol=22))
colnames(random_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  for(sex in c('male', 'female')){
    # limit to this locus
    locus <- unique(table_all$locus)[i]
    input <- table_all[which(table_all$locus == locus & table_all$sex == sex),]
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, sex), raw)
    
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, sex), scaled)
    
    input <- input[which(!is.na(input$se)),]
    if(nrow(input) <= 1){
      next
    }
    
    input_intercept <- input[,c("locus", "count", "intercept","intercept_se")]
    colnames(input_intercept) <- c("locus", "count", "estimate", "se")
    
    # Fixed effects
    fixed_slope <- fixed_effect_ma(input)
    fixed_r2 <- meta_r2(input$r2, input$count, 'FE')
    fixed_r2_adj <- meta_r2(input$r2_adj, input$count, 'FE')
    fixed_intercept <- fixed_effect_ma(input_intercept)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, sex), fixed_slope, fixed_r2, fixed_r2_adj, fixed_intercept)
    
    # Random effects
    random_slope <- random_effect_ma(input)
    random_r2 <- meta_r2(input$r2, input$count, 'REML')
    random_r2_adj <- meta_r2(input$r2_adj, input$count, 'REML')
    random_intercept <- random_effect_ma(input_intercept)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, sex), random_slope, random_r2, random_r2_adj, random_intercept)
    
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))


####################################
## 13. CNV by sex by PRS quartile ##
####################################



# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "13_CNV_by_sex_by_PRS_quartile.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('locus', 'PRSquartile', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('locus', 'PRSquartile', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('locus', 'PRSquartile', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('locus', 'PRSquartile', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  for(j in (1:4)){
    for(sex in c('male', 'female')){
      input <- table_all[which(table_all$locus == locus & table_all$quartile == j & table_all$sex == sex),]
      if(nrow(input[which(!is.na(input$estimate)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, j, sex), fixed)
      # Random effects
      random <- random_effect_ma(input)
      random_effects[nrow(random_effects) + 1,] <- c(c(locus, j, sex), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, j, sex), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, j, sex), scaled)
      
    }
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))

cohorts_all <- c('UKBB.EUR', 'UKBB.AFR', 'UKBB.ASN', 
                 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 
                 'MyCode.EUR', 'MyCode.AFR', 'MyCode.LATNAT', 
                 'EstBB', 
                 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT')
cohorts <- c('UKBB.EUR', 'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT', 'MyCode.EUR', 'EstBB', 'MVP.EUR', 'MVP.AFR', 'MVP.LATNAT')
cohorts_EUR <- c('UKBB.EUR', 'AoU.EUR', 'MyCode.EUR', 'EstBB', 'MVP.EUR')

############################
## 14. Variance explained ##
############################


anova_all <- data.frame(matrix(nrow=0, ncol=6))
colnames(anova_all) <- c("model", "Rsq", "Rsq_adj", "n", "m", "cohort")

anova_fixed <- data.frame(matrix(nrow=0, ncol=7))
colnames(anova_fixed) <- c("model", "Rsq", "ci_lower", "ci_upper", "Rsq_adj", "ci_lower_adj", "ci_upper_adj")
anova_random <- data.frame(matrix(nrow=0, ncol=7))
colnames(anova_random) <- c("model", "Rsq", "ci_lower", "ci_upper", "Rsq_adj", "ci_lower_adj", "ci_upper_adj")


dir.create('Height/fixed_effects/27_variance_explained')
dir.create('Height/random_effects/27_variance_explained')
dir.create('Height/tables/27_variance_explained')

models <- c("CNV_ANOVA.csv",
            "CNV_PRS_interaction_ANOVA.csv", 
            "PRS_ANOVA.csv",
            "covariates_ANOVA.csv",
            "CNV_PRS_additive_ANOVA.csv")
models_EUR <- c("CNV_PRS_interaction_ANOVA.csv", 
                "PRS_ANOVA.csv",
                "CNV_PRS_additive_ANOVA.csv")



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
    anova_df <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/14_variance_explained/', model))
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

write.csv(anova_all, 'Height.v2/tables/14_variance_explained.csv')
write.csv(anova_fixed, 'Height.v2/fixed_effects/14_variance_explained.csv')
write.csv(anova_random, 'Height.v2/random_effects/14_variance_explained.csv')


#####################################
## 15. multibreakpoint interaction ##
#####################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "15_multibreakpoint_int.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=37))
colnames(fixed_effects) <-  c("locusA","locusB","locusAB","genotype",
                              "locusA_count_union", "locusA_estimate","locusA_se","locusA_z", "locusA_p", "locusA_ci_lower", "locusA_ci_upper",
                              "locusB_count_union", "locusB_estimate","locusB_se","locusB_z", "locusB_p", "locusB_ci_lower", "locusB_ci_upper",
                              "locusAB_count_union", "locusAB_estimate","locusAB_se","locusAB_z", "locusAB_p", "locusAB_ci_lower", "locusAB_ci_upper",
                              "additive_R2", "additive_R2_ci_lower", "additive_R2_ci_upper",
                              "additive_R2_adj", "additive_R2_adj_ci_lower", "additive_R2_adj_ci_upper",
                              "int_R2", "int_R2_ci_lower", "inte_R2_ci_upper",
                              "int_R2_adj", "int_R2_adj_ci_lower", "int_R2_adj_ci_upper")

random_effects <- data.frame(matrix(nrow=0, ncol=37))
colnames(random_effects) <- c("locusA","locusB","locusAB","genotype",
                              "locusA_count_union", "locusA_estimate","locusA_se","locusA_z", "locusA_p", "locusA_ci_lower", "locusA_ci_upper",
                              "locusB_count_union", "locusB_estimate","locusB_se","locusB_z", "locusB_p", "locusB_ci_lower", "locusB_ci_upper",
                              "locusAB_count_union", "locusAB_estimate","locusAB_se","locusAB_z", "locusAB_p", "locusAB_ci_lower", "locusAB_ci_upper",
                              "additive_R2", "additive_R2_ci_lower", "additive_R2_ci_upper",
                              "additive_R2_adj", "additive_R2_adj_ci_lower", "additive_R2_adj_ci_upper",
                              "int_R2", "int_R2_ci_lower", "inte_R2_ci_upper",
                              "int_R2_adj", "int_R2_adj_ci_lower", "int_R2_adj_ci_upper")


########################################
### 15. Multibreakpoint int Ltest.   ###
########################################

N_mapping <- c('UKBB.EUR' = 321245,
               'UKBB.AFR' = 6438,
               'UKBB.ASN' = 13785,
               'MyCode.EUR' = 106361,
               'MyCode.AFR' = 2796,
               'MyCode.LATNAT' = 1340,
               'MVP.EUR' = 420569,
               'MVP.AFR' = 108179,
               'MVP.LATNAT' = 58996,
               'AoU.EUR' = 215995,
               'AoU.LATNAT' = 47250,
               'AoU.AAM' = 80830,
               'EstBB' = 63217)


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "15_multibreakpoint_int.dummy.Ltest.txt"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  lines <- readLines(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  for(line in lines){
    if(substr(line, 1, 6) == '[1] "l'){
      loci_clean <- gsub("^locus\\.", "", strsplit(gsub('\\[1\\]|"', '', line), " ")[[1]])
      A <- loci_clean[2]
      B <- loci_clean[3]
      AB <- loci_clean[4]
    }
    if(substr(line, 1, 6) == '[1] "D'){
      gt <- substr(line, 6, 8)
      rowname <- paste0(A,B,AB,gt, sep='_')
    }
    if(substr(line, 1, 2) == '2 '){
      p <- strsplit(line, "\\s+")[[1]][5]
      record <- c(cohort, rowname, A, B, AB, gt, p)
      table_all <- rbind(table_all, record)
    }
  }
}

colnames(table_all) <- c('cohort', 'test', 'locusA', 'locusB', 'locusAB', 'genotype', 'p')
table_all$p <- as.numeric(table_all$p)
table_all <- table_all[which(!is.na(table_all$p)),]
table_all$N <- N_mapping[table_all$cohort]
lTest_out <- data.frame(matrix(nrow=0,ncol=7))
colnames(lTest_out) <- c('test', 'locusA', 'locusB', 'locusAB', 'genotype', 'p_Fisher', 'p_Stouffer')

for(i in seq_along(unique(table_all$test))){
  test <- unique(table_all$test)[i]
  input <- table_all[which(table_all$test == test),]
  p_vals <- input$p
  weights <- input$N
  
  # Fisher's
  X2 <- -2 * sum(log(p_vals))
  p_Fisher <- 1 - pchisq(X2, 2 * length(p_vals))
  
  # Stouffer's
  z <- qnorm(1 - p_vals)
  z_weighted <- sum(weights * z) / sqrt(sum(weights^2))
  p_Stouffer <- 2 * (1 - pnorm(abs(z_weighted)))
  
  lTest_out[nrow(lTest_out) +1,] <- c(test, input$locusA[1], input$locusB[1], input$locusAB[1],
                                      input$genotype[1], p_Fisher, p_Stouffer)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(lTest_out, paste0('Height.v2/fixed_effects/', out_string))




########################################
### 16. CNV group main effect        ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "16_CNVgroup_main_effect.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_raw) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_scaled) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
  # Raw quantiles
  raw <- meta_linear_quantiles(input)
  quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus), raw)
  # Scaled quantiles
  scaled <- meta_normal_quantiles(input)
  quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus), scaled)
  
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))

########################################
### 17. CNV group x PRS              ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "17_CNVgroupxPRS_interaction.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


########################################
### 18. CNV group by PRSquartile     ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "18_CNVgroup_by_PRSquartile.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('locus', 'PRSquartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('locus', 'PRSquartile', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'PRSquartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'PRSquartile', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  if(is.na(locus)){
    next
  }
  for(j in (1:4)){
    input <- table_all[which(table_all$locus == locus & table_all$PRSquartile == j),]
    
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, j), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, j), scaled)
    
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, j), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, j), random)

  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))

########################################
### 19. CNV group effect on PRS      ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "19_CNVgroup_effect_on_PRS.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


########################################
### 20. PRS by CNV group             ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "20_PRS_by_CNVgroup.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

colnames(table_all) <- c("","locus","count","estimate","se","p","r2","r2_adj","intercept","intercept_se","intercept_p", "cohort")

fixed_effects <- data.frame(matrix(nrow=0, ncol=21))
colnames(fixed_effects) <-  c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

random_effects <- data.frame(matrix(nrow=0, ncol=21))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

for(i in seq_along(unique(table_all$locus))){
  locus <- unique(table_all$locus)[i]
  # limit to this locus
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  if(nrow(input[which(!is.na(input$se)),]) ==0){
    next
  }
  input <- input[which(!is.na(input$se)),]
  if(nrow(input) == 1){
    next
  }
  
  input_intercept <- input[,c("locus", "count", "intercept","intercept_se")]
  colnames(input_intercept) <- c("locus", "count", "estimate", "se")
  
  # Fixed effects
  fixed_slope <- fixed_effect_ma(input)
  fixed_r2 <- meta_r2(input$r2, input$count, 'FE')
  fixed_r2_adj <- meta_r2(input$r2_adj, input$count, 'FE')
  fixed_intercept <- fixed_effect_ma(input_intercept)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(locus, fixed_slope, fixed_r2, fixed_r2_adj, fixed_intercept)
  
  # Random effects
  random_slope <- random_effect_ma(input)
  random_r2 <- meta_r2(input$r2, input$count, 'REML')
  random_r2_adj <- meta_r2(input$r2_adj, input$count, 'REML')
  random_intercept <- random_effect_ma(input_intercept)
  random_effects[nrow(random_effects) + 1,] <- c(locus, random_slope, random_r2, random_r2_adj, random_intercept)
  
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))

#############################
## 21. CNVgroup x Sex      ##
#############################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "21_CNVgroupxSex_interaction.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}


fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))


#############################
## 22. CNVgroup by Sex     ##
#############################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "22_CNVgroup_by_Sex.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(fixed_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=9))
colnames(random_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  for(sex in c('male', 'female')){
    # limit to this locus
    locus <- unique(table_all$locus)[i]
    input <- table_all[which(table_all$locus == locus & table_all$sex == sex),]
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Fixed effects
    fixed <- fixed_effect_ma(input)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, sex), fixed)
    # Random effects
    random <- random_effect_ma(input)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, sex), random)
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, sex), raw)
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, sex), scaled)
    
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))


########################################
### 23. PRS by CNV group by sex      ###
########################################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "23_PRS_by_CNVgroup_by_Sex.csv"

for(i in seq_along(cohorts_EUR)){
  cohort <- cohorts_EUR[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=22))
colnames(fixed_effects) <-  c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

random_effects <- data.frame(matrix(nrow=0, ncol=22))
colnames(random_effects) <- c('locus', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper', 
                              "r2", "r2_ci_lower", "r2_ci_upper",
                              "r2_adj", "r2_adj_ci_lower", "r2_adj_ci_upper",
                              "count2", "intercept","intercept_se", "intercept_z", "intercept_p", "intercept_ci_lower", "intercept_ci_upper")

quantiles_raw <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_raw) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=8))
colnames(quantiles_scaled) <- c('locus', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  for(sex in c('male', 'female')){
    # limit to this locus
    locus <- unique(table_all$locus)[i]
    input <- table_all[which(table_all$locus == locus & table_all$sex == sex),]
    if(nrow(input[which(!is.na(input$estimate)),]) ==0){
      next
    }
    
    # Raw quantiles
    raw <- meta_linear_quantiles(input)
    quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, sex), raw)
    
    # Scaled quantiles
    scaled <- meta_normal_quantiles(input)
    quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, sex), scaled)
    
    input <- input[which(!is.na(input$se)),]
    if(nrow(input) <= 1){
      next
    }
    
    input_intercept <- input[,c("locus", "count", "intercept","intercept_se")]
    colnames(input_intercept) <- c("locus", "count", "estimate", "se")
    
    # Fixed effects
    fixed_slope <- fixed_effect_ma(input)
    fixed_r2 <- meta_r2(input$r2, input$count, 'FE')
    fixed_r2_adj <- meta_r2(input$r2_adj, input$count, 'FE')
    fixed_intercept <- fixed_effect_ma(input_intercept)
    fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, sex), fixed_slope, fixed_r2, fixed_r2_adj, fixed_intercept)
    
    # Random effects
    random_slope <- random_effect_ma(input)
    random_r2 <- meta_r2(input$r2, input$count, 'REML')
    random_r2_adj <- meta_r2(input$r2_adj, input$count, 'REML')
    random_intercept <- random_effect_ma(input_intercept)
    random_effects[nrow(random_effects) + 1,] <- c(c(locus, sex), random_slope, random_r2, random_r2_adj, random_intercept)
    
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))



############################################
### 24. CNV group by sex by PRS quartile ###
############################################


# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "24_CNVgroup_by_sex_by_PRS_quartile.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixed_effects) <- c('locus', 'PRSquartile', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=10))
colnames(random_effects) <- c('locus', 'PRSquartile', 'sex', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_raw) <- c('locus', 'PRSquartile', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=9))
colnames(quantiles_scaled) <- c('locus', 'PRSquartile', 'sex', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  for(j in (1:4)){
    for(sex in c('male', 'female')){
      input <- table_all[which(table_all$locus == locus & table_all$quartile == j & table_all$sex == sex),]
      if(nrow(input[which(!is.na(input$estimate)),]) ==0){
        next
      }
      
      # Fixed effects
      fixed <- fixed_effect_ma(input)
      fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus, j, sex), fixed)
      # Random effects
      random <- random_effect_ma(input)
      random_effects[nrow(random_effects) + 1,] <- c(c(locus, j, sex), random)
      # Raw quantiles
      raw <- meta_linear_quantiles(input)
      quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus, j, sex), raw)
      # Scaled quantiles
      scaled <- meta_normal_quantiles(input)
      quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus, j, sex), scaled)
      
    }
  }
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))




########################
## 25. 22q11.2 union  ##
########################

# Read in main effects and put into a single dataframe
table_all <- data.frame()

out_string <-  "25_22q_union_main_effects.csv"

for(i in seq_along(cohorts)){
  cohort <- cohorts[i]
  table <- read.csv(paste0('Inputs_v2/', cohort, '/Height/AnalysisOutput/', out_string))
  table$cohort <- cohort
  table_all <- rbind(table_all, table)
}

fixed_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(fixed_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

random_effects <- data.frame(matrix(nrow=0, ncol=8))
colnames(random_effects) <- c('locus', 'count', 'estimate', 'se', 'z', 'p', 'ci_lower', 'ci_upper')

quantiles_raw <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_raw) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")

quantiles_scaled <- data.frame(matrix(nrow=0, ncol=7))
colnames(quantiles_scaled) <- c('locus', 'count', "pct025", "pct25", "pct50", "pct75", "pct975")


for(i in seq_along(unique(table_all$locus))){
  # limit to this locus
  locus <- unique(table_all$locus)[i]
  input <- table_all[which(table_all$locus == locus),]
  if(nrow(input[which(!is.na(input$estimate)),]) ==0){
    next
  }
  
  # Fixed effects
  fixed <- fixed_effect_ma(input)
  fixed_effects[nrow(fixed_effects) + 1,] <- c(c(locus), fixed)
  # Random effects
  random <- random_effect_ma(input)
  random_effects[nrow(random_effects) + 1,] <- c(c(locus), random)
  # Raw quantiles
  raw <- meta_linear_quantiles(input)
  quantiles_raw[nrow(quantiles_raw) + 1,] <- c(c(locus), raw)
  # Scaled quantiles
  scaled <- meta_normal_quantiles(input)
  quantiles_scaled[nrow(quantiles_scaled) + 1,] <- c(c(locus), scaled)
  
}

write.csv(table_all, paste0('Height.v2/tables/', out_string))
write.csv(fixed_effects, paste0('Height.v2/fixed_effects/', out_string))
write.csv(random_effects, paste0('Height.v2/random_effects/', out_string))
write.csv(quantiles_raw, paste0('Height.v2/raw_quantiles/', out_string))
write.csv(quantiles_scaled, paste0('Height.v2/scaled_quantiles/', out_string))




