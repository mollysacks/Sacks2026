#######################################
#######################################
### Psychiatric Medication Analysis ###
#######################################
#######################################

# Updated December 8, 2025

setwd('~/path/to/your/working/directory/BMI')
dir.create('PsychiatricMedicationAnalysisOutput')
library(car)
library(ggplot2)
library(mediation)
library(dplyr)

################################
## 0. Setup                   ##
################################

# read in genotype matrix from old analysis code
genomatrix <- read.csv('path/to/your/genomatrix.csv')

# add in new CNV genotypes
CNV_genotypes <- read.csv('/path/to/your/genomatrix.recurrentCNV.YYYY-MM-DD.csv')
rownames(CNV_genotypes) <- as.character(CNV_genotypes$identifier)
loci <- grep("^X[0-9]+", names(CNV_genotypes), value = TRUE)
loci <- loci[!grepl("oneway$", loci)]   
loci <- setdiff(loci, c("X22q11.21_AD", "X22q11.21_DH", "X16p11.2", "X15q35", "X15q13"))
loci <- gsub("^X", "locus.", loci)
colnames(CNV_genotypes) <- gsub("^X", "locus.", colnames(CNV_genotypes))
CNV_genotypes <- CNV_genotypes[as.character(genomatrix$sampleID),c('identifier',loci)]
genomatrix <- cbind(genomatrix, CNV_genotypes)

# derive list of CNVs
CNVs <- c(paste0(loci, '_DEL'), paste0(loci, '_DUP'))


# rename PRS column
genomatrix$PRS <- genomatrix$PRS_BMI_zscore
genomatrix$PRS_quartile <- genomatrix$PRS_BMI_quantile

# make sure sex is a factor with male as ref
genomatrix$sex <- relevel(as.factor(genomatrix$sex), ref='1')


# read in main effects
main_effects <- read.csv('https://docs.google.com/spreadsheets/d/1WvwNQWiTJCkvGqCmX3jtqQjREOoP8Hu8RILZEZPoazg/export?format=csv&gid=430227892')

# add effect / frequency groupings
main_effects$group <- ifelse(main_effects$p < 0.05/nrow(main_effects),
                             ifelse(main_effects$estimate < 0, 'negative', 'positive'), 'neutral')

genomatrix$CNV_genotype_1 <- NA
genomatrix$CNV_genotype_2 <- NA
genomatrix$CNV_genotype_3 <- NA
genomatrix$CNV_genotype_4 <- NA

for(i in seq_along(main_effects$locus)){
  CNV <- main_effects$locus[i] # get name of CNV
  loc <- paste0("locus.", substr(CNV, 1, nchar(CNV) - 4)) # get locus
  cn <- ifelse(substr(CNV,nchar(CNV) - 2, nchar(CNV)) == 'DEL', 1, 3) # get copy number
  
  # get rows with CNV that don't already have another CNV
  rows1 <- which((genomatrix[[loc]] == cn) & (is.na(genomatrix$CNV_genotype_1)))  
  
  # get rows with CNV that already have 1 CNV
  rows2 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_1) & (is.na(genomatrix$CNV_genotype_2)))) 
  
  # get rows with CNV that already have 2 CNV
  rows3 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_2) & (is.na(genomatrix$CNV_genotype_3)))) 
  
  # get rows with CNV that already have 3 CNVs
  rows4 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_3) & (is.na(genomatrix$CNV_genotype_4)))) 
  
  
  # Fill in genotypes
  genomatrix[rows4, "CNV_genotype_4"] <- CNV
  genomatrix[rows3, "CNV_genotype_3"] <- CNV
  genomatrix[rows2, "CNV_genotype_2"] <- CNV
  genomatrix[rows1, "CNV_genotype_1"] <- CNV
  
}

### Create CNV group column 
# individuals with CNVs in multiple groups are NA for this
rownames(main_effects) <- main_effects$locus
genomatrix$group1 <- NA
genomatrix$group1 <- main_effects[genomatrix$CNV_genotype_1, 'group']
genomatrix$group2 <- NA
genomatrix$group2 <- main_effects[genomatrix$CNV_genotype_2, 'group']
genomatrix$group3 <- NA
genomatrix$group3 <- main_effects[genomatrix$CNV_genotype_3, 'group']
genomatrix$group4 <- NA
genomatrix$group4 <- main_effects[genomatrix$CNV_genotype_4, 'group']


genomatrix$group <- apply(genomatrix[, c("group1", "group2", "group3", "group4")], 1, function(x) {
  if (all(is.na(x))) {
    "no CNV"
  } else if ("positive" %in% x && !("negative" %in% x)) {
    "positive"
  } else if ("negative" %in% x && !("positive" %in% x)) {
    "negative"
  } else {
    "neutral"
  }
})


# set no CNV as reference
genomatrix$group <- relevel(as.factor(genomatrix$group), ref='no CNV')



# reformat group columns and add to list of CNVs

genomatrix$locus_group_negative <- ifelse(genomatrix$group == 'negative', "3",
                                          ifelse(genomatrix$group == 'no CNV', "2", NA))
genomatrix$locus_group_positive <- ifelse(genomatrix$group == 'positive', "3",
                                          ifelse(genomatrix$group == 'no CNV', "2", NA))
genomatrix$locus_group_neutral <- ifelse(genomatrix$group == 'neutral', "3",
                                         ifelse(genomatrix$group == 'no CNV', "2", NA))
CNVs <- c(CNVs, 
          c("locus_group_negative____", "locus_group_positive____", "locus_group_neutral____"))

# define list of medications and format columns
# medications <- c("all_medication", "mood_stabilizer", "antipsychotic", "antidepressant", 
#                  "ssri", "snri", "tca", "maoi", "atypical")
medications <- c("all_medication", "mood_stabilizer", "antipsychotic", "antidepressant")


numeric_meds <- sapply(genomatrix[ , medications], function(col) as.numeric(as.character(col)))

genomatrix$total_meds <- rowSums(numeric_meds, na.rm = TRUE)
genomatrix$no_meds <- ifelse(genomatrix$total_meds == 0, 1, 0)

medications_all <- c(medications, "no_meds")

for (med in medications) {
  genomatrix[[med]] <- relevel(as.factor(genomatrix[[med]]), ref = '0')
}

genomatrix$meds <- ifelse(genomatrix$all_medication == 0, 'none',
                          ifelse(genomatrix$antipsychotic == 1, 'antipsychotic',
                                 ifelse(genomatrix$mood_stabilizer == 1, 'mood_stabilizer', 'antidepressant')))
genomatrix$meds <- relevel(as.factor(genomatrix$meds), ref='none')

################################
## 1. Meds main effect        ##
################################

out <- data.frame(matrix(nrow=0, ncol=19))
colnames(out) <- c('medication', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(m in seq_along(medications)){
  # format columns
  med <- medications[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$med == "1"),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(med, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(BMI ~ 
                med +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- "med1"
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # get raw distribution
  sumstats_raw <- genomatrix %>%
    filter(med == "1") %>%
    summarise(
      mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
      p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
      p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
      p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
      p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
      p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
    )
  
  # get scaled distribution
  sumstats <- genomatrix %>%
    filter(med == "1") %>%
    summarise(
      mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
      p025  = quantile(BMI, 0.025, na.rm = TRUE),
      p25  = quantile(BMI, 0.25, na.rm = TRUE),
      p50  = quantile(BMI, 0.50, na.rm = TRUE),
      p75  = quantile(BMI, 0.75, na.rm = TRUE),
      p975  = quantile(BMI, 0.975, na.rm = TRUE)
    )
  
  # write output
  out[nrow(out) +1,] <- c(med, count, estimate, se, p, 
                          sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                          sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/1_med_main_effect.csv')


#######################
## 2. Med X CNV      ##
#######################

out <- data.frame(matrix(nrow=0, ncol=6))
colnames(out) <- c('locus', 'medication', 'count', 
                   'estimate', 'se', 'p')

out2 <- data.frame(matrix(nrow=0, ncol=7))
colnames(out2) <- c('Medication_type', 'CNV', 'count', 'proportion_CNV', 'proportion_no_CNV', 'Chi2', 'p')



for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications)){
    # format columns
    med <- medications[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    
    
    # count CNV carriers on medication
    count <- nrow(genomatrix[which(genomatrix$gen == cn & genomatrix$med == "1"),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(BMI ~
                  gen*med +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('gen', cn, ':med1')
    if(cov %in% rownames(results)){
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
    }
    else{
      estimate <- NA
      se <- NA
      p <- NA
    }
    
    
    # write output
    out[nrow(out) +1,] <- c(genotype, med, count, estimate, se, p)
    
    # Chi2 test to compare medication distribution between CNV groups
    contingency_table <- table(genomatrix$gen, genomatrix$med)[c('2', cn),]
    chi2 <- as.numeric(chisq.test(contingency_table, correct = FALSE)$statistic)
    chi2p <- as.numeric(chisq.test(contingency_table, correct = FALSE)$p.value)
    proportion_CNV <- contingency_table[2,2]/(contingency_table[2,1] + contingency_table[2,2])
    proportion_no_CNV <- contingency_table[1,2]/(contingency_table[1,1] + contingency_table[1,2])
    med_CNV <- contingency_table[2,2]
    
    out2[nrow(out2) + 1,] <- c(med, genotype, med_CNV, proportion_CNV, proportion_no_CNV, chi2, chi2p)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/2_CNVxmed.csv')
write.csv(out2, 'PsychiatricMedicationAnalysisOutput/2_CNV_med_chi2.csv')

#######################
## 3. Med by CNV     ##
#######################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'medication', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  genomatrix1 <- genomatrix[which(genomatrix$gen == cn),]
  for(m in seq_along(medications)){
    # format columns
    med <- medications[m]
    genomatrix1$med <- genomatrix1[[med]]
    genomatrix1$med <- relevel(as.factor(genomatrix1$med), ref="0")
    
    # count CNV carriers on medication
    count <- nrow(genomatrix1[which(genomatrix1$med == "1"),])
    if (count < 5){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    if ((count / nrow(genomatrix1)) == 1){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(BMI ~
                  med +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix1, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('med1')
    if(cov %in% rownames(results)){
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
    }
    else{
      estimate <- NA
      se <- NA
      p <- NA
    }
    
    # get raw distribution
    sumstats_raw <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
        p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
        p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
        p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
        p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI, 0.025, na.rm = TRUE),
        p25  = quantile(BMI, 0.25, na.rm = TRUE),
        p50  = quantile(BMI, 0.50, na.rm = TRUE),
        p75  = quantile(BMI, 0.75, na.rm = TRUE),
        p975  = quantile(BMI, 0.975, na.rm = TRUE)
      )
    
    
    # write output
    out[nrow(out) +1,] <- c(genotype, med, count, estimate, se, p,
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/3_med_by_CNV.csv')


#######################
## 4. CNV by Med     ##
#######################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'medication', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications_all)){
    # format columns
    med <- medications_all[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    genomatrix1 <- genomatrix[which(genomatrix$med == "1"),]
    
    # count CNV carriers on medication
    count <- nrow(genomatrix1[which(genomatrix1$gen == cn),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA, NA, NA, NA, 
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    if ((count / nrow(genomatrix1)) == 1){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(BMI ~
                  gen +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix1, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('gen', cn)
    if(cov %in% rownames(results)){
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
    }
    else{
      estimate <- NA
      se <- NA
      p <- NA
    }
    
    # get raw distribution
    sumstats_raw <- genomatrix1 %>%
      filter(gen == cn) %>%
      summarise(
        mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
        p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
        p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
        p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
        p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix1 %>%
      filter(gen == cn) %>%
      summarise(
        mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI, 0.025, na.rm = TRUE),
        p25  = quantile(BMI, 0.25, na.rm = TRUE),
        p50  = quantile(BMI, 0.50, na.rm = TRUE),
        p75  = quantile(BMI, 0.75, na.rm = TRUE),
        p975  = quantile(BMI, 0.975, na.rm = TRUE)
      )
    
    
    # write output
    out[nrow(out) +1,] <- c(genotype, med, count, estimate, se, p,
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/4_CNV_by_med.csv')

#######################
## 5. Med X PRS      ##
#######################

out <- data.frame(matrix(nrow=0, ncol=5))
colnames(out) <- c('medication', 'count', 
                   'estimate', 'se', 'p')

for(m in seq_along(medications)){
  # format columns
  med <- medications[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  count <- nrow(genomatrix[which(genomatrix$med == "1"),])
  
  # fit interaction model
  model <- lm(BMI ~ 
                med*prs +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- "med1:prs"
  if(cov %in% rownames(results)){
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
  }
  else{
    estimate <- NA
    se <- NA
    p <- NA
  }
  
  
  # write output
  out[nrow(out) +1,] <- c(med, count, estimate, se, p)
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/5_MedXPRS.csv')

#######################
## 6. Med by PRS     ##
#######################

out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('medication', 'quartile', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(m in seq_along(medications)){
  # format columns
  med <- medications[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  for(i in 1:4){
    genomatrix1 <- genomatrix[which(genomatrix$PRS_quartile == i),]
    count <- nrow(genomatrix1[which(genomatrix1$med == "1"),])
    
    if(count < 5){
      out[nrow(out) +1,] <- c(med, i, count, NA, NA, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit model
    model <- lm(BMI ~ 
                  med +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix1, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- "med1"
    if(cov %in% rownames(results)){
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
    }
    else{
      estimate <- NA
      se <- NA
      p <- NA
    }
    
    # get raw distribution
    sumstats_raw <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
        p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
        p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
        p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
        p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI, 0.025, na.rm = TRUE),
        p25  = quantile(BMI, 0.25, na.rm = TRUE),
        p50  = quantile(BMI, 0.50, na.rm = TRUE),
        p75  = quantile(BMI, 0.75, na.rm = TRUE),
        p975  = quantile(BMI, 0.975, na.rm = TRUE)
      )
    
    
    # write output
    out[nrow(out) +1,] <- c(med, i, count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/6_Med_by_PRS.csv')

#######################
## 7. PRS by Med     ##
#######################

out <- data.frame(matrix(nrow=0, ncol=8))
colnames(out) <- c('medication', 'count', 
                   'estimate', 'se', 'p', 'intercept_estimate', 'intercept_se', 'intercept_p')

for(m in seq_along(medications_all)){
  # format columns
  med <- medications_all[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  genomatrix1 <- genomatrix[which(genomatrix$med == "1"),]
  count <- nrow(genomatrix1)
  
  # fit interaction model
  model <- lm(BMI ~ 
                PRS +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix1, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- "PRS"
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  intercept_estimate <- as.numeric(results["(Intercept)", 1])
  intercept_se <- as.numeric(results["(Intercept)", 2])
  intercept_p <- as.numeric(results["(Intercept)", 4])
  
  # write output
  out[nrow(out) +1,] <- c(med, count, estimate, se, p, intercept_estimate, intercept_se, intercept_p)
  
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/7_PRS_by_Med.csv')

##########################
## 8. Med X Sex         ##
##########################

out <- data.frame(matrix(nrow=0, ncol=6))
colnames(out) <- c('medication', 'count_males', 'count_females',
                   'estimate', 'se', 'p')

for(m in seq_along(medications)){
  # format columns
  med <- medications[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  # count CNV carriers on medication
  count_males <- nrow(genomatrix[which(genomatrix$sex == "1" & genomatrix$med == "1"),])
  count_females <- nrow(genomatrix[which(genomatrix$sex == "2" & genomatrix$med == "1"),])
  if ((count_males < 5 )| (count_females < 5)){
    out[nrow(out) +1,] <- c(med, count_males, count_females, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(BMI ~
                sex*med +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- 'sex2:med1'
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  
  # write output
  out[nrow(out) +1,] <- c(med, count_males, count_females, estimate, se, p)
  
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/8_MedXSex.csv')

##########################
## 8. Med by Sex        ##
##########################

out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('medication', 'count', 'sex',
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(m in seq_along(medications)){
  # format columns
  med <- medications[m]
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    
    # split genomatrix
    genomatrix1 <- genomatrix[which(genomatrix$sex == i),]
    # count CNV carriers on medication
    count <- nrow(genomatrix1[which(genomatrix1$med == "1"),])
    
    if (count< 5){
      out[nrow(out) +1,] <- c(med, sex, count, NA, NA, NA, NA, NA, NA, NA, NA, 
                              NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(BMI ~
                  med +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix1, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- 'med1'
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
    
    # get raw distribution
    sumstats_raw <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
        p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
        p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
        p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
        p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix1 %>%
      filter(med == "1") %>%
      summarise(
        mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
        p025  = quantile(BMI, 0.025, na.rm = TRUE),
        p25  = quantile(BMI, 0.25, na.rm = TRUE),
        p50  = quantile(BMI, 0.50, na.rm = TRUE),
        p75  = quantile(BMI, 0.75, na.rm = TRUE),
        p975  = quantile(BMI, 0.975, na.rm = TRUE)
      )
    
    
    # write output
    out[nrow(out) +1,] <- c(med, sex, count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/9_Med_by_Sex.csv')

###########################
## 10. Med by CNV by Sex ##
###########################


out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'medication', 'sex', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')
for(i in c("1", "2")){
  sex <- ifelse(i == "1", "male", "female")
  
  for(c in seq_along(CNVs)){
    # format columns
    CNV <- CNVs[c]
    genotype <- substr(CNV, 7, nchar(CNV))
    locus <- substr(CNV, 1, nchar(CNV) - 4)
    cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
    genomatrix$gen <- genomatrix[[locus]]
    genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
    genomatrix1 <- genomatrix[which(genomatrix$gen == cn & genomatrix$sex == i),]
    for(m in seq_along(medications)){
      # format columns
      med <- medications[m]
      genomatrix1$med <- genomatrix1[[med]]
      genomatrix1$med <- relevel(as.factor(genomatrix1$med), ref="0")
      
      # count CNV carriers on medication
      count <- nrow(genomatrix1[which(genomatrix1$med == "1"),])
      if (count < 5){
        out[nrow(out) +1,] <- c(genotype, med, sex, count, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      if ((count / nrow(genomatrix1)) == 1){
        out[nrow(out) +1,] <- c(genotype, med, sex, count, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit main effects model
      model <- lm(BMI ~
                    med +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix1, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- 'med1'
      if(cov %in% rownames(results)){
        estimate <- as.numeric(results[cov, 1])
        se <- as.numeric(results[cov, 2])
        p <- as.numeric(results[cov, 4])
      }
      else{
        estimate <- NA
        se <- NA
        p <- NA
      }
      
      # get raw distribution
      sumstats_raw <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
          p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
          p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
          p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
          p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
        )
      
      # get scaled distribution
      sumstats <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI, 0.025, na.rm = TRUE),
          p25  = quantile(BMI, 0.25, na.rm = TRUE),
          p50  = quantile(BMI, 0.50, na.rm = TRUE),
          p75  = quantile(BMI, 0.75, na.rm = TRUE),
          p975  = quantile(BMI, 0.975, na.rm = TRUE)
        )
      
      
      # write output
      out[nrow(out) +1,] <- c(genotype, med, sex, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
      
    }
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/10_Med_by_CNV_by_Sex.csv')

###########################
## 11. CNV by med by Sex ##
###########################


out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'medication', 'sex', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')
for(i in c("1", "2")){
  sex <- ifelse(i == "1", "male", "female")
  
  for(c in seq_along(CNVs)){
    # format columns
    CNV <- CNVs[c]
    genotype <- substr(CNV, 7, nchar(CNV))
    locus <- substr(CNV, 1, nchar(CNV) - 4)
    cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
    genomatrix$gen <- genomatrix[[locus]]
    genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
    for(m in seq_along(medications_all)){
      # format columns
      med <- medications_all[m]
      genomatrix$med <- genomatrix[[med]]
      genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
      
      genomatrix1 <- genomatrix[which(genomatrix$med == "1" & genomatrix$sex == i),]
      
      # count CNV carriers on medication
      count <- nrow(genomatrix1[which(genomatrix1$gen == cn),])
      if (count < 5){
        out[nrow(out) +1,] <- c(genotype, med, sex, count, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      if ((count / nrow(genomatrix1)) == 1){
        out[nrow(out) +1,] <- c(genotype, med, sex, count, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit main effects model
      model <- lm(BMI ~
                    gen +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix1, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- paste0('gen', cn)
      if(cov %in% rownames(results)){
        estimate <- as.numeric(results[cov, 1])
        se <- as.numeric(results[cov, 2])
        p <- as.numeric(results[cov, 4])
      }
      else{
        estimate <- NA
        se <- NA
        p <- NA
      }
      
      # get raw distribution
      sumstats_raw <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
          p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
          p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
          p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
          p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
        )
      
      # get scaled distribution
      sumstats <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI, 0.025, na.rm = TRUE),
          p25  = quantile(BMI, 0.25, na.rm = TRUE),
          p50  = quantile(BMI, 0.50, na.rm = TRUE),
          p75  = quantile(BMI, 0.75, na.rm = TRUE),
          p975  = quantile(BMI, 0.975, na.rm = TRUE)
        )
      
      
      # write output
      out[nrow(out) +1,] <- c(genotype, med, sex, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
      
    }
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/11_CNV_by_Med_by_Sex.csv')

############################
## 12. PRS by Meds by CNV ##
############################

out <- data.frame(matrix(nrow=0, ncol=9))
colnames(out) <- c('locus', 'medication', 'count', 
                   'estimate', 'se', 'p', 'intercept_estimate', 'intercept_se', 'intercept_p')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications_all)){
    # format columns
    med <- medications_all[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    genomatrix1 <- genomatrix[which((genomatrix$med == "1") & (genomatrix$gen == cn)),]
    
    # count CNV carriers on medication
    count <- nrow(genomatrix1)
    if ((count < 10)){
      out[nrow(out) +1,] <- c(genotype, med, count, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(BMI ~
                  PRS +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix1, na.action=na.exclude)
    results <- summary(model)$coefficients
    results <- summary(model)$coefficients
    cov <- "PRS"
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
    
    intercept_estimate <- as.numeric(results["(Intercept)", 1])
    intercept_se <- as.numeric(results["(Intercept)", 2])
    intercept_p <- as.numeric(results["(Intercept)", 4])
    
    # write output
    out[nrow(out) +1,] <- c(genotype, med, count, estimate, se, p, intercept_estimate, intercept_se, intercept_p)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/12_PRS_by_CNV_by_med.csv')


############################
## 13. Meds by PRS by CNV ##
############################

out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'medication', 'quartile', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications)){
    # format columns
    med <- medications[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    for(i in 1:4){
      genomatrix1 <- genomatrix[which((genomatrix$gen == cn) & (genomatrix$PRS_quartile == i)),]
      
      # count CNV carriers on medication
      count <- nrow(genomatrix1[which(genomatrix1$med == "1"),])
      if ((count < 5 | nrow(genomatrix1) < 10)){
        out[nrow(out) +1,] <- c(genotype, med, i, count, NA, NA, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit main effects model
      model <- lm(BMI ~
                    med +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix1, na.action=na.exclude)
      results <- summary(model)$coefficients
      results <- summary(model)$coefficients
      cov <- "med1"
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
      # get raw distribution
      sumstats_raw <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
          p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
          p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
          p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
          p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
        )
      
      # get scaled distribution
      sumstats <- genomatrix1 %>%
        filter(med == "1") %>%
        summarise(
          mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI, 0.025, na.rm = TRUE),
          p25  = quantile(BMI, 0.25, na.rm = TRUE),
          p50  = quantile(BMI, 0.50, na.rm = TRUE),
          p75  = quantile(BMI, 0.75, na.rm = TRUE),
          p975  = quantile(BMI, 0.975, na.rm = TRUE)
        )
      
      # write output
      out[nrow(out) +1,] <- c(genotype, med, i, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    }
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/13_Meds_by_PRS_by_CNV.csv')


############################
## 14. CNV by PRS by Meds ##
############################


out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'medication', 'quartile', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw', 'pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications_all)){
    # format columns
    med <- medications_all[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    for(i in 1:4){
      genomatrix1 <- genomatrix[which((genomatrix$med == "1") & (genomatrix$PRS_quartile == i)),]
      
      # count CNV carriers on medication
      count <- nrow(genomatrix1[which(genomatrix1$gen == cn),])
      if ((count < 5 | nrow(genomatrix1) < 5)){
        out[nrow(out) +1,] <- c(genotype, med, i, count, NA, NA, NA, NA, NA, NA, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit main effects model
      model <- lm(BMI ~
                    gen +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix1, na.action=na.exclude)
      results <- summary(model)$coefficients
      results <- summary(model)$coefficients
      cov <- paste0("gen", cn)
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
      # get raw distribution
      sumstats_raw <- genomatrix1 %>%
        filter(gen == cn) %>%
        summarise(
          mean = mean(BMI_raw, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI_raw, 0.025, na.rm = TRUE),
          p25  = quantile(BMI_raw, 0.25, na.rm = TRUE),
          p50  = quantile(BMI_raw, 0.50, na.rm = TRUE),
          p75  = quantile(BMI_raw, 0.75, na.rm = TRUE),
          p975  = quantile(BMI_raw, 0.975, na.rm = TRUE)
        )
      
      # get scaled distribution
      sumstats <- genomatrix1 %>%
        filter(gen == cn) %>%
        summarise(
          mean = mean(BMI, na.rm = TRUE), sd = sd(BMI, na.rm = TRUE),
          p025  = quantile(BMI, 0.025, na.rm = TRUE),
          p25  = quantile(BMI, 0.25, na.rm = TRUE),
          p50  = quantile(BMI, 0.50, na.rm = TRUE),
          p75  = quantile(BMI, 0.75, na.rm = TRUE),
          p975  = quantile(BMI, 0.975, na.rm = TRUE)
        )
      
      # write output
      out[nrow(out) +1,] <- c(genotype, med, i, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    }
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/14_CNV_by_PRS_by_meds.csv')


###################################
## 15. PRS by meds by CNV by sex ##
###################################

out <- data.frame(matrix(nrow=0, ncol=10))
colnames(out) <- c('locus', 'medication', 'sex', 'count', 
                   'estimate', 'se', 'p', 'intercept_estimate', 'intercept_se', 'intercept_p')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  for(m in seq_along(medications_all)){
    # format columns
    med <- medications_all[m]
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref="0")
    for(i in c("1", "2")){
      sex <- ifelse(i == "1", "male", "female")
      
      # split genomatrix
      genomatrix1 <- genomatrix[which((genomatrix$med == "1") & (genomatrix$gen == cn) & (genomatrix$sex == i)),]
      # count CNV carriers on medication
      count <- nrow(genomatrix1)
      if ((count < 10)){
        out[nrow(out) +1,] <- c(genotype, med, sex, count, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit main effects model
      model <- lm(BMI ~
                    PRS +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix1, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- "PRS"
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
      intercept_estimate <- as.numeric(results["(Intercept)", 1])
      intercept_se <- as.numeric(results["(Intercept)", 2])
      intercept_p <- as.numeric(results["(Intercept)", 4])
      
      # write output
      out[nrow(out) +1,] <- c(genotype, med, sex, count, estimate, se, p, intercept_estimate, intercept_se, intercept_p)
    }
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/15_PRS_by_CNV_by_med_by_sex.csv')

##################################
## 16. CNV main effect + meds   ##
##################################

out <- data.frame(matrix(nrow=0, ncol=5))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$gen == cn),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(genotype, count, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(BMI ~ 
                gen + meds
              + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('gen', cn)
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(genotype, count, estimate, se, p)
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/16_CNV_main_effect_meds_factor.csv')


##################################
## 17. Variance explained.      ##
##################################

dir.create('PsychiatricMedicationAnalysisOutput/17_variance_explained')

# change all loci to factors
for(i in seq_along(loci)){
  genomatrix[[loci[i]]] <- relevel(as.factor(genomatrix[[loci[i]]]), ref="2")
}

loci_valid <- loci[sapply(genomatrix[loci], function(x) length(unique(na.omit(x))) > 1)]

# create vector of interaction terms
loci_int <- paste0(loci_valid, ':PRS')
loci_int_meds <- paste0(loci_valid, ':meds')

# initialize model
model_init <- lm(BMI ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix, na.action=na.exclude)
# medication only
model_meds <- lm(BMI ~ meds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix, na.action=na.exclude)

# CNVs only

# all loci
model_all_loci <- update(model_meds, as.formula(paste("~ . +", paste(loci_valid, collapse = " + "))))

# PRS only
model_PRS <- lm(BMI ~ meds + PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix, na.action=na.exclude)

# Additive
model_all_loci_additive <- update(model_all_loci, ~ . + PRS)

# Interaction

# all loci int
model_all_loci_int <- update(model_all_loci_additive, as.formula(paste("~ . +", paste(loci_int, collapse = " + "))))
model_all_loci_int <- update(model_all_loci_int, as.formula(paste("~ . +", paste(loci_int_meds, collapse = " + "))))
model_all_loci_int <- update(model_all_loci_int, ~ . + PRS:meds)

# Write model summaries

write.csv(data.frame(summary(model_all_loci)$coefficients), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_meds_summary.csv')
write.csv(data.frame(summary(model_all_loci_additive)$coefficients), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_PRS_meds_additive_summary.csv',)
write.csv(data.frame(summary(model_all_loci_int)$coefficients), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_PRS_meds_interaction_summary.csv')
write.csv(data.frame(summary(model_PRS)$coefficients), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/PRS_meds_summary.csv')
write.csv(data.frame(summary(model_meds)$coefficients), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/covariates_meds_summary.csv')

# Write ANOVAs
write.csv(data.frame(anova(model_all_loci)), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_meds_ANOVA.csv')
write.csv(data.frame(anova(model_all_loci_additive)), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_PRS_meds_additive_ANOVA.csv')
write.csv(data.frame(anova(model_all_loci_int)), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/CNV_PRS_meds_interaction_ANOVA.csv')
write.csv(data.frame(anova(model_PRS)), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/PRS_meds_ANOVA.csv')
write.csv(data.frame(anova(model_meds)), 'PsychiatricMedicationAnalysisOutput/17_variance_explained/covariates_meds_ANOVA.csv')






##################################
## 18. CNV Causal mediation     ##
##################################

extract_mediation_results <- function(med_summary) {
  if(med_summary$d1 == FALSE){
    row_ACME <- c('pure causal mediation effect', med_summary$d0, med_summary$d0.ci[1], med_summary$d0.ci[2], med_summary$d0.p)
    row_ADE <- c('pure direct effect', med_summary$z0, med_summary$z0.ci[1], med_summary$z0.ci[2], med_summary$z0.p)
    df <- rbind(row_ACME, row_ADE)
  }
  else{
    
    # Extract simulated draws for interaction calculation
    acme0.sims <- med_summary$d0.sims
    acme1.sims <- med_summary$d1.sims
    ade0.sims  <- med_summary$z0.sims
    ade1.sims  <- med_summary$z1.sims
    
    # Compute interaction distributions
    mediated_int.sims  <- acme1.sims - acme0.sims
    reference_int.sims <- ade1.sims  - ade0.sims
    
    # Summarize
    mediated_int_est  <- mean(mediated_int.sims)
    reference_int_est <- mean(reference_int.sims)
    
    mediated_int_ci  <- quantile(mediated_int.sims, c(0.025, 0.975))
    reference_int_ci <- quantile(reference_int.sims, c(0.025, 0.975))
    
    # Build output rows
    row_ACME <- c('pure causal mediation effect', med_summary$d0, med_summary$d0.ci[1], med_summary$d0.ci[2], med_summary$d0.p)
    row_ADE <- c('pure direct effect', med_summary$z0, med_summary$z0.ci[1], med_summary$z0.ci[2], med_summary$z0.p)
    row_med_int <- c('mediated interaction', mediated_int_est, mediated_int_ci[1], mediated_int_ci[2], NA)
    row_ref_int <- c('reference interaction', reference_int_est, reference_int_ci[1], reference_int_ci[2], NA)
    
    df <- rbind(row_ACME, row_ADE, row_med_int, row_ref_int)
  }
  colnames(df) <- c('Effect', 'estimate', 'ci_lower', 'ci_upper', 'p')
  return(data.frame(df))
}

groups <- c("locus_group_negative____", "locus_group_positive____", "locus_group_neutral____")

out <- data.frame(matrix(nrow=0, ncol=11))
colnames(out) <- c('CNV', 'medication', 'Effect', 'estimate', 'ci_lower', 'ci_upper', 'p','nCNV_med', 'nCNV_nomed', 'n_noCNV_med', 'n_noCNV_nomed') 


for(c in seq_along(groups)){
  # format columns
  CNV <- groups[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(m in seq_along(medications)){
    med <- medications[m]
    if(med == 'no_meds'){
      next
    }
    print(med)
    print(genotype)
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref='0')
    genomatrix_m <- genomatrix[which(genomatrix$gen %in% c(cn, '2')),]
    genonatrix_m <- genomatrix_m[which(genomatrix_m$med == 1 | genomatrix_m$all_medication == 0),]
    
    # count CNV carriers
    counts <- table(genomatrix_m$gen, genomatrix_m$med)
    if(!(cn %in% genomatrix$gen)){
      next
    }
    if ((counts[cn, '1'] < 10)){
      out[nrow(out) +1,] <- c(genotype, med, NA, NA, NA, NA, NA, 
                              counts[cn, '1'], counts[cn, '0'], counts['2', '1'], counts['2', '0'])
      next
    }
    
    
    # Mediator model: med ~ Genotype
    model_m <- glm(med ~ gen + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix_m, family='binomial')
    
    # Outcome model: BMI ~ Genotype + med
    model_y <- lm(BMI ~ gen + med + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix_m)
    
    # Mediation analysis
    med.out <- mediate(
      model.m = model_m,
      model.y = model_y,
      treat = "gen",
      mediator = "med",
      sims = 50,
      robustSE = TRUE,
      keep.data=FALSE
    )
    res1 <- extract_mediation_results(summary(med.out))
    res1$CNV <- genotype
    res1$medication <- med
    res1$nCNV_med <- counts[cn, '1']
    res1$nCNV_nomed <- counts[cn, '0']
    res1$n_noCNV_med <- counts['2', '1']
    res1$n_noCNV_nomed <- counts['2', '0']
    res1 <- res1[,colnames(out)]
    out <- rbind(out, res1)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/18_CNV_causal_mediation.csv')


##################################
## 19. PRS Causal mediation     ##
##################################

out <- data.frame(matrix(nrow=0, ncol=8))
colnames(out) <- c('medication', 'Effect', 'estimate', 'ci_lower', 'ci_upper', 'p','n_med', 'n_nomed')

for(m in seq_along(medications)){
  med <- medications[m]
  if(med == 'no_meds'){
    next
  }
  print(med)
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref='0')
  genomatrix_m <- genomatrix[which(genomatrix$med == 1 | genomatrix$all_medication == 0),]
  
  # count CNV carriers
  counts <- table(genomatrix_m$med)
  
  if ((counts['1'] < 10)){
    out[nrow(out) +1,] <- c(genotype, med, NA, NA, NA, NA, NA, 
                            counts['1'], counts['0'])
    next
  }
  
  
  # Mediator model: med ~ Genotype
  model_m <- glm(med ~ prs + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix_m, family='binomial')
  
  # Outcome model: BMI ~ Genotype + med
  model_y <- lm(BMI ~ prs + med + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix_m)
  
  # Mediation analysis
  med.out <- mediate(
    model.m = model_m,
    model.y = model_y,
    treat = "prs",
    mediator = "med",
    sims = 50,
    robustSE = TRUE,
    keep.data=FALSE
  )
  res1 <- extract_mediation_results(summary(med.out))
  res1$medication <- med
  res1$n_med <- counts['1']
  res1$n_nomed <- counts['0']
  res1 <- res1[,colnames(out)]
  out <- rbind(out, res1)
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/19_PRS_causal_mediation.csv')


#############################################
## 20. CNV causal mediation w/ interaction ##
#############################################

groups <- c("locus_group_negative____", "locus_group_positive____", "locus_group_neutral____")

out <- data.frame(matrix(nrow=0, ncol=11))
colnames(out) <- c('CNV', 'medication', 'Effect', 'estimate', 'ci_lower', 'ci_upper', 'p','nCNV_med', 'nCNV_nomed', 'n_noCNV_med', 'n_noCNV_nomed') 


for(c in seq_along(groups)){
  # format columns
  CNV <- groups[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(m in seq_along(medications)){
    med <- medications[m]
    if(med == 'no_meds'){
      next
    }
    print(med)
    print(genotype)
    genomatrix$med <- genomatrix[[med]]
    genomatrix$med <- relevel(as.factor(genomatrix$med), ref='0')
    genomatrix_m <- genomatrix[which(genomatrix$gen %in% c(cn, '2')),]
    genonatrix_m <- genomatrix_m[which(genomatrix_m$med == 1 | genomatrix_m$all_medication == 0),]
    
    # count CNV carriers
    counts <- table(genomatrix_m$gen, genomatrix_m$med)
    if(!(cn %in% genomatrix$gen)){
      next
    }
    if ((counts[cn, '1'] < 10)){
      out[nrow(out) +1,] <- c(genotype, med, NA, NA, NA, NA, NA, 
                              counts[cn, '1'], counts[cn, '0'], counts['2', '1'], counts['2', '0'])
      next
    }
    
    
    # Mediator model: med ~ Genotype
    model_m <- glm(med ~ gen + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix_m, family='binomial')
    
    # Outcome model: BMI ~ Genotype + med
    model_y <- lm(BMI ~ gen*med + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                  data = genomatrix_m)
    
    # Mediation analysis
    med.out <- mediate(
      model.m = model_m,
      model.y = model_y,
      treat = "gen",
      mediator = "med",
      sims = 50,
      robustSE = TRUE,
      keep.data=FALSE
    )
    res1 <- extract_mediation_results(summary(med.out))
    res1$CNV <- genotype
    res1$medication <- med
    res1$nCNV_med <- counts[cn, '1']
    res1$nCNV_nomed <- counts[cn, '0']
    res1$n_noCNV_med <- counts['2', '1']
    res1$n_noCNV_nomed <- counts['2', '0']
    res1 <- res1[,colnames(out)]
    out <- rbind(out, res1)
  }
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/20_CNV_causal_mediation_interaction.csv')


##########################################
## 21. PRS Causal mediation interaction ##
##########################################

out <- data.frame(matrix(nrow=0, ncol=8))
colnames(out) <- c('medication', 'Effect', 'estimate', 'ci_lower', 'ci_upper', 'p','n_med', 'n_nomed')

for(m in seq_along(medications)){
  med <- medications[m]
  if(med == 'no_meds'){
    next
  }
  print(med)
  genomatrix$med <- genomatrix[[med]]
  genomatrix$med <- relevel(as.factor(genomatrix$med), ref='0')
  genomatrix_m <- genomatrix[which(genomatrix$med == 1 | genomatrix$all_medication == 0),]
  
  # count CNV carriers
  counts <- table(genomatrix_m$med)
  
  if ((counts['1'] < 10)){
    out[nrow(out) +1,] <- c(genotype, med, NA, NA, NA, NA, NA, 
                            counts['1'], counts['0'])
    next
  }
  
  
  # Mediator model: med ~ Genotype
  model_m <- glm(med ~ prs + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix_m, family='binomial')
  
  # Outcome model: BMI ~ Genotype + med
  model_y <- lm(BMI ~ prs*med + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                data = genomatrix_m)
  
  # Mediation analysis
  med.out <- mediate(
    model.m = model_m,
    model.y = model_y,
    treat = "prs",
    mediator = "med",
    sims = 50,
    robustSE = TRUE,
    keep.data=FALSE
  )
  res1 <- extract_mediation_results(summary(med.out))
  res1$medication <- med
  res1$n_med <- counts['1']
  res1$n_nomed <- counts['0']
  res1 <- res1[,colnames(out)]
  out <- rbind(out, res1)
}

write.csv(out, 'PsychiatricMedicationAnalysisOutput/21_PRS_causal_mediation_interaction.csv')
