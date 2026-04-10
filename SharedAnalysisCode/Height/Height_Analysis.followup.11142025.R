##################################
##################################
### Height Analysis            ###
##################################
##################################

# Updated December 8, 2025

setwd('~/path/to/your/working/directory/Height')
dir.create('AnalysisOutput')
library(car)
library(ggplot2)
library(lmtest)
library(sandwich)
library(dplyr)

################################
## 0. Setup                   ##
################################

# read in genotype matrix from old analysis code
genomatrix <- read.csv('path/to/your/genomatrix.csv')

# add in new CNV genotypes
CNV_genotypes <- read.csv('path/to/your/genomatrix.recurrentCNV.YYYY-MM-DD.csv')
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
genomatrix$PRS <- genomatrix$PRS_HEIGHT_zscore
genomatrix$PRS_quartile <- genomatrix$PRS_HEIGHT_quantile

# make sure sex is a factor with male as ref
genomatrix$sex <- relevel(as.factor(genomatrix$sex), ref='1')


##################################
## 1. Age and Height distributions ##
##################################

# plot Height, age, and sex distributions
age_hist <- ggplot(genomatrix, aes(x=age)) +
  geom_histogram() +
  theme_classic() +
  theme(text = element_text(size=24)) +
  xlim(c(18,108))

ggsave('AnalysisOutput/age_distribution.png', age_hist, height = 8, width= 8)

Height_hist <- ggplot(genomatrix, aes(x=HEIGHT_raw, color=sex, fill=sex)) +
  geom_histogram(position = "dodge",binwidth = 1) +
  theme_classic() +
  scale_fill_manual(values = c("1" = "skyblue3", "2" = "salmon"), 
                    labels = c("1" = "Male", "2" = "Female")) +
  scale_color_manual(values = c("1" = "skyblue3", "2" = "salmon"), 
                     labels = c("1" = "Male", "2" = "Female")) +
  theme(text = element_text(size=24)) +
  xlim(c(130, 230))

ggsave('AnalysisOutput/raw_Height_distribution.png', Height_hist, height = 8, width= 8)


##################
### R2.         ##
##################

# re-write summary to include R2

model <- lm(HEIGHT ~ PRS +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

coeff <- summary(model)$coeff

# write text output to include R2
sink("AnalysisOutput/3_PRS_main_effect.txt")
summary(model)
sink()

model <- lm(HEIGHT ~ PRS*sex +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

# write text output to include R2
sink("AnalysisOutput/10_PRSxSex.txt")
summary(model)
sink()

model_male <- lm(HEIGHT ~ PRS +
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix[which(genomatrix$sex == "1"),], na.action=na.exclude)
model_female <- lm(HEIGHT ~ PRS +
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix[which(genomatrix$sex == "2"),], na.action=na.exclude)


# write text output to include R2
sink("AnalysisOutput/11_PRS_by_Sex.txt")
summary(model_male)
summary(model_female)
sink()




#######################################
## 15b. multi breakpoint interaction ##
#######################################




# redo old model, but record variance/covariance matrix

out <- data.frame(matrix(nrow=0, ncol=22))
colnames(out) <- c('locusA', 'locusB', 'locusAB', 'genotype',
                   'locusA_count_union', 'locusB_count_union', 'locusA_locusB_count',
                   'locusA_estimate', 'locusA_se', 'locusA_p', 'locusA_additive_VIF',
                   'locusB_estimate', 'locusB_se', 'locus2_p', 'locusB_additive_VIF',
                   'interaction_estimate', 'interaction_se', 'interaction_p',
                   'additive_R2', 'additive_R2_adj', 'int_R2', 'int_R2_adj')

loci_trios <- list(c('locus.22q11.21_A_B', 'locus.22q11.21_B_D', 'locus.22q11.21_A_D'),
                   c('locus.22q11.21_A_C', 'locus.22q11.21_C_D', 'locus.22q11.21_A_D'),
                   c('locus.22q11.21_B_C', 'locus.22q11.21_C_D', 'locus.22q11.21_B_D'))


for(t in seq_along(loci_trios)){
  trio <- loci_trios[[t]]
  genomatrix$locusA_specific <- relevel(as.factor(genomatrix[[trio[1]]]), ref="2")
  genomatrix$locusB_specific <- relevel(as.factor(genomatrix[[trio[2]]]), ref="2")
  genomatrix$locusAB_specific <- relevel(as.factor(genomatrix[[trio[3]]]), ref="2")
  # D_A, D_B, D_AB are 0/1; covariates in data.frame covars
  fit <- lm(HEIGHT ~ locusA_specific  + locusB_specific + locusAB_specific +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 , data=genomatrix)
  # Robust SEs
  coeftest(fit, vcov = vcovHC(fit, type="HC3"))
  write.csv(summary(fit)$coefficients, 'AnalysisOutput/15_multibreakpoint_int.dummy.modelsummary.csv')
  
  # Contrast: L = alpha_AB - alpha_A - alpha_B
  if(('locusA_specific1' %in% rownames(summary(fit)$coefficients)) && ('locusB_specific1' %in% rownames(summary(fit)$coefficients)) && ('locusAB_specific1' %in% rownames(summary(fit)$coefficients))){
    coef_names <- names(coef(fit))
    # start with a zero vector whose names match the model coefficients
    zeros <- setNames(rep(0, length(coef_names)), coef_names)
    del_target <- c("locusA_specific1" = -1,
                    "locusB_specific1" = -1,
                    "locusAB_specific1" =  1)
    # only keep those that actually exist in this model
    del_present <- intersect(names(del_target), coef_names)
    L_DEL_vec <- zeros
    if(length(del_present) > 0) L_DEL_vec[del_present] <- del_target[del_present]
    L_DEL <- rbind(epistasis = L_DEL_vec)
    sink("AnalysisOutput/15_multibreakpoint_int.dummy.Ltest.txt", append = TRUE)
    testDEL <- linearHypothesis(fit, L_DEL, vcov = vcovHC(fit, type="HC3")) 
    print(trio)
    print("DEL")
    print(testDEL)
    sink()
  }
  
  if(('locusA_specific3' %in% rownames(summary(fit)$coefficients)) && ('locusB_specific3' %in% rownames(summary(fit)$coefficients)) && ('locusAB_specific3' %in% rownames(summary(fit)$coefficients))){
    coef_names <- names(coef(fit))
    # start with a zero vector whose names match the model coefficients
    zeros <- setNames(rep(0, length(coef_names)), coef_names)
    dup_target <- c("locusA_specific3" = -1,
                    "locusB_specific3" = -1,
                    "locusAB_specific3" =  1)
    # only keep those that actually exist in this model
    dup_present <- intersect(names(dup_target), coef_names)
    L_DUP_vec <- zeros
    if(length(dup_present) > 0) L_DUP_vec[dup_present] <- dup_target[dup_present]
    L_DUP <- rbind(epistasis = L_DUP_vec)
    sink("AnalysisOutput/15_multibreakpoint_int.dummy.Ltest.txt", append = TRUE)
    testDUP <- linearHypothesis(fit, L_DUP, vcov = vcovHC(fit, type="HC3")) 
    print(trio)
    print("DUP")
    print(testDUP)
    sink()
  }
}
###########################
### Rerun everything    ###
### with CNV group.     ###
###########################


# read in main effects
main_effects <- read.csv('https://docs.google.com/spreadsheets/d/1WvwNQWiTJCkvGqCmX3jtqQjREOoP8Hu8RILZEZPoazg/export?format=csv&gid=1494408220')

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



groups_all <- c('positive', 'no CNV', 'negative', 'neutral')
groups <- c('positive', 'negative', 'neutral')




###############################
## 16. CNV group main effect ##
###############################

out <- data.frame(matrix(nrow=0, ncol=19))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(groups_all)){
  # format columns
  group <- groups_all[c]
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$group == group),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    next
  }
  
  if(group == 'no CNV'){
    estimate <- NA
    se <- NA
    p <- NA
  }
  else{
    # fit main effects model
    model <- lm(HEIGHT ~ 
                  group +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('group', group)
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
  }
  
  # get raw distribution
  sumstats_raw <- genomatrix %>%
    filter(group == groups_all[c]) %>%
    summarise(
      mean = mean(HEIGHT_raw, na.rm = TRUE),
      sd = sd(HEIGHT_raw, na.rm = TRUE),
      p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
      p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
      p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
      p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
      p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
    )
  
  # get scaled distribution
  sumstats <- genomatrix %>%
    filter(group == groups_all[c]) %>%
    summarise(
      mean = mean(HEIGHT, na.rm = TRUE),
      sd = sd(HEIGHT, na.rm = TRUE),
      p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
      p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
      p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
      p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
      p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
    )
  
  # write output
  out[nrow(out) +1,] <- c(group, count, estimate, se, p, 
                          sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                          sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
}

write.csv(out, 'AnalysisOutput/16_CNVgroup_main_effect.csv')

###################################
## 17. CNVgroupxPRS              ##
###################################

out <- data.frame(matrix(nrow=0, ncol=5))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p')

for(c in seq_along(groups)){
  # format columns
  group <- groups[c]
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$group == group),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(HEIGHT ~ 
                PRS*group +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('PRS:group', group)
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(group, count, estimate, se, p)
}

write.csv(out, 'AnalysisOutput/17_CNVgroupxPRS_interaction.csv')

###################################
## 18. CNV group by PRS quartile ##
###################################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'PRSquartile', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(groups_all)){
  # format columns
  group <- groups_all[c]
  
  for(i in (1:4)){
    # subset to PRS quartile i
    genomatrix_q <- genomatrix[which(genomatrix$PRS_quartile == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_q[which(genomatrix_q$group == group),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(group, i, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    if(group == 'no CNV'){
      estimate <- NA
      se <- NA
      p <- NA
    }
    else{
      # fit main effects model
      model <- lm(HEIGHT ~ 
                    group +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = genomatrix_q, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- paste0('group', group)
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
    }
    
    
    # get raw distribution
    sumstats_raw <- genomatrix_q %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT_raw, na.rm = TRUE),
        sd = sd(HEIGHT_raw, na.rm = TRUE),
        p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix_q %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT, na.rm = TRUE),
        sd = sd(HEIGHT, na.rm = TRUE),
        p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
      )
    
    # write output
    out[nrow(out) +1,] <- c(group, i,count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
  
}

write.csv(out, 'AnalysisOutput/18_CNVgroup_by_PRSquartile.csv')

#################################
## 19. CNV group effect on PRS ##
#################################

out <- data.frame(matrix(nrow=0, ncol=5))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p')

for(c in seq_along(groups)){
  # format columns
  group <- groups[c]
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$group == group),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(PRS ~ 
                group +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('group', group)
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(group, count, estimate, se, p)
}

write.csv(out, 'AnalysisOutput/19_CNVgroup_effect_on_PRS.csv')

############################
### 20. PRS by CNV group ###
############################

out <- data.frame(matrix(nrow=0, ncol=10))
colnames(out) <- c('locus', 'count', 'estimate', 'se', 'p', 'r2', 'r2_adj', 'intercept', 'intercept_se', 'intercept_p')

for(c in seq_along(groups_all)){
  # format columns
  group <- groups_all[c]
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$group == group),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA, NA, NA, NA, NA, NA)
    next
  }
  
  genomatrix_c <- genomatrix[which(genomatrix$group == group),]
  
  # fit PRS model
  model <- lm(HEIGHT ~ 
                PRS +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix_c, na.action=na.exclude)
  results <- summary(model)$coefficients
  estimate <- as.numeric(results['PRS', 1])
  se <- as.numeric(results['PRS', 2])
  p <- as.numeric(results['PRS', 4])
  
  intercept <- as.numeric(results['(Intercept)', 1])
  intercept_se <- as.numeric(results['(Intercept)', 2])
  intercept_p <- as.numeric(results['(Intercept)', 4])
  
  r2 <- summary(model)$r.squared
  r2_adj <- summary(model)$adj.r.squared
  
  # write output
  out[nrow(out) +1,] <- c(group, count, estimate, se, p, r2, r2_adj, intercept, intercept_se, intercept_p)
}

write.csv(out, 'AnalysisOutput/20_PRS_by_CNVgroup.csv')

############################
## 21. CNV groupxSex      ##
############################


out <- data.frame(matrix(nrow=0, ncol=5))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p')

for(c in seq_along(groups)){
  # format columns
  group <- groups[c]
  
  # count CNV carriers
  count <- nrow(genomatrix[which(genomatrix$group == group),])
  if ((count < 5)){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(HEIGHT ~ 
                sex*group +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('sex2:group', group)
  if(!(cov %in% rownames(results))){
    out[nrow(out) +1,] <- c(group, count, NA, NA, NA)
    next
  }
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(group, count, estimate, se, p)
}

write.csv(out, 'AnalysisOutput/21_CNVgroupxSex_interaction.csv')


############################
## 22. CNVgroup by sex    ##
############################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'sex', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')


for(c in seq_along(groups_all)){
  # format columns
  group <- groups_all[c]
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    
    # split genomatrix
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_s[which(genomatrix_s$group == group),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(group, sex, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    if(group == 'no CNV'){
      estimate <- NA
      se <- NA
      p <- NA
    }
    else{
      # fit main effects model
      model <- lm(HEIGHT ~ 
                    group +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = genomatrix_s, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- paste0('group', group)
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
    }
    
    # get raw distribution
    sumstats_raw <- genomatrix_s %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT_raw, na.rm = TRUE),
        sd = sd(HEIGHT_raw, na.rm = TRUE),
        p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix_s %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT, na.rm = TRUE),
        sd = sd(HEIGHT, na.rm = TRUE),
        p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
      )
    
    # write output
    out[nrow(out) +1,] <- c(group, sex ,count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'AnalysisOutput/22_CNVgroup_by_Sex.csv')


##################################
## 23. PRS by CNV group by sex  ##
##################################

out <- data.frame(matrix(nrow=0, ncol=25))
colnames(out) <- c('locus', 'sex', 'count', 'estimate', 'se', 'p', 'r2', 'r2_adj', 'intercept', 'intercept_se', 'intercept_p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(groups_all)){
  # format columns
  group <- groups_all[c]
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_s[which(genomatrix_s$group == group),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(group, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    genomatrix_c <- genomatrix_s[which(genomatrix_s$group == group),]
    
    # fit PRS model
    model <- lm(HEIGHT ~ 
                  PRS +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix_c, na.action=na.exclude)
    results <- summary(model)$coefficients
    estimate <- as.numeric(results['PRS', 1])
    se <- as.numeric(results['PRS', 2])
    p <- as.numeric(results['PRS', 4])
    
    intercept <- as.numeric(results['(Intercept)', 1])
    intercept_se <- as.numeric(results['(Intercept)', 2])
    intercept_p <- as.numeric(results['(Intercept)', 4])
    
    r2 <- summary(model)$r.squared
    r2_adj <- summary(model)$adj.r.squared
    
    
    # get raw distribution
    sumstats_raw <- genomatrix_c %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT_raw, na.rm = TRUE),
        sd = sd(HEIGHT_raw, na.rm = TRUE),
        p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
      )
    
    # get scaled distribution
    sumstats <- genomatrix_c %>%
      filter(group == groups_all[c]) %>%
      summarise(
        mean = mean(HEIGHT, na.rm = TRUE),
        sd = sd(HEIGHT, na.rm = TRUE),
        p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
        p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
        p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
        p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
        p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
      )
    
    # write output
    out[nrow(out) +1,] <- c(group, sex, count, estimate, se, p, r2, r2_adj, intercept, intercept_se, intercept_p,
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'AnalysisOutput/23_PRS_by_CNVgroup_by_sex.csv')


#########################################
## 24. CNV group by sex by PRSquartile ##
#########################################

out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'sex', 'quartile', 'count', 'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(groups)){
  # format columns
  group <- groups[c]
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    for(j in 1:4){
      genomatrix_q <- genomatrix_s[which(genomatrix_s$PRS_quartile == j),]
      # count CNV carriers
      count <- nrow(genomatrix_q[which(genomatrix_q$group == group),])
      if ((count < 5)){
        out[nrow(out) +1,] <- c(group, sex, j, count, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit PRS model
      model <- lm(HEIGHT ~ 
                    group +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = genomatrix_q, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- paste0('group', group)
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
      
      # get raw distribution
      sumstats_raw <- genomatrix_q %>%
        filter(group == groups[c]) %>%
        summarise(
          mean = mean(HEIGHT_raw, na.rm = TRUE),
          sd = sd(HEIGHT_raw, na.rm = TRUE),
          p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
          p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
          p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
          p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
          p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
        )
      
      # get scaled distribution
      sumstats <- genomatrix_q %>%
        filter(group == groups[c]) %>%
        summarise(
          mean = mean(HEIGHT, na.rm = TRUE),
          sd = sd(HEIGHT, na.rm = TRUE),
          p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
          p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
          p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
          p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
          p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
        )
      
      # write output
      out[nrow(out) +1,] <- c(group, sex, j, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
      
    }
  }
}

write.csv(out, 'AnalysisOutput/24_CNVgroup_by_sex_by_PRS_quartile.csv')




AB_loci <- c("locus.22q11.21_A_B", "locus.22q11.21_A_C",  "locus.22q11.21_A_D")
BC_loci <- c("locus.22q11.21_B_C", "locus.22q11.21_A_C",  "locus.22q11.21_B_D", "locus.22q11.21_A_D")
CD_loci <- c("locus.22q11.21_C_D", "locus.22q11.21_B_D",  "locus.22q11.21_A_D")
AC_loci <- c("locus.22q11.21_A_C",  "locus.22q11.21_A_D")
BD_loci <- c("locus.22q11.21_B_D",  "locus.22q11.21_A_D")


genomatrix$locus.22q11.21_A_B_union <- with(genomatrix,
                                            ifelse(rowSums(genomatrix[AB_loci] == 1) > 0, 1,
                                                   ifelse(rowSums(genomatrix[AB_loci] == 3) > 0, 3, 2)))

genomatrix$locus.22q11.21_B_C_union <- with(genomatrix,
                                            ifelse(rowSums(genomatrix[BC_loci] == 1) > 0, 1,
                                                   ifelse(rowSums(genomatrix[BC_loci] == 3) > 0, 3, 2)))

genomatrix$locus.22q11.21_C_D_union <- with(genomatrix,
                                            ifelse(rowSums(genomatrix[CD_loci] == 1) > 0, 1,
                                                   ifelse(rowSums(genomatrix[CD_loci] == 3) > 0, 3, 2)))

genomatrix$locus.22q11.21_A_C_union <- with(genomatrix,
                                            ifelse(rowSums(genomatrix[AC_loci] == 1) > 0, 1,
                                                   ifelse(rowSums(genomatrix[AC_loci] == 3) > 0, 3, 2)))

genomatrix$locus.22q11.21_B_D_union <- with(genomatrix,
                                            ifelse(rowSums(genomatrix[BD_loci] == 1) > 0, 1,
                                                   ifelse(rowSums(genomatrix[BD_loci] == 3) > 0, 3, 2)))




out <- data.frame(matrix(nrow=0, ncol=19))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

CNVs <- c('locus.22q11.21_A_B_union_DEL', 
          'locus.22q11.21_B_C_union_DEL', 
          'locus.22q11.21_C_D_union_DEL', 
          'locus.22q11.21_A_C_union_DEL', 
          'locus.22q11.21_B_D_union_DEL',
          'locus.22q11.21_A_B_union_DUP', 
          'locus.22q11.21_B_C_union_DUP', 
          'locus.22q11.21_C_D_union_DUP', 
          'locus.22q11.21_A_C_union_DUP', 
          'locus.22q11.21_B_D_union_DUP')

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
    out[nrow(out) +1,] <- c(genotype, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    next
  }
  
  # fit main effects model
  model <- lm(HEIGHT ~ 
                gen +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('gen', cn)
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # get raw distribution
  sumstats_raw <- genomatrix %>%
    filter(gen == cn) %>%
    summarise(
      mean = mean(HEIGHT_raw, na.rm = TRUE),
      sd = sd(HEIGHT_raw, na.rm = TRUE),
      p025  = quantile(HEIGHT_raw, 0.025, na.rm = TRUE),
      p25  = quantile(HEIGHT_raw, 0.25, na.rm = TRUE),
      p50  = quantile(HEIGHT_raw, 0.50, na.rm = TRUE),
      p75  = quantile(HEIGHT_raw, 0.75, na.rm = TRUE),
      p975  = quantile(HEIGHT_raw, 0.975, na.rm = TRUE)
    )
  
  # get scaled distribution
  sumstats <- genomatrix %>%
    filter(gen == cn) %>%
    summarise(
      mean = mean(HEIGHT, na.rm = TRUE),
      sd = sd(HEIGHT, na.rm = TRUE),
      p025  = quantile(HEIGHT, 0.025, na.rm = TRUE),
      p25  = quantile(HEIGHT, 0.25, na.rm = TRUE),
      p50  = quantile(HEIGHT, 0.50, na.rm = TRUE),
      p75  = quantile(HEIGHT, 0.75, na.rm = TRUE),
      p975  = quantile(HEIGHT, 0.975, na.rm = TRUE)
    )
  
  # write output
  out[nrow(out) +1,] <- c(genotype, count, estimate, se, p, 
                          sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                          sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
}


write.csv(out, 'AnalysisOutput/25_22q_union_main_effects.csv')
