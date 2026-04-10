##################################
##################################
### Height Analysis               ###
##################################
##################################

# Updated July 28, 2025

setwd('~/path/to/your/working/directory/Height')
dir.create('AnalysisOutput')
library(car)
library(ggplot2)

################################
## 0. Setup                   ##
################################

# read in genotype matrix
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
# 
# # reformat group columns and add to list of CNVs
# genomatrix$locus_group_negative <- ifelse(genomatrix$group_c == 'negative', "3",
#                                     ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# genomatrix$locus_group_positive <- ifelse(genomatrix$group_c == 'positive', "3",
#                                     ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# genomatrix$locus_group_negative_rare <- ifelse(genomatrix$group == 'negative_rare', "3",
#                                     ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# genomatrix$locus_group_positive_rare <- ifelse(genomatrix$group == 'positive_rare', "3",
#                                          ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# genomatrix$locus_group_negative_common <- ifelse(genomatrix$group == 'negative', "3",
#                                          ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# genomatrix$locus_group_positive_common <- ifelse(genomatrix$group == 'positive', "3",
#                                          ifelse(genomatrix$group_c == 'no CNV', "2", NA))
# CNVs <- c(CNVs_no_groups, 
#           c("locus_group_negative____", "locus_group_positive____", 
#             "locus_group_negative_rare____", "locus_group_positive_rare____", 
#             "locus_group_negative_common____", "locus_group_positive_common____"))
##################################
## 1. Age and Height distributions ##
##################################

# plot Height, age, and sex distributions
age_hist <- ggplot(genomatrix, aes(x=age)) +
  geom_histogram() +
  theme_classic() +
  theme(text = element_text(size=24)) +
  xlim(c(18,85))

ggsave('AnalysisOutput/age_distribution.png', age_hist, height = 8, width= 8)

Height_hist <- ggplot(genomatrix, aes(x=HEIGHT_raw, color=sex, fill=sex)) +
  geom_histogram(position = "dodge",binwidth = 1) +
  theme_classic() +
  scale_fill_manual(values = c("1" = "skyblue3", "2" = "salmon"), 
                    labels = c("1" = "Male", "2" = "Female")) +
  scale_color_manual(values = c("1" = "skyblue3", "2" = "salmon"), 
                     labels = c("1" = "Male", "2" = "Female")) +
  theme(text = element_text(size=24)) +
  xlim(c(15, 65))

ggsave('AnalysisOutput/raw_Height_distribution.png', Height_hist, height = 8, width= 8)

#############################
## 2. CNV main effects     ##
#############################

out <- data.frame(matrix(nrow=0, ncol=19))
colnames(out) <- c('locus', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

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

write.csv(out, 'AnalysisOutput/2_CNV_main_effect.csv')


#############################
## 3. PRS main effect.     ##
#############################

model <- lm(HEIGHT ~ PRS +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

coeff <- summary(model)$coeff
write.csv(coeff, 'AnalysisOutput/3_PRS_main_effect.csv')

# write text output to include R2
sink("AnalysisOutput/3_PRS_main_effect.txt")
summary(summary(model))
sink()

#############################
## 4. CNVxPRS              ##
#############################

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
  model <- lm(HEIGHT ~ 
                PRS*gen +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('PRS:gen', cn)
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(genotype, count, estimate, se, p)
}

write.csv(out, 'AnalysisOutput/4_CNVxPRS_interaction.csv')

############################
## 5. CNV by PRS quartile ##
############################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'PRSquartile', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(i in (1:4)){
    # subset to PRS quartile i
    genomatrix_q <- genomatrix[which(genomatrix$PRS_quartile == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_q[which(genomatrix_q$gen == cn),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(genotype, i, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(HEIGHT ~ 
                  gen +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix_q, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('gen', cn)
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
    
    # get raw distribution
    sumstats_raw <- genomatrix_q %>%
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
    sumstats <- genomatrix_q %>%
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
    out[nrow(out) +1,] <- c(genotype, i,count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
  
}

write.csv(out, 'AnalysisOutput/5_CNV_by_PRSquartile.csv')

##########################
## 6. CNV effect on PRS ##
##########################

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
  model <- lm(PRS ~ 
                gen +
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

write.csv(out, 'AnalysisOutput/6_CNV_effect_on_PRS.csv')

##########################
### 7. PRS by CNV      ###
##########################

out <- data.frame(matrix(nrow=0, ncol=10))
colnames(out) <- c('locus', 'count', 'estimate', 'se', 'p', 'r2', 'r2_adj', 'intercept', 'intercept_se', 'intercept_p')

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
    out[nrow(out) +1,] <- c(genotype, count, NA, NA, NA, NA, NA, NA, NA, NA)
    next
  }
  
  genomatrix_c <- genomatrix[which(genomatrix$gen == cn),]
  
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
  out[nrow(out) +1,] <- c(genotype, count, estimate, se, p, r2, r2_adj, intercept, intercept_se, intercept_p)
}

write.csv(out, 'AnalysisOutput/7_PRS_by_CNV.csv')

############################
## 8. CNVxSex             ##
############################


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
  model <- lm(HEIGHT ~ 
                sex*gen +
                + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
  results <- summary(model)$coefficients
  cov <- paste0('sex2:gen', cn)
  if(!(cov %in% rownames(results))){
    out[nrow(out) +1,] <- c(genotype, count, NA, NA, NA)
    next
  }
  estimate <- as.numeric(results[cov, 1])
  se <- as.numeric(results[cov, 2])
  p <- as.numeric(results[cov, 4])
  
  # write output
  out[nrow(out) +1,] <- c(genotype, count, estimate, se, p)
}

write.csv(out, 'AnalysisOutput/8_CNVxSex_interaction.csv')


############################
## 9. CNV by sex          ##
############################


out <- data.frame(matrix(nrow=0, ncol=20))
colnames(out) <- c('locus', 'sex', 'count', 
                   'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')


for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    
    # split genomatrix
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_s[which(genomatrix_s$gen == cn),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(genotype, sex, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # fit main effects model
    model <- lm(HEIGHT ~ 
                  gen +
                  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix_s, na.action=na.exclude)
    results <- summary(model)$coefficients
    cov <- paste0('gen', cn)
    estimate <- as.numeric(results[cov, 1])
    se <- as.numeric(results[cov, 2])
    p <- as.numeric(results[cov, 4])
    
    
    # get raw distribution
    sumstats_raw <- genomatrix_s %>%
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
    sumstats <- genomatrix_s %>%
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
    out[nrow(out) +1,] <- c(genotype, sex ,count, estimate, se, p, 
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'AnalysisOutput/9_CNV_by_Sex.csv')

######################
## 10. PRSxSex      ##
######################

model <- lm(HEIGHT ~ PRS*sex +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

coeff <- summary(model)$coeff
write.csv(coeff, 'AnalysisOutput/10_PRSxSex.csv')

# write text output to include R2
sink("AnalysisOutput/10_PRSxSex.txt")
summary(summary(model))
sink()

######################
## 11. PRS by Sex   ##
######################

model_male <- lm(HEIGHT ~ PRS +
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix[which(genomatrix$sex == "1"),], na.action=na.exclude)
model_female <- lm(HEIGHT ~ PRS +
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix[which(genomatrix$sex == "2"),], na.action=na.exclude)

coeff_male <- data.frame(summary(model_male)$coeff)
coeff_male$sex <- "male"
coeff_female <- data.frame(summary(model_female)$coeff)
coeff_female$sex <- "female"

coeff <- rbind(coeff_male, coeff_female)
write.csv(coeff, 'AnalysisOutput/11_PRS_by_Sex.csv')

# write text output to include R2
sink("AnalysisOutput/11_PRS_by_Sex.txt")
summary(summary(model_male))
summary(summary(model_female))
sink()


############################
## 12. PRS by CNV by sex  ##
############################

out <- data.frame(matrix(nrow=0, ncol=25))
colnames(out) <- c('locus', 'sex', 'count', 'estimate', 'se', 'p', 'r2', 'r2_adj', 'intercept', 'intercept_se', 'intercept_p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    # count CNV carriers
    count <- nrow(genomatrix_s[which(genomatrix_s$gen == cn),])
    if ((count < 5)){
      out[nrow(out) +1,] <- c(genotype, count, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    genomatrix_c <- genomatrix_s[which(genomatrix_s$gen == cn),]
    
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
    sumstats <- genomatrix_c %>%
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
    out[nrow(out) +1,] <- c(genotype, sex, count, estimate, se, p, r2, r2_adj, intercept, intercept_se, intercept_p,
                            sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                            sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
    
  }
}

write.csv(out, 'AnalysisOutput/12_PRS_by_CNV_by_sex.csv')


###################################
## 13. CNV by sex by PRSquartile ##
###################################

out <- data.frame(matrix(nrow=0, ncol=21))
colnames(out) <- c('locus', 'sex', 'quartile', 'count', 'estimate', 'se', 'p',
                   'mean_raw', 'sd_raw','pct025_raw', 'pct25_raw', 'pct50_raw', 'pct75_raw', 'pct975_raw',
                   'mean', 'sd', 'pct025', 'pct25', 'pct50', 'pct75', 'pct975')

for(c in seq_along(CNVs)){
  # format columns
  CNV <- CNVs[c]
  genotype <- substr(CNV, 7, nchar(CNV))
  locus <- substr(CNV, 1, nchar(CNV) - 4)
  cn <- ifelse(substr(CNV, nchar(CNV) - 2, nchar(CNV)) == 'DEL', '1', '3')
  genomatrix$gen <- genomatrix[[locus]]
  genomatrix$gen <- relevel(as.factor(genomatrix$gen), ref="2")
  
  for(i in c("1", "2")){
    sex <- ifelse(i == "1", "male", "female")
    genomatrix_s <- genomatrix[which(genomatrix$sex == i),]
    
    for(j in 1:4){
      genomatrix_q <- genomatrix_s[which(genomatrix_s$PRS_quartile == j),]
      # count CNV carriers
      count <- nrow(genomatrix_q[which(genomatrix_q$gen == cn),])
      if ((count < 5)){
        out[nrow(out) +1,] <- c(genotype, sex, j, count, 
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        next
      }
      
      # fit PRS model
      model <- lm(HEIGHT ~ 
                    gen +
                    + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = genomatrix_q, na.action=na.exclude)
      results <- summary(model)$coefficients
      cov <- paste0('gen', cn)
      estimate <- as.numeric(results[cov, 1])
      se <- as.numeric(results[cov, 2])
      p <- as.numeric(results[cov, 4])
      
      
      # get raw distribution
      sumstats_raw <- genomatrix_q %>%
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
      sumstats <- genomatrix_q %>%
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
      out[nrow(out) +1,] <- c(genotype, sex, j, count, estimate, se, p,
                              sumstats_raw$mean, sumstats_raw$sd, sumstats_raw$p025, sumstats_raw$p25, sumstats_raw$p50, sumstats_raw$p75, sumstats_raw$p975, 
                              sumstats$mean, sumstats$sd, sumstats$p025, sumstats$p25, sumstats$p50, sumstats$p75, sumstats$p975)
      
    }
  }
}

write.csv(out, 'AnalysisOutput/13_CNV_by_sex_by_PRS_quartile.csv')

############################
## 14. Variance explained ##
############################

dir.create('AnalysisOutput/14_variance_explained')

# change all loci to factors
for(i in seq_along(loci)){
  genomatrix[[loci[i]]] <- relevel(as.factor(genomatrix[[loci[i]]]), ref="2")
}

loci_valid <- loci[sapply(genomatrix[loci], function(x) length(unique(na.omit(x))) > 1)]

# create vector of interaction terms
loci_int <- paste0(loci_valid, ':PRS')

# initialize model
model_init <- lm(HEIGHT ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix, na.action=na.exclude)
# CNVs only

# all loci
model_all_loci <- update(model_init, as.formula(paste("~ . +", paste(loci_valid, collapse = " + "))))

# PRS only
model_PRS <- lm(HEIGHT ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = genomatrix, na.action=na.exclude)

# Additive
model_all_loci_additive <- update(model_all_loci, ~ . + PRS)

# Interaction

# all loci int
model_all_loci_int <- update(model_all_loci_additive, as.formula(paste("~ . +", paste(loci_int, collapse = " + "))))


# Write model summaries

write.csv(data.frame(summary(model_all_loci)$coefficients), 'AnalysisOutput/14_variance_explained/CNV_summary.csv')
write.csv(data.frame(summary(model_all_loci_additive)$coefficients), 'AnalysisOutput/14_variance_explained/CNV_PRS_additive_summary.csv',)
write.csv(data.frame(summary(model_all_loci_int)$coefficients), 'AnalysisOutput/14_variance_explained/CNV_PRS_interaction_summary.csv')
write.csv(data.frame(summary(model_PRS)$coefficients), 'AnalysisOutput/14_variance_explained/PRS_summary.csv')
write.csv(data.frame(summary(model_init)$coefficients), 'AnalysisOutput/14_variance_explained/covariates_summary.csv')

# Write ANOVAs
write.csv(data.frame(anova(model_all_loci)), 'AnalysisOutput/14_variance_explained/CNV_ANOVA.csv')
write.csv(data.frame(anova(model_all_loci_additive)), 'AnalysisOutput/14_variance_explained/CNV_PRS_additive_ANOVA.csv')
write.csv(data.frame(anova(model_all_loci_int)), 'AnalysisOutput/14_variance_explained/CNV_PRS_interaction_ANOVA.csv')
write.csv(data.frame(anova(model_PRS)), 'AnalysisOutput/14_variance_explained/PRS_ANOVA.csv')
write.csv(data.frame(anova(model_init)), 'AnalysisOutput/14_variance_explained/covariates_ANOVA.csv')



#####################################
## 15. Multi-breakpoint analysis   ##
#####################################

out <- data.frame(matrix(nrow=0, ncol=22))
colnames(out) <- c('locusA', 'locusB', 'locusAB', 'genotype',
                   'locusA_count_union', 'locusB_count_union', 'locusA_locusB_count',
                   'locusA_estimate', 'locusA_se', 'locusA_p', 'locusA_additive_VIF',
                   'locusB_estimate', 'locusB_se', 'locus2_p', 'locusB_additive_VIF',
                   'interaction_estimate', 'interaction_se', 'interaction_p',
                   'additive_R2', 'additive_R2_adj', 'int_R2', 'int_R2_adj')

loci_trios <- list(c('locus.22q11.21_A_B', 'locus.22q11.21_B_D', 'locus.22q11.21_A_D'))
# c('locus_1q21.1_TAR_BP1_BP2', 'locus_1q21.1_TAR_BP2_BP3', 'locus_1q21.1_TAR_BP1_BP3'))
# c('locus_16p12.1_BP1_BP2', 'locus_16p12.1_BP2_BP3', 'locus_16p12.1_BP1_BP3'),
# c('locus_16p13.11_BP1_BP2', 'locus_16p13.11_BP2_BP3', 'locus_16p13.11_BP1_BP3'))


for(t in seq_along(loci_trios)){
  trio <- loci_trios[[t]]
  genomatrix$locusA_specific <- relevel(as.factor(genomatrix[[trio[1]]]), ref="2")
  genomatrix$locusB_specific <- relevel(as.factor(genomatrix[[trio[2]]]), ref="2")
  genomatrix$locusAB_specific <- relevel(as.factor(genomatrix[[trio[3]]]), ref="2")
  genomatrix <- genomatrix %>% 
    mutate(locusA = case_when(
      locusA_specific == "1" | locusAB_specific == "1" ~ "1",
      locusA_specific == "3" | locusAB_specific == "3" ~ "3",
      TRUE ~ "2" 
    ))
  genomatrix <- genomatrix %>% 
    mutate(locusB = case_when(
      locusB_specific == "1" | locusAB_specific == "1" ~ "1",
      locusB_specific == "3" | locusAB_specific == "3" ~ "3",
      TRUE ~ "2" 
    ))
  
  genomatrix$locusA <- relevel(as.factor(genomatrix$locusA), ref="2")
  genomatrix$locusB <- relevel(as.factor(genomatrix$locusB), ref="2")
  
  model_additive <- lm(HEIGHT ~ locusA + locusB + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                       data = genomatrix, na.action=na.exclude)
  model_int <- lm(HEIGHT ~ locusA*locusB + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = genomatrix, na.action=na.exclude)
  
  for(cn in c("1", "3")){
    counts <- table(genomatrix$locusA, genomatrix$locusB)
    countA_union <- counts[cn, "2"] + counts[cn, cn]
    countB_union <- counts["2", cn] + counts[cn, cn]
    countAB <- counts[cn, cn]
    
    if(countA_union < 5 | countB_union < 5 | countAB < 5){
      out[nrow(out) + 1,] <- c(trio[1], trio[2], trio[3], cn, 
                               countA_union, countB_union, countAB,
                               NA, NA, NA, NA,
                               NA, NA, NA, NA,
                               NA, NA, NA,
                               NA, NA, NA, NA)
      next
    }
    
    res <- data.frame(summary(model_int)$coefficients)
    estimate_a <- res[paste0('locusA', cn), 1]
    se_a <- res[paste0('locusA', cn), 2]
    p_a <- res[paste0('locusA', cn), 4]
    vif_a <- data.frame(vif(model_additive))["locusA", 1]
    
    estimate_b <- res[paste0('locusB', cn), 1]
    se_b <- res[paste0('locusB', cn), 2]
    p_b <- res[paste0('locusB', cn), 4]
    vif_b <- data.frame(vif(model_additive))["locusB", 1]
    
    estimate_int <- res[paste0('locusA', cn, ':locusB', cn), 1]
    se_int <- res[paste0('locusA', cn, ':locusB', cn), 2]
    p_int <- res[paste0('locusA', cn, ':locusB', cn), 4]
    
    additive_R2 <- summary(model_additive)$r.squared
    additive_R2_adj <- summary(model_additive)$adj.r.squared
    int_R2 <- summary(model_int)$adj.r.squared
    int_R2_adj <- summary(model_int)$adj.r.squared
    
    out[nrow(out) + 1,] <- c(trio[1], trio[2], trio[3], cn, 
                             countA_union, countB_union, countAB,
                             estimate_a, se_a, p_a, vif_a,
                             estimate_b, se_b, p_b, vif_b,
                             estimate_int, se_int, p_int,
                             additive_R2, additive_R2_adj, int_R2, int_R2_adj)
  }
}

write.csv(out, 'AnalysisOutput/15_multibreakpoint_int.csv')

