##############################
### Create Genotype Matrix ###
##############################

# This script demonstrates how phenotypes and genotypes were tabulated

library(MASS)     # for boxcox()
library(forecast) # for BoxCox()
library(dplyr)    # for %>% and mutate()
library(purrr)    # for map_lgl()

####################
### Demographics ###
####################

demographics <- read.csv('path/to/your/demographics')


demographics$sex <- as.factor(demographics$sex, ref='1') # PLINK notation, 1 for male, 2 for female
demographics$age <- as.numeric(demongraphics$age)
demographics$age_sq <- demographics$age ** 2

rownames(demographics) <- demographics$sampleID

genomatrix <- demographics

##################
### Phenotypes ###
##################

pheno <- read.csv('path/to/your/height_bmi_measurements.csv')

rownames(pheno) <- pheno$sampleID

pheno$BMI_raw <- as.numeric(pheno$bmi) #make sure height is in cm
pheno$HEIGHT_raw <- as.numeric(pheno$height) #make sure height is in cm

genomatrix <- cbind(genomatrix, pheno)


# BoxCox scaling by sex

### Scale BMI column ###
# calculate age and sex corrected BMI
genomatrix <- genomatrix[which(!(is.na(genomatrix$BMI_raw)) & !(is.na(genomatrix$sex)) & !(is.na(genomatrix$age))),]
model_resid <- lm(BMI_raw ~ age + age_sq + sex, data=genomatrix, na.action=na.exclude)
genomatrix$bmi_residuals <- residuals(model_resid) + 1 - min(residuals(model_resid)) # add 1 +- (min residual) to ensure all residuals are non-negative for boxcox scaling

# perform BoxCox transformation for males and females separately
genomatrix_male <- genomatrix[which(genomatrix$sex == 1),]
bc <- boxcox(lm(genomatrix_male$bmi_residuals ~ 1))
genomatrix_male$BMI <- BoxCox(genomatrix_male$bmi_residuals, bc$x[which.max(bc$y)])
genomatrix_male$BMI <- scale(genomatrix_male$BMI) # scale to N(0,1)

genomatrix_female <- genomatrix[which(genomatrix$sex == 2),]
bc <- boxcox(lm(genomatrix_female$bmi_residuals ~ 1))
genomatrix_female$BMI <- BoxCox(genomatrix_female$bmi_residuals, bc$x[which.max(bc$y)])
genomatrix_female$BMI <- scale(genomatrix_female$BMI)
genomatrix <- rbind(genomatrix_male, genomatrix_female) # scale to N(0,1)


### Scale HEIGHT column ###
# calculate age and sex corrected HEIGHT
genomatrix <- genomatrix[which(!(is.na(genomatrix$HEIGHT_raw)) & !(is.na(genomatrix$sex)) & !(is.na(genomatrix$age))),]
model_resid <- lm(HEIGHT_raw ~ age + age_sq + sex, data=genomatrix, na.action=na.exclude)
genomatrix$height_residuals <- residuals(model_resid) + 1 - min(residuals(model_resid)) # add 1 +- (min residual) to ensure all residuals are non-negative for boxcox scaling

# perform BoxCox transformation for males and females separately
genomatrix_male <- genomatrix[which(genomatrix$sex == 1),]
bc <- boxcox(lm(genomatrix_male$height_residuals ~ 1))
genomatrix_male$HEIGHT <- BoxCox(genomatrix_male$height_residuals, bc$x[which.max(bc$y)])
genomatrix_male$HEIGHT <- scale(genomatrix_male$HEIGHT) # scale to N(0,1)

genomatrix_female <- genomatrix[which(genomatrix$sex == 2),]
bc <- boxcox(lm(genomatrix_female$height_residuals ~ 1))
genomatrix_female$HEIGHT <- BoxCox(genomatrix_female$height_residuals, bc$x[which.max(bc$y)])
genomatrix_female$HEIGHT <- scale(genomatrix_female$HEIGHT)
genomatrix <- rbind(genomatrix_male, genomatrix_female) # scale to N(0,1)





###########
### PCs ###
###########

pcs <- read.csv('path/to/your/ancestry_pcs.csv')
rownames(pcs) <- pcs$sampleID

# merge PCs with genomatrix
genomatrix <- merge(genomatrix, pcs[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')],
                    by = 'row.names', all.x = TRUE)
rownames(genomatrix) <- genomatrix$Row.names
genomatrix$Row.names <- NULL


##############
### PGSs  ###
##############

pgs <- read.csv('path/to/your/pgs_scores.csv')
rownames(pgs) <- pgs$sampleID

# merge PGS with genomatrix and scale to N(0,1)
genomatrix <- merge(genomatrix, pgs[, c('PGS_HEIGHT', 'PGS_BMI')],
                    by = 'row.names', all.x = TRUE)
rownames(genomatrix) <- genomatrix$Row.names
genomatrix$Row.names <- NULL

genomatrix$PGS_HEIGHT <- scale(genomatrix$PGS_HEIGHT) # scale to N(0,1)
genomatrix$PGS_BMI <- scale(genomatrix$PGS_BMI) # scale to N(0,1)



###################
### Medications ###
###################

# Add medication information
# I.E. A table of a list of medication names for each sampleID (column is called 'meds_names').
medications <- read.csv('path/to/your/medications.csv')
rownames(medications) <- medications$sampleID
genomatrix <- cbind(genomatrix, medications)

# This code matches medication names to categories:
# mood_stabilizer, antipsychotic, ssri, snri, tca, maoi, and atypical

medication_index <- read.csv('../BMI/PsychiatricMedicationsRXNorm.csv') 

mood_stabilizer_ids <- unique(medication_index[which(medication_index$drug_class == 'mood_stabilizer'), 'drug_name'])
antipsychotic_ids <- unique(medication_index[which(medication_index$drug_class == 'antipsychotic'), 'drug_name'])
ssri_ids <- unique(medication_index[which(medication_index$drug_class == 'ssri'), 'drug_name'])
snri_ids <- unique(medication_index[which(medication_index$drug_class == 'snri'), 'drug_name'])
tca_ids <- unique(medication_index[which(medication_index$drug_class == 'tca'), 'drug_name'])
maoi_ids <- unique(medication_index[which(medication_index$drug_class == 'maoi'), 'drug_name'])
atypical_ids <- unique(medication_index[which(medication_index$drug_class == 'atypical'), 'drug_name'])

### use this if you have RXNormIDs
# mood_stabilizer_ids <- unique(medication_index[which(medication_index$drug_class == 'mood_stabilizer'), 'rxcui'])
# antipsychotic_ids <- unique(medication_index[which(medication_index$drug_class == 'antipsychotic'), 'rxcui'])
# ssri_ids <- unique(medication_index[which(medication_index$drug_class == 'ssri'), 'rxcui'])
# snri_ids <- unique(medication_index[which(medication_index$drug_class == 'snri'), 'rxcui'])
# tca_ids <- unique(medication_index[which(medication_index$drug_class == 'tca'), 'rxcui'])
# maoi_ids <- unique(medication_index[which(medication_index$drug_class == 'maoi'), 'rxcui'])
# atypical_ids <- unique(medication_index[which(medication_index$drug_class == 'atypical'), 'rxcui'])



medication_df <- read.csv('path/to/your/medication/table.csv') # .csv where medication info is stored
rownames(medication_df) <- medication_df$IID # rename rows to sampleID
medication_df <- medication_df[genomatrix$sampleID,] # match rows to genomatrix
medication_df$meds_names <- medication_df$meds_names # whatever your column of lists is called

# convert string to list if necessary
# medication_df <- medication_df %>%
#   mutate(meds_names = strsplit(substr(meds_names, 3, nchar(meds_names) - 2), "', '"))

medication_df <- medication_df %>% 
  mutate(mood_stabilizer = as.integer(map_lgl(meds_names, ~ any(.x %in% mood_stabilizer_ids))))
medication_df <- medication_df %>% 
  mutate(antipsychotic = as.integer(map_lgl(meds_names, ~ any(.x %in% antipsychotic_ids))))
medication_df <- medication_df %>% 
  mutate(ssri = as.integer(map_lgl(meds_names, ~ any(.x %in% ssri_ids))))
medication_df <- medication_df %>% 
  mutate(snri = as.integer(map_lgl(meds_names, ~ any(.x %in% snri_ids))))
medication_df <- medication_df %>% 
  mutate(tca = as.integer(map_lgl(meds_names, ~ any(.x %in% tca_ids))))
medication_df <- medication_df %>% 
  mutate(maoi = as.integer(map_lgl(meds_names, ~ any(.x %in% maoi_ids))))
medication_df <- medication_df %>% 
  mutate(atypical = as.integer(map_lgl(meds_names, ~ any(.x %in% atypical_ids))))


genomatrix$ssri <- medication_df$ssri 
genomatrix$snri <- medication_df$snri 
genomatrix$tca <- medication_df$tca 
genomatrix$maoi <- medication_df$maoi
genomatrix$atypical <- medication_df$atypical
genomatrix$mood_stabilizer <- medication_df$mood_stabilizer # add mood stabilizers
genomatrix$antipsychotic <- medication_df$antipsychotic # add anti-psychotics
genomatrix$antidepressant <- ifelse(rowSums(genomatrix[, c("ssri", "snri", "tca", "maoi", "atypical")]) > 0, 1, 0)
genomatrix$all_medication <- ifelse(rowSums(genomatrix[, c("antidepressant", "antipsychotic", "mood_stabilizer")]) > 0, 1, 0)


