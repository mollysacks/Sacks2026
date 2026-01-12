##############################
### Psychiatric Meds Plots ###
##############################

setwd("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByMolly/MetaAnalysis")
library("ggplot2")
library("dplyr")
library(data.table)
library(ggpmisc)
library(ggrepel)
library(ggpubr)
library(GenomicRanges)
library(biomaRt)
library(ggh4x)
library(ggpattern)
library(tidyverse)
library(dplyr)
library(scales)
library(stringr)


change_CNV_names <- function(CNV_main_effect){
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_A_D", " A-D", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_B_D", " B-D", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_C_D", " C-D", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_A_B", " A-B", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_E_F", " E-F", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_F_H", " F-H", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP4_BP5", " BP4-BP5", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP1_BP4", " BP1-BP4", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP1_BP3", " BP1-BP3", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP2_BP4", " BP2-BP4", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP1_BP2", " BP1-BP2", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP1_BP3", " BP1-BP3", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_BP2_BP3", " BP2-BP3", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_proximal", " proximal", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_distal", " distal", locus))
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_WBS", " WBS", locus))
  
  
  CNV_main_effect <- CNV_main_effect  %>% 
    mutate(locus = gsub("_D", " D", locus))
  
  return(CNV_main_effect)
}




########################################
### Chi2 meds/ CNV                   ###
########################################



medications <- c("antidepressant", "antipsychotic", "mood_stabilizer")
lifetime_use <- c("MyCode.EUR", 
                  "MVP.EUR", "MVP.AFR", 'MVP.LATNAT',
                  "AoU.EUR", "AoU.AAM", "AoU.LATNAT")
current_use <- c("UKBB.EUR", "EstBB")
options(digits = 13)

Medication_use_table <- read.csv("PsychiatricMedication.v2/tables/18_CNV_causal_mediation.csv")
Medication_use_table <- Medication_use_table[which(Medication_use_table$medication %in% medications),]
Medication_use_table <- Medication_use_table[which(Medication_use_table$cohort %in% lifetime_use),]
Medication_use_table <- Medication_use_table[which(Medication_use_table$Effect == "pure direct effect"),]
Medication_use_table$use_pattern <- "current_use"
Medication_use_table <- Medication_use_table[,c("CNV", "medication", "nCNV_med",
                                                "nCNV_nomed", "n_noCNV_med", "n_noCNV_nomed", "cohort", "use_pattern"),]
#Medication_use_table$CNV_group <- sub("(.+?)(mood|anti).*", "\\1", Medication_use_table$group)
#Medication_use_table$CNV_group <- "CNV"
#Medication_use_table$usage_type <- ifelse(Medication_use_table$cohort %in% lifetime_use, "lifetime_use", "current_use")


colnames(Medication_use_table) <- c("group", "Medication_type", "A", "B", "C", "D", "cohort", "use_pattern")

Medication_use_table <- Medication_use_table %>%
  mutate(Medication_type = gsub("_", " ", Medication_type)) 

Medication_use_table <- Medication_use_table %>%
  group_by(Medication_type, group) %>%
  summarise(
    A = sum(A, na.rm = TRUE),
    B = sum(B, na.rm = TRUE),
    C = sum(C, na.rm = TRUE),
    D = sum(D, na.rm = TRUE),
    .groups = "drop"
  )

# Medication_use_table <- Medication_use_table %>%
#   group_by(Medication_type) %>%
#   summarise(
#     A = sum(A, na.rm = TRUE),
#     B = sum(B, na.rm = TRUE),
#     C = dplyr::first(C),
#     D = dplyr::first(D),
#     .groups = "drop"
#   )

Medication_use_table$chisq_p_value <- sapply(1:nrow(Medication_use_table), function(i) {
  fisher.test(matrix(c(Medication_use_table$A[i], Medication_use_table$B[i], 
                       Medication_use_table$C[i], Medication_use_table$D[i]), nrow = 2))$p.value
})

# Force scientific notation for small p-values
Medication_use_table$chisq_p_value <- format(Medication_use_table$chisq_p_value, scientific = TRUE)
#Medication_use_table <- Medication_use_table[which((Medication_use_table$CNV_group == "no effect") & (Medication_use_table$cohort %in% c('MVP', 'UKBB'))),]

plot_data <- data.frame(matrix(nrow=0, ncol=4))
colnames(plot_data) <- c('medication', 'percent', 'p', 'CNV_group')


for(i in 1:nrow(Medication_use_table)){
  row <- Medication_use_table[i,]
  new_row1 <- c(row$Medication_type, row$A / (row$B + row$A), row$chisq_p_value, row$group)
  plot_data[nrow(plot_data) + 1,] <- new_row1
  
  new_row2 <- c(row$Medication_type, row$C / (row$C + row$D), 1, 'No CNV')
  plot_data[nrow(plot_data) + 1,] <- new_row2
}

plot_data$sig <- ifelse(as.numeric(plot_data$p) < 0.05, ifelse(as.numeric(plot_data$p) < 0.001, '**', '*'), '')


plot_data$CNV_group <- gsub('group_', '', plot_data$CNV_group)
plot_data$CNV_group <- gsub('____', '', plot_data$CNV_group)

custom_colors <- c(
  "CNV"="goldenrod",
  "No CNV"="grey30"
)

custom_colors <- c(
  "positive" = "#D95F02",
  "negative" = "#1B9E77",
  "neutral"="goldenrod",
  "no CNV"="grey30"
  # Add more categories and colors as needed
)

custom_pattern <- c(
  "none" = "none",
  "antidepressant" = "crosshatch",
  "antipsychotic" = "stripe",
  "mood_stabilizer" = "circle"
)

plot_data$percent <- as.numeric(plot_data$percent)
plot_data$medication <- factor(plot_data$medication, levels=c('antidepressant', 'mood stabilizer', 'antipsychotic'))

bar_plot <- ggplot(plot_data, aes(x = medication, y = percent, fill = CNV_group, color=CNV_group)) +
  geom_bar(stat = "identity", position=position_dodge(width = 0.9)) +  # "fill" makes it proportional
  scale_y_continuous(limits = c(0, 1),labels = scales::percent_format()) +  # Show percentages
  scale_fill_manual(values=custom_colors, guide="legend") +
  scale_color_manual(values=custom_colors, guide="none") +
  geom_text(aes(label=sig, y=percent + .05), position =position_dodge(width = 0.9), color='black', size=15) +
  # Merge into a single legend
  labs(y = "CNV group", x = "Proportion", fill = "CNV group") +
  #facet_wrap(~ CNV_group, nrow = 1) +
  # geom_text(aes(label = sig, y= percent),
  #           position = position_fill(vjust = 0.9),  # Center within segments
  #           size = 10, color = "black" )+
  theme_classic() +
  theme(strip.text = element_text(hjust = 0.5, size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size=25),
        legend.title = element_blank(),
        legend.position = "bottom")

bar_plot

ggsave('Plots.v2/PsychiatricMeds/random_effects/meds_chi2.png', bar_plot, height=4, width=10)




########################################
### CNV main effect meds factor      ###
########################################



########################################
### Med effect by CNV group          ###
########################################


drug_by_CNV <- read.csv('PsychiatricMedication.v2/random_effects/3_med_by_CNV.csv')
drug_main_effect <- read.csv('PsychiatricMedication.v2/random_effects/1_med_main_effect.csv')
drug_main_effect$locus <- "no CNV"
drug_by_CNV <- rbind(drug_by_CNV, drug_main_effect)
drug_by_CNV <- drug_by_CNV[which(drug_by_CNV$medication != 'maoi'),]
drug_by_CNV$sig <- ifelse(drug_by_CNV$p < 0.05, 
                          ifelse(drug_by_CNV$p < 0.001, '**', '*'), '')
drug_by_CNV <- drug_by_CNV[which(drug_by_CNV$medication %in% c('all_medication', 'antidepressant', 'mood_stabilizer', 'antipsychotic')),]
drug_by_CNV$rowname <- paste0(drug_by_CNV$group, drug_by_CNV$medication)


drugxCNV <- read.csv('PsychiatricMedication.v2/random_effects/2_CNVxmed.csv')
drugxCNV$sig <- ifelse(drugxCNV$p < 0.05, 
                       ifelse(drugxCNV$p < 0.001, '**', '*'), '')
rownames(drugxCNV) <- drugxCNV$group
# 
# drug_by_CNV$sig_int <- drugxCNV[drug_by_CNV$rowname, 'sig']


custom_colors <- c(
  "positive CNV" = "#D95F02",
  "negative CNV" = "#1B9E77",
  "neutral CNV"="goldenrod",
  "no CNV"="grey30"
  # Add more categories and colors as needed
)



plot_Medication_by_CNV <- function(effects){
  print(effects$group)
  custom_colors <- c(
    "positive" = "#D95F02",
    "negative" = "#1B9E77",
    "neutral"="goldenrod",
    "no CNV"="grey30"
    # Add more categories and colors as needed
  )
  
  plot <- ggplot(effects, aes(x=medication, y=estimate, ymin=ci_lower, ymax=ci_upper, color=group, fill=group)) +
    geom_hline(yintercept = 0) +
    geom_linerange(size=6, alpha=.5, position=position_dodge(width = 0.8)) +
    geom_point(size=8, position=position_dodge(width = 0.8)) +
    geom_text(aes(label=sig), position=position_dodge(width = 0.8), color="white", vjust = 0.8, size=6) +
    #geom_text(aes(y=0.42, label=sig_int), position=position_dodge(width = 0.8), color="black", vjust = 0.8, size=10) +
    theme_classic() +
    scale_color_manual(values=custom_colors) +
    scale_fill_manual(values=custom_colors) +
    #annotate("text", x = Inf, y = -Inf, label = paste("* = p interaction < 0.05"), hjust = 1.1, vjust = -.5, size = 8) +
    # theme(plot.title = element_text(hjust = 0.5, size=25), 
    #       axis.text.x = element_text(size=25),
    #       axis.text.y = element_text(size=25),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_text(size=25),
    #       legend.text = element_text(size=25),
    #       legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size=20), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_blank()) +
    labs(y='Drug effect on BMI (s.d.)')
  return(plot)
}

drug_by_CNV <- drug_by_CNV[which(drug_by_CNV$locus %in% 
                                   c('group_negative____', 'group_neutral____', 'group_positive____', 'no CNV')), ]

drug_by_CNV$group <- drug_by_CNV$locus

drug_by_CNV <- drug_by_CNV %>%
  mutate(medication = gsub("_", " ", medication)) 

drug_by_CNV <- drug_by_CNV %>%
  mutate(group = gsub("group_", "", group)) 

drug_by_CNV <- drug_by_CNV %>%
  mutate(group = gsub("____", "", group)) 

drug_by_CNV <- drug_by_CNV[which(drug_by_CNV$group %in% c("negative", "no CNV", "neutral", "positive")),]
drug_by_CNV$group <- factor(drug_by_CNV$group, levels=c("negative", "no CNV", "neutral", "positive"))

drug_by_CNV$medication <- ifelse(drug_by_CNV$medication == 'all medication', 'combined', drug_by_CNV$medication)
drug_by_CNV$medication <- factor(drug_by_CNV$medication, levels=c('antidepressant', 'mood stabilizer', 'antipsychotic', 'combined'))

plot <- plot_Medication_by_CNV(drug_by_CNV)
print(plot)
#ggsave(paste0('Plots/BMI/random_effects/Medication_by_CNV_group.png'), plot, height=5, width=15)
ggsave(paste0('Plots.v2/PsychiatricMeds/random_effects/Medication_by_CNV_group.png'), plot, height=5, width=11)

####################################
### Medication Main effect       ###
####################################


drug_by_CNV <- read.csv('PsychiatricMedication.v2/tables/1_med_main_effect.csv')
drug_by_CNV <- drug_by_CNV[which(drug_by_CNV$medication %in% c('all_medication', 'antidepressant', 'mood_stabilizer', 'antipsychotic')),]
drug_by_CNV$sig <- ifelse(drug_by_CNV$p < 0.05, 
                          ifelse(drug_by_CNV$p < 0.001, '**', '*'), '')

drug_by_CNV$cohort <- factor(drug_by_CNV$cohort, levels=c('EstBB', 'MyCode.EUR', 
                                                          'MVP.EUR', 'MVP.LATNAT', 'MVP.AFR',
                                                          'UKBB.EUR', 
                                                          'AoU.EUR', 'AoU.AAM', 'AoU.LATNAT'))

med_main_effects <- drug_by_CNV

med_main_effects <- med_main_effects %>%
  mutate(medication = gsub("_", " ", medication))

med_main_effects$medication <- ifelse(med_main_effects$medication == 'all medication', 'combined', med_main_effects$medication)
med_main_effects$medication <- factor(med_main_effects$medication, levels=c('antidepressant', 'mood stabilizer', 'antipsychotic', 'combined'))

med_main_effects$ci_lower <- med_main_effects$estimate - (1.96 * med_main_effects$se)
med_main_effects$ci_upper <- med_main_effects$estimate + (1.96 * med_main_effects$se)

plot_med_main_effects <- ggplot(med_main_effects, aes(x=medication, y=estimate, ymin=ci_lower, ymax=ci_upper)) +
  geom_linerange(size=10, alpha=.5, position=position_dodge(width = 0.1), color='black', fill='black') +
  geom_point(size=12, color='black', fill='black', position=position_dodge(width = 0.1)) +
  geom_text(aes(label=sig), color='white', fill='white', position=position_dodge(width = 0.1), color="white", hjust = 0.5, vjust=0.8, size=10) +
  theme_classic() +
  scale_color_manual(values=custom_colors) +
  scale_fill_manual(values=custom_colors) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~cohort, nrow = 3, ncol = 3) +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.x = element_text(size=25, angle=30, vjust=1, hjust=1),
        axis.text.y = element_text(size=25),
        strip.text = element_text(size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=25),
        legend.position = "none") +
  labs(y='Drug effect on BMI (s.d.)')

plot_med_main_effects

ggsave(paste0('Plots.v2/PsychiatricMeds/random_effects/Medication_main_effects.png'), plot_med_main_effects, height=15, width=15)



########################################
### CNV effect (meds vs. no meds)    ###
########################################

CNV_by_meds <- read.csv('PsychiatricMedication.v2/random_effects/4_CNV_by_med.csv')

CNV_by_meds <- CNV_by_meds[which(CNV_by_meds$medication %in% c("no_meds", "all_medication")),]
                           
CNV_by_meds <- CNV_by_meds[which(!(CNV_by_meds$locus %in% 
                                     c("group_negative____", "group_neutral____", "group_positive____"))),]                   

CNV_by_meds <- change_CNV_names(CNV_by_meds)


CNV_by_meds_no_meds <- CNV_by_meds[which(CNV_by_meds$medication == "no_meds"),]
CNV_by_meds_meds <- CNV_by_meds[which(CNV_by_meds$medication == "all_medication"),]

rownames(CNV_by_meds_no_meds) <- CNV_by_meds_no_meds$locus
rownames(CNV_by_meds_meds) <- CNV_by_meds_meds$locus

CNV_by_meds_meds <- CNV_by_meds_meds[rownames(CNV_by_meds_no_meds),]

plot_data <- data.frame(cbind(CNV_by_meds_no_meds$locus, 
                   CNV_by_meds_no_meds$estimate, CNV_by_meds_no_meds$p,
                   CNV_by_meds_meds$estimate, CNV_by_meds_meds$p))

colnames(plot_data) <- c("locus", "estimate_noMeds", "p_noMeds",
                         "estimate_Meds", "p_Meds")

plot_data$estimate_Meds <- as.numeric(plot_data$estimate_Meds)
plot_data$estimate_noMeds <- as.numeric(plot_data$estimate_noMeds)



plot_data <- plot_data[which(plot_data$p_noMeds < 0.05),]


scatter <- ggplot(plot_data, aes(x = estimate_noMeds, y = estimate_Meds)) +
  geom_point(aes(size=6), color="black", fill="black", shape=21, stroke = 1.5) +
  geom_smooth(show.legend = FALSE, method = lm, se=TRUE, data=plot_data, fullrange=TRUE, color="darkgrey", fill="darkgrey") +
  geom_hline(yintercept=0, linetype='dashed', color = "darkgrey") +
  geom_vline(xintercept=0, linetype='dashed', color = "darkgrey") +
  geom_abline(yintercept=0, slope=1, linetype='dashed', color = "red") +
  labs(x = "Effect Size in non-users (BMI s.d.)",
       y = "Effect Size in medication users (BMI s.d.)") +
  geom_text_repel(aes(label = locus, color="black"), color="black", size = 5, show.legend = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 25),
        legend.position=c(.85,.1),
        plot.title = element_text(hjust = 0.5)) +
  guides(
    size = "none",  # Remove size from the legend
    fill = guide_legend(override.aes = list(size = 5)),  # Increase the size of the legend icons
    color = guide_legend(override.aes = list(size = 5))  # Increase the size of the legend icons
  ) 



ggsave('Plots.v2/PsychiatricMeds/random_effects/CNV_main_effects.meds_vs_nomeds.png', scatter, height=10, width=10)


########################################
### Med effect by PRS quartile.      ###
########################################

drug_by_PRS <- read.csv('PsychiatricMedication.v2/random_effects/6_Med_by_PRS.csv')
drug_by_PRS <- drug_by_PRS[which(drug_by_PRS$medication != 'maoi'),]
drug_by_PRS$sig <- ifelse(drug_by_PRS$p < 0.05, 
                          ifelse(drug_by_PRS$p < 0.001, '**', '*'), '')
drug_by_PRS <- drug_by_PRS[which(drug_by_PRS$medication %in% c('all_medication', 'antidepressant', 'mood_stabilizer', 'antipsychotic')),]
drug_by_PRS$quartile <- as.character(drug_by_PRS$quartile)
rownames(drug_by_PRS) <- paste0(drug_by_PRS$quartile, drug_by_PRS$medication)

drugXPRS <- read.csv('PsychiatricMedication.v2/random_effects/5_MedXPRS.csv')
rownames(drugXPRS) <- paste0('2', drugXPRS$medication)
drugXPRS <- drugXPRS[rownames(drug_by_PRS),]

drug_by_PRS <- drug_by_PRS %>% 
  mutate(medication = gsub("_", " ", medication))

drug_by_PRS$medication <- ifelse(drug_by_PRS$medication == 'all medication', 'combined', drug_by_PRS$medication)
drug_by_PRS$medication <- factor(drug_by_PRS$medication, levels=c('antidepressant', 'mood stabilizer', 'antipsychotic', 'combined'))

drug_by_PRS$sig_int <- ifelse(drugXPRS$p < 0.05, paste0('PRSxMedication \np =', substr(drugXPRS$p, 1,3), substr(drugXPRS$p, nchar(drugXPRS$p) - 4, nchar(drugXPRS$p))), NA)
custom_colors <- c(
  "1" = "#85a9c9",
  "2" = "#466f94",
  "3"="#1c4366",
  "4"="#082640"
  # Add more categories and colors as needed
)


plot <- ggplot(drug_by_PRS, aes(x=medication, y=estimate, ymin=ci_lower, ymax=ci_upper, color=quartile, fill=quartile)) +
  geom_hline(yintercept = 0) +
  geom_linerange(size=10, alpha=.6, position=position_dodge(width = 0.8)) +
  geom_point(size=12, position=position_dodge(width = 0.8)) +
  geom_text(aes(label=sig), position=position_dodge(width = 0.8), color="white", vjust = 0.8, size=10) +
  #geom_text(aes(label=sig_int, y= 0.4), position=position_dodge(width = 0.8), color="black", vjust = 0.8, size=6) +
  scale_color_manual(values=custom_colors) +
  scale_fill_manual(values=custom_colors) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=25), 
        axis.text.x = element_text(size=25, vjust=.7),
        axis.text.y = element_text(size=25),
        axis.title = element_text(size=25),
        axis.title.x = element_blank(),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25)) +
  labs(x='PGS-BMI quartile', y='Drug effect on BMI (s.d.)', color='PGS-BMI\nquartile', fill='PGS-BMI\nquartile')

print(plot)
ggsave(paste0('Plots.v2/PsychiatricMeds/random_effects/Med_by_PRS.png'), plot,  height=5, width=15)



################################################
###  Proportion on meds vs CNV main effect   ###
################################################


main_effects <-  read.csv('BMI.v2/fixed_effects/2_CNV_main_effect.csv')
main_effects <- change_CNV_names(main_effects)
rownames(main_effects) <- main_effects$locus

props <- read.csv('PsychiatricMedication.v2/random_effects/3_med_by_CNV.csv')
props <- change_CNV_names(props)
props$CNV_main <- main_effects[props$locus,'estimate']
props$CNV_main_p <- main_effects[props$locus, 'p']
props$count_total <- main_effects[props$locus, 'count']
props$prop <- props$count/ props$count_total

props <- props[which(props$CNV_main_p < 0.05),]

props <- props[which(props$prop < 0.7),]

scatter <- ggplot(props, aes(x = prop, y = CNV_main)) +
  geom_point(aes(size=6), color="black", fill="black", shape=21, stroke = 1.5) +
  geom_smooth(show.legend = FALSE, method = lm, se=TRUE, data=props, fullrange=TRUE, color="darkgrey", fill="darkgrey") +
  geom_hline(yintercept=0, linetype='dashed', color = "darkgrey") +
  geom_vline(xintercept=0, linetype='dashed', color = "darkgrey") +
  labs(x = "Proportion on medication",
       y = "Effect Size (BMI s.d.)") +
  geom_text_repel(aes(label = locus, color="black"), color="black", size = 5, show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~medication,nrow = 2, ncol=2) +
  theme(text = element_text(size = 25),
        legend.position=c(.85,.1),
        plot.title = element_text(hjust = 0.5)) +
  guides(
    size = "none",  # Remove size from the legend
    fill = guide_legend(override.aes = list(size = 5)),  # Increase the size of the legend icons
    color = guide_legend(override.aes = list(size = 5))  # Increase the size of the legend icons
  ) 

print(scatter)


ggsave('Plots.v2/PsychiatricMeds/random_effects/CNV_main_effects.vs.proportion_on_meds.png', scatter, height=20, width=20)


################################################
###  CNV main effect BMI vs meds factor      ###
################################################

# Read CNV main effect on BMI
CNV_main_BMI <- read.csv('BMI.v2/random_effects/2_CNV_main_effect.csv')
CNV_main_BMI <- change_CNV_names(CNV_main_BMI)
CNV_main_BMI <- CNV_main_BMI[which(!(CNV_main_BMI$locus %in%
                                       c("group_negative____", "group_neutral____", "group_positive____"))),]
rownames(CNV_main_BMI) <- CNV_main_BMI$locus

# Read CNV main effect meds factor
CNV_meds_factor <- read.csv('PsychiatricMedication.v2/random_effects/16_CNV_main_effect_meds_factor.csv')
CNV_meds_factor <- change_CNV_names(CNV_meds_factor)
CNV_meds_factor <- CNV_meds_factor[which(!(CNV_meds_factor$locus %in%
                                             c("group_negative____", "group_neutral____", "group_positive____"))),]
rownames(CNV_meds_factor) <- CNV_meds_factor$locus

# Get common loci
common_loci <- intersect(rownames(CNV_main_BMI), rownames(CNV_meds_factor))
CNV_main_BMI <- CNV_main_BMI[common_loci,]
CNV_meds_factor <- CNV_meds_factor[common_loci,]

# Create plot data
plot_data_BMI_meds <- data.frame(
  locus = CNV_main_BMI$locus,
  estimate_BMI = CNV_main_BMI$estimate,
  p_BMI = CNV_main_BMI$p,
  estimate_meds_factor = CNV_meds_factor$estimate,
  p_meds_factor = CNV_meds_factor$p
)

# Filter to significant CNVs (based on BMI main effect)
plot_data_BMI_meds <- plot_data_BMI_meds[which(plot_data_BMI_meds$p_BMI < 0.05),]

scatter_BMI_meds <- ggplot(plot_data_BMI_meds, aes(x = estimate_BMI, y = estimate_meds_factor)) +
  geom_point(aes(size=6), color="black", fill="black", shape=21, stroke = 1.5) +
  geom_smooth(show.legend = FALSE, method = lm, se=TRUE, data=plot_data_BMI_meds, fullrange=TRUE, color="darkgrey", fill="darkgrey") +
  geom_hline(yintercept=0, linetype='dashed', color = "darkgrey") +
  geom_vline(xintercept=0, linetype='dashed', color = "darkgrey") +
  geom_abline(yintercept=0, slope=1, linetype='dashed', color = "red") +
  labs(x = "CNV main effect on BMI (BMI s.d.)",
       y = "CNV main effect meds factor (BMI s.d.)") +
  geom_text_repel(aes(label = locus, color="black"), color="black", size = 5, show.legend = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 25),
        legend.position=c(.85,.1),
        plot.title = element_text(hjust = 0.5)) +
  guides(
    size = "none",  # Remove size from the legend
    fill = guide_legend(override.aes = list(size = 5)),  # Increase the size of the legend icons
    color = guide_legend(override.aes = list(size = 5))  # Increase the size of the legend icons
  )

print(scatter_BMI_meds)

ggsave('Plots.v2/PsychiatricMeds/random_effects/CNV_main_effects.BMI_vs_meds_factor.png', scatter_BMI_meds, height=10, width=10)


################################################
###  Causal Mediation Analysis Plots         ###
################################################

# Read causal mediation data
causal_mediation <- read.csv('PsychiatricMedication.v2/random_effects/18_CNV_causal_mediation.csv')

# Filter to only direct effect and causal mediation effect (ignore interaction terms)
causal_mediation <- causal_mediation[which(causal_mediation$Effect %in%
                                             c("pure direct effect", "pure causal mediation effect")),]

# Filter to CNV groups
causal_mediation <- causal_mediation[which(causal_mediation$CNV %in%
                                             c("group_negative____", "group_neutral____", "group_positive____")),]

# Clean up names
causal_mediation$CNV_group <- gsub("group_", "", causal_mediation$CNV)
causal_mediation$CNV_group <- gsub("____", "", causal_mediation$CNV_group)
causal_mediation$CNV_group <- factor(causal_mediation$CNV_group, levels = c("negative", "neutral", "positive"))

causal_mediation$medication <- gsub("_", " ", causal_mediation$medication)
causal_mediation$medication <- ifelse(causal_mediation$medication == "all medication", "combined", causal_mediation$medication)

# Rename effects for legend
causal_mediation$Effect <- gsub("pure direct effect", "direct effect", causal_mediation$Effect)
causal_mediation$Effect <- gsub("pure causal mediation effect", "causal mediation effect", causal_mediation$Effect)
causal_mediation$Effect <- factor(causal_mediation$Effect, levels = c("direct effect", "causal mediation effect"))

# Add significance stars
causal_mediation$sig <- ifelse(causal_mediation$p < 0.05,
                                ifelse(causal_mediation$p < 0.001, '**', '*'), '')

# Colors matching the example plot
mediation_colors <- c("direct effect" = "#5B9BD5", "causal mediation effect" = "#BFA56A")

# Create a plot for each medication
medications <- c("antidepressant", "antipsychotic", "mood stabilizer", "combined")

for(med in medications) {
  plot_data_med <- causal_mediation[which(causal_mediation$medication == med),]

  med_plot <- ggplot(plot_data_med, aes(x = CNV_group, y = estimate, fill = Effect)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sig, y = ifelse(estimate >= 0, ci_upper + 0.01, ci_lower - 0.01)),
              position = position_dodge(width = 0.8), size = 8, vjust = ifelse(plot_data_med$estimate >= 0, 0, 1)) +
    scale_fill_manual(values = mediation_colors) +
    labs(title = paste0("Mediation of CNV effect on BMI by ", med),
         x = "CNV group",
         y = "Effect on BMI (s.d.)",
         fill = "") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.position = "bottom"
    )

  print(med_plot)

  # Save plot
  filename <- paste0('Plots.v2/PsychiatricMeds/random_effects/mediation_analysis_',
                     gsub(" ", "_", med), '.png')
  ggsave(filename, med_plot, height = 6, width = 8)
}


################################################
###  PRS Causal Mediation Analysis Plot      ###
################################################

# Read PRS causal mediation data (per-cohort)
prs_mediation <- read.csv('PsychiatricMedication.v2/tables/19_PRS_causal_mediation.csv')

# Filter to only direct effect and causal mediation effect (ignore interaction terms)
prs_mediation <- prs_mediation[which(prs_mediation$Effect %in%
                                       c("pure direct effect", "pure causal mediation effect")),]

# Meta-analyze across cohorts using inverse variance weighting
prs_meta <- prs_mediation %>%
  group_by(medication, Effect) %>%
  summarise(
    # Inverse variance weighted mean
    estimate = sum(estimate / se^2) / sum(1 / se^2),
    # Combined standard error
    se = sqrt(1 / sum(1 / se^2)),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    z = estimate / se,
    p = 2 * pnorm(-abs(z))
  )

# Clean up names
prs_meta$medication <- gsub("_", " ", prs_meta$medication)
prs_meta$medication <- ifelse(prs_meta$medication == "all medication", "combined", prs_meta$medication)
prs_meta$medication <- factor(prs_meta$medication, levels = c("antidepressant", "antipsychotic", "mood stabilizer", "combined"))

# Rename effects for legend
prs_meta$Effect <- gsub("pure direct effect", "direct effect", prs_meta$Effect)
prs_meta$Effect <- gsub("pure causal mediation effect", "causal mediation effect", prs_meta$Effect)
prs_meta$Effect <- factor(prs_meta$Effect, levels = c("direct effect", "causal mediation effect"))

# Add significance stars
prs_meta$sig <- ifelse(prs_meta$p < 0.05,
                        ifelse(prs_meta$p < 0.001, '**', '*'), '')

# Create the plot
prs_med_plot <- ggplot(prs_meta, aes(x = medication, y = estimate, fill = Effect)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sig, y = ifelse(estimate >= 0, ci_upper + 0.005, ci_lower - 0.005)),
            position = position_dodge(width = 0.8), size = 8, vjust = 0) +
  scale_fill_manual(values = mediation_colors) +
  labs(title = "Mediation of PRS effect on BMI by medication",
       x = "Medication type",
       y = "Effect on BMI (s.d.)",
       fill = "") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

print(prs_med_plot)

ggsave('Plots.v2/PsychiatricMeds/random_effects/mediation_analysis_PRS.png', prs_med_plot, height = 6, width = 10)









