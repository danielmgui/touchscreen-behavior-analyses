# Loading Libraries
library(plyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(stringr)
library(rcompanion)
library(rstatix)
library(Rmisc)
library(emmeans)

# Importing data 
raw.file = file.choose()
data = read.csv(raw.file)

# Summarising (mean) performance variables, graphing and analyses
avgAcc = group_by(data, ID, Genotype, Sex) %>%
  summarise(S3 = mean(S3.Overall),
            S1 = mean(S1.Overall))

avgAcc.melted = melt(avgAcc, id.vars = c("ID", "Genotype", "Sex"))
colnames(avgAcc.melted)[4] = "Separation.Level"
colnames(avgAcc.melted)[5] = "Accuracy"

summary = summarySE(avgAcc.melted, measurevar = "Accuracy", groupvars = c("Sex", "Genotype", "Separation.Level"))
summary$Genotype = factor(summary$Genotype, levels = c("WT", "CSF1R+/-"))

# Breaking down into different sexes for graphing
females = summary %>%
  filter(Sex == "Female")
males = summary %>%
  filter(Sex == "Male")

# Plotting female accuracy data and exporting plot as .tiff
g1 = ggplot(females, aes(x = Separation.Level, y = Accuracy, fill = Genotype)) + 
  geom_errorbar(aes(ymin = Accuracy-se, ymax = Accuracy+se), position=position_dodge(.9), width = .2) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 50, linetype = "dashed", size = 1.25) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  ylab("Accuracy (%)") +
  xlab("Separation Level") +
  theme_classic() +
  theme(text = element_text(size = 21), 
        plot.title = element_text(size = 19, face = "bold")) +
  scale_fill_grey() +
  ggtitle("Females")

ggsave(filename = "femalesAccuracySEPARATIONPROBE.tiff", plot = g1)

# Plotting male accuracy data and exporting plot as .tiff
g2 = ggplot(males, aes(x = Separation.Level, y = Accuracy, fill = Genotype)) + 
  geom_errorbar(aes(ymin = Accuracy-se, ymax = Accuracy+se), position=position_dodge(.9), width = .2) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 50, linetype = "dashed", size = 1.25) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  ylab("Accuracy (%)") +
  xlab("Separation Level") +
  theme_classic() +
  theme(text = element_text(size = 21), 
        plot.title = element_text(size = 19, face = "bold")) +
  scale_fill_grey() +
  ggtitle("Males")

ggsave(filename = "malesAccuracySEPARATIONPROBE.tiff", plot = g2)

# Running stats (repeated measures ANOVA + Sidak and Benjamini-Hochberg procedures)
# sex segregated
# Females
statsfemales = avgAcc.melted %>%
  filter(Sex == "Female")

result.anova = aov(formula = Accuracy ~ Genotype * Separation.Level + Error(ID),
                   data = statsfemales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Separation.Level)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")

# Males
statsmales = avgAcc.melted %>%
  filter(Sex == "Male")

result.anova = aov(formula = Accuracy ~ Genotype * Separation.Level + Error(ID),
                   data = statsmales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Separation.Level)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")