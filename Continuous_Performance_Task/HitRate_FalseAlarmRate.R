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

# Importing data and renaming variables
file = file.choose()
probewholesession = read.csv(file, stringsAsFactors = FALSE)

probewholesession$Genotype = ifelse(probewholesession$Genotype == "C57BL/6J", "WT", "CSF1R+/-")
probewholesession$Sex = ifelse(probewholesession$Sex == "f", "Female", "Male")
probewholesession$Age = ifelse(probewholesession$Age == "3_6       ", "3.6", "z")

# Filtering for a specific timepoint
firstcpt = probewholesession %>%
  filter(Age == "3.6")

# Creating Hit Rate and False Alarm Rate variables
newdf = firstcpt %>%
  mutate(HR2 = ifelse(AVG_End.Summary...Miss.at.2s.SD == 0, 
                      1 - 1/(2*AVG_End.Summary...Hits.at.2s.SD), 
                      AVG_End.Summary...Hits.at.2s.SD/(AVG_End.Summary...Hits.at.2s.SD + 
                                                         AVG_End.Summary...Miss.at.2s.SD)),
         HR1 = ifelse(AVG_End.Summary...Miss.at.1s.SD == 0, 
                      1 - 1/(2*AVG_End.Summary...Hits.at.1s.SD), 
                      AVG_End.Summary...Hits.at.1s.SD/(AVG_End.Summary...Hits.at.1s.SD + 
                                                         AVG_End.Summary...Miss.at.1s.SD)),
         HR0.5 = ifelse(AVG_End.Summary...Miss.at.0.5s.SD == 0, 
                        1 - 1/(2*AVG_End.Summary...Hits.at.0.5s.SD), 
                        AVG_End.Summary...Hits.at.0.5s.SD/(AVG_End.Summary...Hits.at.0.5s.SD + 
                                                             AVG_End.Summary...Miss.at.0.5s.SD)),
         HR0.2 = ifelse(AVG_End.Summary...Miss.at.0.2s.SD == 0, 
                        1 - 1/(2*AVG_End.Summary...Hits.at.0.2s.SD), 
                        AVG_End.Summary...Hits.at.0.2s.SD/(AVG_End.Summary...Hits.at.0.2s.SD + 
                                                             AVG_End.Summary...Miss.at.0.2s.SD)))

newdf = newdf %>%
  mutate(FAR2 = ifelse(AVG_End.Summary...Mistake.at.2s.SD == 0,
                       1/(2*AVG_End.Summary...CorrectRejection.at.2s.SD),
                       AVG_End.Summary...Mistake.at.2s.SD/(AVG_End.Summary...Mistake.at.2s.SD + 
                                                             AVG_End.Summary...CorrectRejection.at.2s.SD)),
         FAR1 = ifelse(AVG_End.Summary...Mistake.at.1s.SD == 0,
                       1/(2*AVG_End.Summary...CorrectRejection.at.1s.SD),
                       AVG_End.Summary...Mistake.at.1s.SD/(AVG_End.Summary...Mistake.at.1s.SD + 
                                                             AVG_End.Summary...CorrectRejection.at.1s.SD)),
         FAR0.5 = ifelse(AVG_End.Summary...Mistake.at.0.5s.SD == 0,
                         1/(2*AVG_End.Summary...CorrectRejection.at.0.5s.SD),
                         AVG_End.Summary...Mistake.at.0.5s.SD/(AVG_End.Summary...Mistake.at.0.5s.SD + 
                                                                 AVG_End.Summary...CorrectRejection.at.0.5s.SD)),
         FAR0.2 = ifelse(AVG_End.Summary...Mistake.at.0.2s.SD == 0,
                         1/(2*AVG_End.Summary...CorrectRejection.at.0.2s.SD),
                         AVG_End.Summary...Mistake.at.0.2s.SD/(AVG_End.Summary...Mistake.at.0.2s.SD + 
                                                                 AVG_End.Summary...CorrectRejection.at.0.2s.SD)))

# Summarising (mean) performance variables, graphing and analyses
# Hit Rate
AVGHitRate = group_by(newdf, AnimalID, Genotype, Sex) %>%
  summarise("2" = mean(HR2),
            "1" = mean(HR1),
            "0.5" = mean(HR0.5),
            "0.2" = mean(HR0.2))

AVGHitRate.melted = melt(AVGHitRate, id.vars = c("AnimalID", "Genotype", "Sex"))
colnames(AVGHitRate.melted)[4] = "Duration"
colnames(AVGHitRate.melted)[5] = "HitRate"

summary = summarySE(AVGHitRate.melted, measurevar = "HitRate", groupvars = c("Sex", "Genotype", "Duration"))
summary$Genotype = factor(summary$Genotype, levels = c("WT", "CSF1R+/-"))

# Breaking down into different sexes for graphing
females = summary %>%
  filter(Sex == "Female")
males = summary %>%
  filter(Sex == "Male")

# Plotting female hit rate data and exporting plot as .tiff
summaryf = females %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g1 = ggplot(summaryf, aes(x = Duration, y = HitRate, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = HitRate-se, ymax = HitRate+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(0.25, 1) +
  ylab("Hit Rate") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29), 
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Females")

ggsave(filename = "hitratePROBE1runfemales.tiff", plot = g1)

# Plotting male hit rate data and exporting plot as .tiff
summarym = males %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g2 = ggplot(summarym, aes(x = Duration, y = HitRate, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = HitRate-se, ymax = HitRate+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(0.25, 1) +
  ylab("Hit Rate") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29),
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Males")

ggsave(filename = "hitratePROBE1runmales.tiff", plot = g2)

# Running stats (repeated measures ANOVA + Sidak and Benjamini-Hochberg procedures)
# sex segregated
# Females
statsfemales = AVGHitRate.melted %>%
  filter(Sex == "Female")

result.anova = aov(formula = HitRate ~ Genotype * Duration + Error(AnimalID),
                   data = statsfemales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")

# Males
statsmales = AVGHitRate.melted %>%
  filter(Sex == "Male")

result.anova = aov(formula = HitRate ~ Genotype * Duration + Error(AnimalID),
                   data = statsmales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")

# Summarising (mean) performance variables, graphing and analyses
# False Alarm Rate
AVGFalseAlarmRate = group_by(newdf, AnimalID, Genotype, Sex) %>%
  summarise("2" = mean(AVG_False.Alarm.Rate.at.2s),
            "1" = mean(AVG_False.Alarm.Rate.at.1s),
            "0.5" = mean(AVG_False.Alarm.Rate.at.0.5s),
            "0.2" = mean(AVG_False.Alarm.Rate.at.0.2s))

AVGFalseAlarmRate.melted = melt(AVGFalseAlarmRate, id.vars = c("AnimalID", "Genotype", "Sex"))
colnames(AVGFalseAlarmRate.melted)[4] = "Duration"
colnames(AVGFalseAlarmRate.melted)[5] = "FalseAlarmRate"

summary = summarySE(AVGFalseAlarmRate.melted, measurevar = "FalseAlarmRate", groupvars = c("Sex", "Genotype", "Duration"))
summary$Genotype = factor(summary$Genotype, levels = c("WT", "CSF1R+/-"))

# Breaking down into different sexes for graphing
females = summary %>%
  filter(Sex == "Female")
males = summary %>%
  filter(Sex == "Male")

# Plotting female false alarm rate data and exporting plot as .tiff
summaryf = females %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g3 = ggplot(summaryf, aes(x = Duration, y = FalseAlarmRate, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = FalseAlarmRate-se, ymax = FalseAlarmRate+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(0, 0.2) +
  ylab("False Alarm Rate") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29),
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Females")

ggsave(filename = "falsealarmratePROBE1runfemales.tiff", plot = g3)

# Plotting male false alarm rate data and exporting plot as .tiff
summarym = males %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g4 = ggplot(summarym, aes(x = Duration, y = FalseAlarmRate, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = FalseAlarmRate-se, ymax = FalseAlarmRate+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(0, 0.25) +
  ylab("False Alarm Rate") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29),
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Males")

ggsave(filename = "falsealarmratePROBE1runmales.tiff", plot = g4)

# Running stats (repeated measures ANOVA + Sidak and Benjamini-Hochberg procedures)
# sex segregated
# Females
statsfemales = AVGFalseAlarmRate.melted %>%
  filter(Sex == "Female")

result.anova = aov(formula = FalseAlarmRate ~ Genotype * Duration + Error(AnimalID),
                   data = statsfemales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")

# Males
statsmales = AVGFalseAlarmRate.melted %>%
  filter(Sex == "Male")

result.anova = aov(formula = FalseAlarmRate ~ Genotype * Duration + Error(AnimalID),
                   data = statsmales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")