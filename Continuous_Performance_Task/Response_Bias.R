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

# Using the previously created variables (HR and FAR) to create the Response Bias variable (CB)
newdf = newdf %>%
  mutate(CB2 = -0.5*(qnorm(HR2) + qnorm(FAR2)),
         CB1 = -0.5*(qnorm(HR1) + qnorm(FAR1)),
         CB0.5 = -0.5*(qnorm(HR0.5) + qnorm(FAR0.5)),
         CB0.2 = -0.5*(qnorm(HR0.2) + qnorm(FAR0.2)))

AVGCbias = group_by(newdf, AnimalID, Genotype, Sex) %>%
  summarise("2" = mean(CB2),
            "1" = mean(CB1),
            "0.5" = mean(CB0.5),
            "0.2" = mean(CB0.2))

# Summarising (mean) performance variable, graphing and analyses
AVGCbias.melted = melt(AVGCbias, id.vars = c("AnimalID", "Genotype", "Sex"))
colnames(AVGCbias.melted)[4] = "Duration"
colnames(AVGCbias.melted)[5] = "Cbias"

summary = summarySE(AVGCbias.melted, measurevar = "Cbias", groupvars = c("Sex", "Genotype", "Duration"), na.rm = T)
summary$Genotype = factor(summary$Genotype, levels = c("WT", "CSF1R+/-"))

# Breaking down into different sexes for graphing
females = summary %>%
  filter(Sex == "Female")

males = summary %>%
  filter(Sex == "Male")

# Plotting female Response Bias data and exporting plot as .tiff
summaryf = females %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g1 = ggplot(summaryf, aes(x = Duration, y = Cbias, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = Cbias-se, ymax = Cbias+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(-0.15, 0.8) +
  ylab("Response Bias (c)") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29),
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Females")

ggsave(filename = "cbiasPROBE1runfemales.tiff", plot = g1)

# Plotting male Response Bias data and exporting plot as .tiff
summarym = males %>%
  mutate(s_g = ifelse(Duration == "z", NA, paste(Sex, Genotype, sep = ",")))

g2 = ggplot(summarym, aes(x = Duration, y = Cbias, shape = Genotype, group = s_g)) + 
  geom_errorbar(aes(ymin = Cbias-se, ymax = Cbias+se), width = .2) +
  geom_line () +
  geom_point(size = 9) +
  ylim(-0.03, 0.80) +
  ylab("Response Bias (c)") +
  xlab("Stimulus Duration (s)") +
  scale_x_discrete(breaks = c("2", "1", "0.5", "0.2"), drop = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 29),
        plot.title = element_text(size = 27, face = "bold")) +
  ggtitle("Males")

ggsave(filename = "cbiasPROBE1runmales.tiff", plot = g2)

# Running stats (repeated measures ANOVA + Sidak and Benjamini-Hochberg procedures)
# sex segregated
# Females
statsfemales = AVGCbias.melted %>%
  filter(Sex == "Female")

result.anova = aov(formula = Cbias ~ Genotype * Duration + Error(AnimalID),
                   data = statsfemales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")

# Males
statsmales = AVGCbias.melted %>%
  filter(Sex == "Male")

result.anova = aov(formula = Cbias ~ Genotype * Duration + Error(AnimalID),
                   data = statsmales)
summary(result.anova)
emm = emmeans(result.anova, ~ Genotype | Duration)
contrast(emm, interaction = TRUE, adjust = "Sidak")
pairs(emm, adjust = "BH")