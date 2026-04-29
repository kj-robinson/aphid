#### Plant data analysis ####

library(readxl)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(car)
library(ggpubr)
library(rstatix)
library(lme4)
library(nlme)

# set up theme
theme_tess <- function () {
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.text.x=element_text(size=20))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.x=element_text(size=20))+
    theme(axis.title.y=element_text(size=20))+
    theme(plot.title = element_text(hjust = 0.5,size=30))+
    theme(legend.title = element_text(size = 18))+
    theme(legend.title=element_text(size=16,face="bold"))+
    theme(legend.text=element_text(size = 16))
}
theme_set(theme_tess())

# import data
df <- read_csv("Data/plant.csv") 
view(df)

# change variables to appropriate category and drop NAs
df <- df %>% mutate(height = as.numeric(height))
df <- df %>% mutate(leaves = as.numeric(leaves))
df <- df[!is.na(df$height),]
df <- df[!is.na(df$leaves),]
df$period <- factor(df$period, levels=c("start", "end"))
df$warmingtreatment <- factor(df$warmingtreatment, levels=c("unwarmed", "warmed"))
df$predatortreatment <- factor(df$predatortreatment, levels=c("nopredator", "predator"))

df_means <- df %>%
  group_by(period, predatortreatment, warmingtreatment) %>%
  summarize(mean_leaves = mean(leaves, na.rm = TRUE), 
            n=n(), 
            sd = sd(leaves), 
            se = sd/sqrt(n)) 

# create figure showing number of leaves on plant

custom_labels <- c("Start", "End")

leavesfig <- ggplot(df_means,
                    aes(x = period,
                        y = mean_leaves,
                        color = warmingtreatment,
                        shape = predatortreatment)) +
  labs(x = "Period",
       y = "Number of living leaves on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming",
                     labels = c("No", "Yes"),
                     values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", 
                     labels = c("No", "Yes"),
                     values = c(1, 16)) +
  geom_point(aes(group = interaction(warmingtreatment, predatortreatment)),
             position = position_dodge(width = 0.5),
             size = 2) +
  geom_errorbar(aes(ymin = mean_leaves - se,
                    ymax = mean_leaves + se,
                    group = interaction(warmingtreatment, predatortreatment)),
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(data = df,
             aes(x = period,
                 y = leaves,
                 group = cage),
             position = position_dodge(width = 0.5),
             alpha = 1/10) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_line(data = df,
            aes(x = period,
                y = leaves,
                group = cage,
                linetype = predatortreatment),
            alpha = 1/10,
            position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator",
                        labels = c("No", "Yes"),
                        values = c("dashed", "solid")) +
  scale_x_discrete(labels = custom_labels) +
  geom_line(data = df_means,
            aes(x = period,
                y = mean_leaves,
                linetype = predatortreatment,
                group = interaction(warmingtreatment, predatortreatment)),
            position = position_dodge(width = 0.5))

leavesfig

# summary stats for leaves
df %>%
  group_by(warmingtreatment, predatortreatment, period) %>% 
  summarize(mean(leaves))

# look at general outliers
df %>%
  group_by(period) %>%
  identify_outliers(leaves)

# run 3 way repeated measures anova
model <- lmer(leaves ~ warmingtreatment * predatortreatment * period + (1 | cage), data = df)

Anova(model)
