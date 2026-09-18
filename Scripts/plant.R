# load packages
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

# import aphid data to look at maximum populations
aphid_df <- read_csv("Data/aphid.csv")
aphid_df <- aphid_df[!is.na(aphid_df$focal_wingless),]
aphid_df <- aphid_df[!is.na(aphid_df$focal_winged),]

# find maximum number of aphids per cage and use this as a 
#covariate in plant analysis; copy this column over to plant
aphid_df$total_aphids <- aphid_df$focal_winged + aphid_df$focal_wingless

aphid_max <- aphid_df %>% 
  group_by(cage) %>% 
  slice_max(total_aphids)

aphid_max <- aphid_max %>% 
  select(total_aphids, cage)

df <- left_join(df, aphid_max, by = "cage")

# change variables to appropriate category and drop NAs
df <- df %>% mutate(height = as.numeric(height))
df <- df %>% mutate(leaves = as.numeric(leaves))
df <- df[!is.na(df$height),]
df <- df[!is.na(df$leaves),]
df$period <- factor(df$period, levels=c("start", "end"))
df$warmingtreatment <- factor(df$warmingtreatment, levels=c("unwarmed", "warmed"))
df$predatortreatment <- factor(df$predatortreatment, levels=c("nopredator", "predator"))

#get means
df_means <- df %>%
  group_by(period, predatortreatment, warmingtreatment) %>%
  summarize(mean_leaves = mean(leaves, na.rm = TRUE), 
            n=n(), 
            sd = sd(leaves), 
            se = sd/sqrt(n),
            mean_height = mean(height, na.rm = TRUE),
            sd_height = sd(height),
            se_height = sd_height/sqrt(n))


#find threshold value for splitting dataset (looking at median # of aphids) and save as variable
threshold <- median(df$total_aphids, na.rm = TRUE)

df <- df %>%
  mutate(aphid_group = ifelse(total_aphids <= threshold, "Low aphid abundance", "High aphid abundance"))

df_means <- df %>%
  group_by(period, warmingtreatment, predatortreatment, aphid_group) %>%
  summarise(mean_height = mean(height, na.rm = TRUE), se = sd(height, na.rm = TRUE) / sqrt(n()), .groups = "drop")

custom_labels <- c("Start", "End")

#### Plant height figure (Fig. 3) #####
heightfig <- ggplot(df_means,
                    aes(x = period,
                        y = mean_height,
                        color = warmingtreatment,
                        shape = predatortreatment,
                        linetype = predatortreatment)) +
  labs(x = "Time point",
       y = "Height of focal plant (cm)") +
  theme_tess() +
  scale_color_manual(name = "Warming",
                     labels = c("No", "Yes"),
                     values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator",
                     labels = c("No", "Yes"),
                     values = c(1, 16)) +
  scale_linetype_manual(name = "Predator",
                        labels = c("No", "Yes"),
                        values = c("dashed", "solid")) +
  geom_point(aes(group = interaction(warmingtreatment, predatortreatment)),
             position = position_dodge(width = 0.5),
             size = 2) +
  geom_errorbar(aes(ymin = mean_height - se,
                    ymax = mean_height + se,
                    group = interaction(warmingtreatment, predatortreatment)),
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(data = df,
             aes(x = period,
                 y = height,
                 group = cage),
             position = position_dodge(width = 0.5),
             alpha = 1/10) +
  geom_line(data = df,
            aes(x = period,
                y = height,
                group = cage,
                linetype = predatortreatment),
            alpha = 1/10,
            position = position_dodge(width = 0.5)) +
  geom_line(aes(group = interaction(warmingtreatment, predatortreatment),
                linetype = predatortreatment),
            position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  scale_x_discrete(labels = custom_labels) +
  facet_wrap(~factor(aphid_group, levels = c("Low aphid abundance", "High aphid abundance")))+
  theme(strip.text = element_text(size = 20)) +
  coord_cartesian(ylim = c(65, 110))

ggsave(file="Figures/Fig 3.pdf", heightfig, width = 22, 
       height = 16, units = "cm")

#### Plant height analysis ####

# run 3 way repeated measures anova for plant height including max aphids reached in cage. Note that "total_aphids" = max aphids, see code above for getting this from the aphid dataset
# scale aphids first because huge numbers
df$aphids_sc <- scale(df$total_aphids)

model <- lmer(height ~ warmingtreatment * predatortreatment * 
                period * aphids_sc + (1 | cage), data = df)
Anova(model)
# interaction between warming*time*aphids, significant effect of max aphids
fixef(model)
