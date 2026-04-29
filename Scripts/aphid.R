#### DATAFRAME SETUP + LIBRARY ####

# check working directory
getwd()

# load packages
library(dplyr)
library(ggplot2)
library(cowplot)
library(lubridate)
library(usethis)
library(car)
library(ggpubr)
library(rstatix)
library(lme4)
library(stats)
library(scales)
library(DHARMa)
library(glmmTMB)
library(tidyverse)
library(gridExtra)
library(gtable)
library(grid)

# open data file
countdata <- read.csv("./Data/aphid.csv")
countdata

# change all count data into numeric variable
countdata <- countdata %>% mutate(sticky_winged = as.numeric(sticky_winged))
countdata <- countdata %>% mutate(sticky_wingless = as.numeric(sticky_wingless))
countdata <- countdata %>% mutate(sentinel_winged = as.numeric(sentinel_winged))
countdata <- countdata %>% mutate(sentinel_wingless = as.numeric(sentinel_wingless))
countdata <- countdata %>% mutate(focal_winged = as.numeric(focal_winged))
countdata <- countdata %>% mutate(focal_wingless = as.numeric(focal_wingless))
countdata <- countdata %>% mutate(date = as.Date(date))

# predatortreatment and warmingtreatment as factors
countdata <- countdata %>% mutate(predatortreatment = as.factor(predatortreatment))
countdata <- countdata %>% mutate(warmingtreatment = as.factor(warmingtreatment))

countdata <- countdata[!is.na(countdata$focal_wingless),]
countdata <- countdata[!is.na(countdata$focal_winged),]

# remove excluded cages (WP4, WN3, WN9, UP2, UP7)
countdata <- countdata %>%
  filter(!cage %in% c("WP4", "WN3", "WN9", "UP2", "UP7"))

# add column converting date to julian date
countdata$jd <- yday(countdata$date)
  
# set up theme
theme_tess <- function () {
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)))+
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

# make a totalled column of aphids found on focal plant
countdata <- countdata %>%
  mutate(focal_aphids = focal_wingless + focal_winged, 
         focal_proportion = focal_winged/focal_aphids)

View(countdata)

# summarize counts by means, grouping by date observed, predatortreatment and warmingtreatment
countdata_means <- countdata %>%
  group_by(date, predatortreatment, warmingtreatment) %>%
  summarize(mean_aphids = mean(focal_aphids, na.rm = TRUE), 
            n=n(),
            sd = sd(focal_aphids), 
            se = sd/sqrt(n))

countdata_means$date <- as.Date(countdata_means$date)
countdata$date <- as.Date(countdata$date)

#### TEMPERATURE PLOTS ####
# load data
tempdata <- read.csv("./Data/temperature.csv")
View(tempdata)
tempdata <- tempdata[!is.na(tempdata$temperature),]

# separate data and time into two columns
df <- tempdata %>%
  mutate(datetime = str_squish(datetime)) %>% #remove extra white spaces
  separate(datetime, into = c("date", "time"), sep = " ") 
view(df)

# make a new column that specifies if a a row is in the day or night - day is 7:00 to 19:00
df <- df %>%
  mutate(
    time_parsed = hms(time),
    hour = hour(time_parsed),
    minute = minute(time_parsed),
    period = case_when(
      hour > 7 & hour < 19 ~ "day",
      hour == 7 & minute > 0 ~ "day",
      hour == 19 & minute == 0 ~ "day",
      TRUE ~ "night"))

# make a column for date that goes 1 --> last date
df <- df %>%
  mutate(date = mdy(date),                # convert string to Date format (month-day-year)
         day_number = dense_rank(date))

# find average daytime temperatures for each cage per day
df_dt <- df %>% 
  group_by(cage, day_number, date, warmingtreatment) %>% 
  summarise(meantemp = mean(temperature),
            n=n(),
            sd = sd(temperature),
            setemp = sd / sqrt(n))

# find average daytime temperatures for warmed and unwarmed on each day of experiment
df_dt_means <- df_dt %>% 
  group_by(day_number, date, warmingtreatment) %>% 
  summarise(meantempdaily = mean(meantemp))

# find maximum daytime temperatures for each cage per day
df_dt_max <- df %>% 
  group_by(cage, day_number, date, warmingtreatment) %>% 
  summarise(cagemax = max(temperature))

# find average maximum daytime temperatures for warmed and unwarmed on each day of experiment
df_dt_maxavg <- df_dt_max %>% 
  group_by(day_number, date, warmingtreatment) %>% 
  summarise(daymax = max(cagemax))

# carries over last sampled number of lady beetles for days that were not sampled
countdata_filled <- countdata %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day")
  ) %>%
  fill(ladybug_total, .direction = "up")

dailymeansdn <- df %>%
  group_by(day_number, cage, warmingtreatment) %>% #put warmingtreatment here so it keeps that column
  summarise(meantemp = mean(temperature),
            n=n(),
            sd = sd(temperature),
            setemp = sd / sqrt(n))

treatmentmeansdn <- dailymeansdn %>%
  group_by(warmingtreatment) %>%
  summarise(mean = mean(meantemp),
            n=n(),
            sd = sd(meantemp),
            se = sd / sqrt(n))

# creating a data point that just tells my plot where to put the mean points
summary_date <- as.Date("2025-08-31")

# jittering means
treatmentmeansdn <- treatmentmeansdn |>
  mutate(
    x_plot = case_when(
      warmingtreatment == "unwarmed"  ~ summary_date - 2,
      warmingtreatment == "warmed" ~ summary_date + 2))

# plot temp means throughout growing season
avgtemp <- ggplot(df_dt_means, aes(x = date, y = meantempdaily)) +
  geom_point(
    aes(colour = warmingtreatment),
    size = 2) +
  geom_line(
    aes(
      colour = warmingtreatment,
      linetype = warmingtreatment,
      group = warmingtreatment),
    alpha = 0.5) +
  geom_step(
    data = countdata_filled,
    aes(
      x = date,
      y = (ladybug_total - 2) * 0.4 + 12),
    inherit.aes = FALSE,
    colour = "forestgreen",
    linewidth = 1,
    na.rm = TRUE) +
  labs(x = "Date") +
  theme_tess() +
  scale_color_manual(
    name = "Warming",
    labels = c("No", "Yes"),
    values = c("blue", "red")) +
  scale_linetype_manual(
    name = "Warming",
    labels = c("No", "Yes"),
    values = c("solid", "solid")) +
  scale_x_date(
    breaks = seq(
      from = as.Date("2025-07-15"),
      to   = as.Date("2025-09-02"),
      by   = "1 week"),
    date_labels = "%b %d",
    limits = ymd("2025-07-15", "2025-09-02")) +
  scale_y_continuous(
    name = "Mean daily temperature (°C)",
    sec.axis = sec_axis(
      ~ (. - 12) / 0.4 + 2,
      name = "Number of lady beetles")) +
  theme(axis.line.y.right = element_line(colour = "forestgreen"), 
         axis.ticks.y.right = element_line(colour = "forestgreen"),
         axis.text.y.right = element_text(colour = "forestgreen")) +
  geom_pointrange(
    data = treatmentmeansdn,
    aes(
      x = x_plot,
      y = mean,
      ymin = mean - se,
      ymax = mean + se,
      colour = warmingtreatment),
    size = 0.75,
    linewidth = 1,
    shape = 17,
    show.legend = FALSE)
  
avgtemp

# plot treatment maximums thoughout growing season
maxtemp <- ggplot(df_dt_maxavg, aes(x = date, y = daymax)) +
  geom_point(
    aes(colour = warmingtreatment),
    size = 2) +
  geom_line(
    aes(
      colour = warmingtreatment,
      linetype = warmingtreatment,
      group = warmingtreatment),
    alpha = 0.5) +
  labs(x = "Date", y = "Mean maximum daily temperature (°C)") +
  theme_tess() +
  scale_color_manual(
    name = "Warming",
    labels = c("No", "Yes"),
    values = c("blue", "red")) +
  scale_linetype_manual(
    name = "Warming",
    labels = c("No", "Yes"),
    values = c("solid", "solid")) +
  scale_x_date(
    breaks = seq(
      from = as.Date("2025-07-15"),
      to   = as.Date("2025-09-02"),
      by   = "1 week"),
    date_labels = "%b %d",
    limits = ymd("2025-07-15",
                 "2025-09-02"))

maxtemp

#### GLOBAL APHID ABUNDANCE OVER TIME ####

# make plot showing number of aphids on focal plant over time
focalplantcount <- ggplot(countdata_means, aes(x = date, y = mean_aphids, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Number of aphids on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = mean_aphids - se, ymax = mean_aphids + se), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata, aes(x = date, y = focal_aphids, group = cage, linetype = warmingtreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "1 week"),
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_means, aes(x = date, y = mean_aphids, linetype = predatortreatment))
  
focalplantcount
windows();focalplantcount

## multipanelled temp/abundance plot ##
grid.arrange(focalplantcount, avgtemp, nrow = 2)
grid.arrange(avgtemp, maxtemp, nrow = 2)

g2 <- ggplotGrob(focalplantcount)
g3 <- ggplotGrob(avgtemp)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

# filter to early count dates to look at predator effects
countdata_zoom <- countdata %>%
  filter(date >= as.Date("2025-07-15") & date <= as.Date("2025-07-21"))

countdata_means_zoom <- countdata_means %>%
  filter(date >= as.Date("2025-07-15") & date <= as.Date("2025-07-21"))

zoomedcount <- ggplot(countdata_means_zoom, aes(x = date, y = mean_aphids, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Number of aphids on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = mean_aphids - se, ymax = mean_aphids + se), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_zoom, aes(x = date, y = focal_aphids, group = cage, linetype = warmingtreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "3 days"),
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_means_zoom, aes(x = date, y = mean_aphids, linetype = predatortreatment))
zoomedcount

windows();focalplantcount

# ensure jd is numeric and scaled for ease of analysis
countdata$jd <- as.numeric(countdata$jd)
countdata$jd_sc <- scale(countdata$jd)

# model
lmfullmodel<-lmer(focal_aphids ~ predatortreatment * warmingtreatment * jd_sc + predatortreatment * warmingtreatment * I(jd_sc^2) + (1 | cage), data = countdata)
Anova(lmfullmodel, type=2)

summary(lmfullmodel)

#Response: focal_aphids
#Chisq Df Pr(>Chisq)    
#predatortreatment              0.0537  1    0.81681    
#warmingtreatment                3.2335  1    0.07215 .  
#jd                   116.3792  1    < 2e-16 ***
#predatortreatment:warmingtreatment      0.0053  1    0.94202    
#predatortreatment:jd           3.5685  1    0.05889 .  
#warmingtreatment:jd             8.5685  1    0.00342 ** 
#predatortreatment:warmingtreatment:jd   0.1533  1    0.69542    


#significant interaction between time and warmingtreatment
#significant effect of date
#significant effect of warmingtreatment


#### DAY BY DAY PRED/WARM ANALYSIS ####

# grouping by date for individual day-by-day ANOVAs
jd196 <- countdata %>% filter(jd == 196)
jd196

jd198 <- countdata %>% filter(jd == 198)
jd198

jd202 <- countdata %>% filter(jd == 202)
jd202

jd205 <- countdata %>% filter(jd == 205)
jd205

jd209 <- countdata %>% filter(jd == 209)
jd209

jd212 <- countdata %>% filter(jd == 212)
jd212

jd217 <- countdata %>% filter(jd == 217)
jd217

jd219 <- countdata %>% filter(jd == 219)
jd219

jd223 <- countdata %>% filter(jd == 223)
jd223

jd226 <- countdata %>% filter(jd == 226)
jd226

jd230 <- countdata %>% filter(jd == 230)
jd230

jd233 <- countdata %>% filter(jd == 233)
jd233

jd237 <- countdata %>% filter(jd == 237)
jd237

jd240 <- countdata %>% filter(jd == 240)
jd240

jd245 <- countdata %>% filter(jd == 245)
jd245

###day by day analysis###

countdata$predatortreatment <- factor(countdata$predatortreatment)
countdata$warmingtreatment <- factor(countdata$warmingtreatment)

#July 17
lm198 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd198)
Anova(lm198, type=2)
#significant effect of warmingtreatment and an interaction with warmingxpred

#July 21
lm202 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd202)
Anova(lm202, type=2)
#nothing significant

#July 24
lm205 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd205)
Anova(lm205, type=2)
#nothing significant

#July 28
lm209 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd209)
Anova(lm209, type=2)
#predatortreatment weakly significant

# July 31
lm212 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd212)
Anova(lm212, type=2)
#predatortreatment marginally significant

#August 5
lm217 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd217)
Anova(lm217, type=2)
#predatortreatment marginally significant

#August 7
lm219 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd219)
Anova(lm219, type=2)
#nothing significant

#August 11
lm223 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd223)
Anova(lm223, type=2)
#nothing significant

#August 14
lm226 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd226)
Anova(lm226, type=2)
#nothing significant

#August 18
lm230 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd230)
Anova(lm230, type=2)
#strong significant effect of warmingtreatment

#August 21
lm233 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd233)
Anova(lm233, type = 2)
#significant effect of warmingtreatment

#August 25
lm237 <- lm(focal_aphids ~ predatortreatment * warmingtreatment, data = jd237)
Anova(lm237, type = 2)
#strong significant effect of warmingtreatment

#August 28
lm240 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd240)
Anova(lm240, type = 2)
#nothing significant

#Sept 2
lm245 <- lm(focal_aphids ~ predatortreatment * warmingtreatment,data = jd245)
Anova(lm245, type = 2)
#weak significant effect of predatortreatment


#### ABUNDANCE PEAKS ####

# find peak populations per treatment
countdata_treatmentmax <- countdata %>%
  group_by(predatortreatment, warmingtreatment) %>% 
  slice_max(focal_winged + focal_wingless)

View(countdata_treatmentmax)

# find peak populations per cage to see when aphids begin to decline
countdata_cagemax <- countdata %>% 
  group_by(cage) %>% 
  slice_max(focal_winged + focal_wingless)

View(countdata_cagemax)

#### GLOBAL PROPORTION WINGED ####

# find means to plot
countdata_means <- countdata %>%
  mutate(prop = (focal_winged/focal_aphids))%>%
  group_by(date, predatortreatment, warmingtreatment) %>%
  summarize(prop_winged = mean(prop, na.rm = TRUE), 
            n=n(),
            sd_prop = sd(prop), 
            se_prop = sd_prop/sqrt(n))

countdata_means <- countdata %>%
  group_by(date, predatortreatment, warmingtreatment) %>%
  summarize(prop_winged = mean(focal_proportion, na.rm = TRUE), 
            n=sum(!is.na(focal_proportion)),
            sd_prop = sd(focal_proportion,na.rm = TRUE), 
            se_prop = sd_prop/sqrt(n))

summary(countdata_means)

# plot winged/total
propwinged <- ggplot(countdata_means, aes(x = date, y = prop_winged, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Proportion of winged aphids on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = prop_winged - se_prop, ymax = prop_winged + se_prop), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata, aes(x = date, y = (focal_winged/focal_aphids), group = cage, linetype = warmingtreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "1 week"),
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_means, aes(x = date, y = prop_winged, linetype = predatortreatment))

propwinged

windows();propwinged

## analysis winged

countdata$cage <- as.factor(countdata$cage)

dispersal_model_poly <- glmmTMB(
  cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment * jd_sc + warmingtreatment * predatortreatment * I(jd_sc^2) + (1 | cage),
  family = betabinomial(),
  data = countdata)
Anova(dispersal_model_poly, type = 2)

# checking dispersion
testDispersion(dispersal_model_poly)


#### DAY BY DAY PROPORTION WINGED ANALYSIS ####

## analysis winged by julian day

#July 17th 
glmm198 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                   family = betabinomial(), data = jd198)
Anova(glmm198,type = 2)
windows();boxplot(focal_proportion~warmingtreatment, data=jd198)
# weak effect of predatortreatment, strong effect of warmingtreatment (more winged in warmed)
# checking for overdispersion
testDispersion(glmm198)
# not overdispersed

#July 21st 
glmm202 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
               family = betabinomial(),  data = jd202)
Anova(glmm202,type = 2)
# nothing significant
# checking for overdispersion
testDispersion(glmm202)
# not overdispersed

#July 24th 
glmm205 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd205)
Anova(glmm205,type = 2)
# nothing significant
# checking for overdispersion
testDispersion(glmm205)
# overdispersed

#July 28th
glmm209 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd209)
Anova(glmm209,type = 2)
windows();boxplot(focal_proportion~warmingtreatment, data=jd209)
windows();boxplot(focal_proportion~predatortreatment, data=jd209)
# strong effect of warmingtreatment (higher in warmed), effect of predatortreatment (lower in pred), no interaction
# checking for overdispersion
testDispersion(glmm209)
# overdispersed

#July 31st
glmm212 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd212)
Anova(glmm212,type = 2)
windows();boxplot(focal_proportion~warmingtreatment, data=jd209)
windows();boxplot(focal_proportion~predatortreatment, data=jd209)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm212)
# overdispersed

#August 5th 
glmm217 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd217)
Anova(glmm217,type = 2)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm217)
# overdispersed

#August 7th 
glmm219 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd219)
Anova(glmm219,type = 2)
windows();boxplot(focal_proportion~warmingtreatment*predatortreatment, data=jd219)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm219)
# overdispersed

#August 11
glmm223 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd223)
Anova(glmm223, type = 2)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm223)
# overdispersed

#August 14
glmm226 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd226)
Anova(glmm226,type = 2)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm226)
# overdispersed

# August 18
glmm230 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd230)
Anova(glmm230,type = 2)
# strong interactin and predator effect
# checking for overdispersion
testDispersion(glmm230)
# overdispersed

# August 21
glmm233 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd233)
Anova(glmm233,type = 2)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm233)
# overdispersed

# August 25
glmm237 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd237)
Anova(glmm237,type = 2)
# weak effect pred, strong effect warmingtreatment
# checking for overdispersion
testDispersion(glmm237)
# overdispersed

# August 28
glmm240 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd240)
Anova(glmm240,type = 2)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm240)
# overdispersed

# Sept 2
glmm245 <- glmmTMB(cbind(focal_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd245)
Anova(glmm245,type = 2)
# checking for overdispersion
testDispersion(glmm245)
# overdispersed


#### GLOBAL DISPERSAL ANALYSIS ####

# prepare data frames and limit survey dates to before sentinel plant is removed from cages
countdata_dispersed <- countdata %>%
  filter(date >= as.Date("2025-07-24") & date <= as.Date("2025-08-05")) %>% 
  mutate(sentinel_aphids = sentinel_winged + sentinel_wingless,
         sticky_aphids = sticky_winged + sticky_wingless,
         dispersed = (sentinel_aphids + sticky_aphids)/(sentinel_aphids + sticky_aphids + focal_aphids),
         sentinel_dispersed = (sentinel_aphids / (sentinel_aphids + sticky_aphids)))

# stop producing na values
countdata_dispersed_means <- countdata_dispersed %>%
  group_by(date, predatortreatment, warmingtreatment) %>%
  summarize(
    n = sum(!is.na(dispersed)),
    dispersed_mean = mean(dispersed, na.rm = TRUE),
    sd_dispersed = sd(dispersed, na.rm = TRUE),
    se_dispersed = sd_dispersed / sqrt(n),
    
    sentinel_dispersed_mean = mean(sentinel_dispersed, na.rm = TRUE),
    sd_sentinel_dispersed = sd(sentinel_dispersed, na.rm = TRUE),
    se_sentinel_dispersed = sd_sentinel_dispersed / sqrt(n),
    .groups = "drop"
  )

summary(countdata_dispersed_means)

 # plot total aphids dispersed over total aphids on focal plant
dispersedplot <- ggplot(countdata_dispersed_means, aes(x = date, y = dispersed_mean, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Proportion of dispersed aphids from focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = dispersed_mean - se_dispersed, ymax = dispersed_mean + se_dispersed), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_dispersed, aes(x = date, y = dispersed, group = cage, linetype = warmingtreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_dispersed_means, aes(x = date, y = dispersed_mean, linetype = predatortreatment), position = position_dodge(width = 0.5))

dispersedplot

windows();dispersedplot

# run lmm
dispersed_model <- lmer(dispersed ~ warmingtreatment * predatortreatment * jd + (1 | cage), data = countdata_dispersed)
Anova(dispersed_model, type = 2)

testDispersion(dispersed_model)

# plot dispersal to sentinel plant (success) vs to sticky card (failure)
sentineldispersedplot <- ggplot(countdata_dispersed_means, aes(x = date, y = sentinel_dispersed_mean, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Proportion of dispersed aphids to sentinel plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("blue", "red")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = sentinel_dispersed_mean - se_sentinel_dispersed, ymax = sentinel_dispersed_mean + se_sentinel_dispersed), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_dispersed, aes(x = date, y = sentinel_dispersed, group = cage, linetype = warmingtreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_dispersed_means, aes(x = date, y = sentinel_dispersed_mean, linetype = predatortreatment), position = position_dodge(width = 0.5))

sentineldispersedplot

windows();sentineldispersedplot

## analysis dispersed to sentinel plant vs sticky card

# rescale jd and focal_aphids because R would not recognize it (because it's so different from the other predictors that are factors)
countdata_dispersed$jd_sc <- scale(countdata_dispersed$jd)
countdata_dispersed$focal_aphids_sc <- scale(countdata_dispersed$focal_aphids)

# run glmm
sentinel_dispersed_model <- glmer(cbind(sentinel_aphids, sticky_aphids) ~ warmingtreatment * predatortreatment * jd_sc + (1 | cage),
                         family = binomial, data = countdata_dispersed)
Anova(sentinel_dispersed_model, type = 2)

testDispersion(sentinel_dispersed_model)

#### DAY BY DAY DISPERSAL ANALYSIS ####

# make dataframe subset incl. dispersal data

djd205 <- countdata_dispersed[which(countdata_dispersed$jd == '205'), ]
djd205

djd209 <- countdata_dispersed[which(countdata_dispersed$jd == '209'), ]
djd209

djd212 <- countdata_dispersed[which(countdata_dispersed$jd == '212'), ]
djd212

djd217 <- countdata_dispersed[which(countdata_dispersed$jd == '217'), ]
djd217


### analysis dispersed to sentinel by julian day ###

djd205 <- djd205[(djd205$sentinel_aphids + djd205$sticky_aphids) > 0, ]
glmm205 <- glm(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment * focal_aphids_sc,
                   family = binomial, data = djd205)
Anova(glmm205, type = 2)
testDispersion((glmm205))
# not overdispersed

djd209 <- djd209[(djd209$sentinel_aphids + djd209$sticky_aphids) > 0, ]
glmm209 <- glm(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment * focal_aphids_sc, 
                   family = binomial, data = djd209)
Anova(glmm209, type = 2)
testDispersion((glmm209))
# overdispersed

djd212 <- djd212[(djd212$sentinel_aphids + djd212$sticky_aphids) > 0, ]
glmm212 <- glm(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment * focal_aphids_sc,
                   family = binomial, data = djd212)
Anova(glmm212, type = 2)
testDispersion((glmm212))
# overdispersed

djd217 <- djd217[(djd217$sentinel_aphids + djd217$sticky_aphids) > 0, ]
glmm217 <- glm(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment * focal_aphids_sc,
                   family = binomial, data = djd217)
Anova(glmm217, type = 2)
testDispersion((glmm217))
# overdispersed

### using glmmTMB package which allows for quasi/betabinomials and fixes most overdispersion
glmm205 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd205)
Anova(glmm205, type = 2)

glmm209 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment, 
                 family = betabinomial, data = djd209)
Anova(glmm209, type = 2)

glmm212 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd212)
Anova(glmm212, type = 2)

glmm217 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd217)
Anova(glmm217, type = 2)

# total dispersal
glmm205 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd205)
Anova(glmm205, type = 2)

glmm209 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd209)
Anova(glmm209, type = 2)

glmm212 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd212)
Anova(glmm212, type = 2)

glmm217 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd217)
Anova(glmm217, type = 2)


