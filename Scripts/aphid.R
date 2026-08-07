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
library(mgcv)

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
countdata$net_winged[is.na(countdata$net_winged)] <- 0

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
  mutate(focal_aphids = focal_wingless + focal_winged + net_winged, 
         focal_proportion = (focal_winged + net_winged)/focal_aphids)

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
view(df_dt)

# find average daytime temperatures for warmed and unwarmed on each day of experiment
df_dt_means <- df_dt %>% 
  group_by(day_number, date, warmingtreatment) %>% 
  summarise(meantempdaily = mean(meantemp),
            n=n(),
            sd = sd(meantemp),
            setemp = sd/sqrt(n))
view(df_dt_means)

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
    date = seq(min(date), max(date), by = "day")) %>%
  fill(ladybug_total, .direction = "up")

dailymeansdn <- df %>%
  group_by(day_number, cage, warmingtreatment) %>% #put warmingtreatment here so it keeps that column
  summarise(meantemp = mean(temperature),
            n=n(),
            sd = sd(temperature),
            setemp = sd / sqrt(n))
view(dailymeansdn)

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
  geom_line(aes(
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
    values = c("steelblue1", "red3")) +
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
    show.legend = FALSE) +
  geom_errorbar(data = df_dt_means,
    aes(x = date, y = meantempdaily,
      ymin = meantempdaily - setemp,
      ymax = meantempdaily + setemp,
      colour = warmingtreatment),
    width = 0) +
  geom_point(data = df_dt,aes(x = date, y = meantemp, colour = warmingtreatment),alpha = 0.2)
  
avgtemp

# plot treatment maximums thoughout growing season
maxtemp <- ggplot(df_dt_maxavg, aes(x = date, y = daymax)) +
  geom_point(
    aes(colour = warmingtreatment),
    size = 2) +
  geom_line(aes(colour = warmingtreatment,
      linetype = warmingtreatment,
      group = warmingtreatment),
    alpha = 0.5) +
  labs(x = "Date", y = "Mean maximum daily temperature (°C)") +
  theme_tess() +
  scale_color_manual(
    name = "Warming",
    labels = c("No", "Yes"),
    values = c("steelblue1", "red3")) +
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
                 "2025-09-02")) +
  geom_point(data = df_dt_max,aes(x = date, y = cagemax, colour = warmingtreatment),alpha = 0.2)

maxtemp

#### GLOBAL APHID ABUNDANCE OVER TIME ####

# make plot showing number of aphids on focal plant over time
focalplantcount <- ggplot(countdata_means, aes(x = date, y = mean_aphids, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Number of aphids on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = mean_aphids - se, ymax = mean_aphids + se), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata, aes(x = date, y = focal_aphids, group = cage, linetype = predatortreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "1 week"),
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_means, aes(x = date, y = mean_aphids, linetype = predatortreatment))
  
focalplantcount
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
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = mean_aphids - se, ymax = mean_aphids + se), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_zoom, aes(x = date, y = focal_aphids, group = cage, linetype = predatortreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "3 days"),
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_means_zoom, aes(x = date, y = mean_aphids, linetype = predatortreatment))
zoomedcount

# ensure jd is numeric and scaled for ease of analysis
countdata$jd <- as.numeric(countdata$jd)
countdata$jd_sc <- scale(countdata$jd)
# separate count into factor to use in ar1 (ordered variable)
countdata$count <- factor(countdata$count)

# model (original)
lmfullmodel<-lmer(focal_aphids ~ predatortreatment * warmingtreatment * jd_sc + predatortreatment * warmingtreatment * I(jd_sc^2) + (1 | cage), data = countdata)
abundancetable <- Anova(lmfullmodel, type=2)
abundancetable
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

#redoing glmm with ar1 as i believe ar1 is the correct way to do this now

glmm_nb <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment * jd_sc + predatortreatment * warmingtreatment * I(jd_sc^2) +
                     ar1(count + 0 | cage),
                   family = nbinom2,
                   data = countdata)

Anova(glmm_nb, type = 2)
# seeing sig date/inverse date, interactions with pred*date and warm*date

#trying GAMMs
countdata$trt <- interaction(countdata$warmingtreatment,
                             countdata$predatortreatment,
                             drop = TRUE)
countdata$cage <- factor(countdata$cage)

mod <- gam(focal_aphids ~ warmingtreatment * predatortreatment +
    s(jd_sc, by = trt, k = 6) +
    s(cage, bs = "re"),
  family = nb(), # figure out which is the best fit for this
  method = "REML",
  data = countdata)

summary(mod)
# sig effect of date, interaction date*pred and date*warm

# model comparisons/AICs for glmmTMBs
AIC(glmm_g, glmm_nb, glmm_p)
# negative binomial has lowest AIC value

# model comparisons/AICs for lmer and gamms
AIC(mod, lmfullmodel)
# gamms gives lower value

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
lm198 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd198)
Anova(lm198, type=2)
coef(summary(lm198))
#significant effect of warmingtreatment and an interaction with warmingxpred

#July 21
lm202 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd202)
Anova(lm202, type=2)
coef(summary(lm202))
#nothing significant

#July 24
lm205 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd205)
Anova(lm205, type=2)
coef(summary(lm205))
#nothing significant **changed to sig pred

#July 28
lm209 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2,  data = jd209)
Anova(lm209, type=2)
coef(summary(lm209))
#predatortreatment weakly significant **changed to sig pred

# July 31
lm212 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd212)
Anova(lm212, type=2)
coef(summary(lm212))
#predatortreatment marginally significant

#August 5
lm217 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd217)
Anova(lm217, type=2)
coef(summary(lm217))
#predatortreatment marginally significant

#August 7
lm219 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2,  data = jd219)
Anova(lm219, type=2)
coef(summary(lm219))
#nothing significant

#August 11
lm223 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd223)
Anova(lm223, type=2)
coef(summary(lm223))
#nothing significant

#August 14
lm226 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd226)
Anova(lm226, type=2)
coef(summary(lm226))
#nothing significant

#August 18
lm230 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd230)
Anova(lm230, type=2)
coef(summary(lm230))
#strong significant effect of warmingtreatment

#August 21
lm233 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd233)
Anova(lm233, type = 2)
coef(summary(lm233))
#significant effect of warmingtreatment

#August 25
lm237 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd237)
Anova(lm237, type = 2)
coef(summary(lm237))
#strong significant effect of warmingtreatment **changed + marginal pred

#August 28
lm240 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd240)
Anova(lm240, type = 2)
coef(summary(lm240))
#nothing significant

#Sept 2
lm245 <- glmmTMB(focal_aphids ~ predatortreatment * warmingtreatment, family = nbinom2, data = jd245)
Anova(lm245, type = 2)
coef(summary(lm245))
#weak significant effect of predatortreatment

# p-value correction (Benjamini-Hochberg)
# vector of raw p-values
raw_p_pred <- c(0.11323, 0.1553, 0.03073, 0.03498, 0.09931, 0.08799,0.1771, 0.9295, 0.9078, 0.418899, 0.39175, 0.066269, 0.3056, 0.1690)
raw_p_warm <- c(0.2964, 0.1736, 0.27136, 0.88243, 0.43058, 0.35704, 0.1165, 0.5420, 0.5874, 0.001401, 0.01598, 0.003379, 0.1266, 0.3271)

# Benjamini-Hochberg correction
adjusted_p_pred <- p.adjust(raw_p_pred, method = "BH")
adjusted_p_warm <- p.adjust(raw_p_warm, method = "BH")

# view the results
print(adjusted_p_pred)
print(adjusted_p_warm)


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
  mutate(focal_proportion = ((focal_winged + net_winged)/focal_aphids))%>%
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

view(countdata_means)

# plot winged/total
propwinged <- ggplot(countdata_means,
  aes(x = date,
      y = prop_winged,
      color = warmingtreatment,
      shape = predatortreatment)) +
  labs(x = "Date", y = "Proportion of winged aphids on focal plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = prop_winged - se_prop, ymax = prop_winged + se_prop),
    width = 0,
    position = position_dodge(width = 0.5)) +
  # lady beetle numbers added in (scaled to fit nicely with graph)
  #geom_step(data = countdata_filled,
   #         aes(x = date, y = ladybug_total * (0.5 / 16)), inherit.aes = FALSE, colour = "forestgreen", linewidth = 1, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = countdata,
    aes(x = date, y = focal_winged / focal_aphids, group = cage, linetype = predatortreatment),
    alpha = 0.1,
    position = position_dodge(width = 0.5)) +
  geom_line(data = countdata_means,
    aes(x = date, y = prop_winged, linetype = predatortreatment)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(breaks = seq(from = as.Date("2025-07-15"), to = as.Date("2025-09-02"), by = "1 week"), date_labels = "%b %d") #+
  # adding secondary axis + scaling
  #scale_y_continuous(limits = c(0, 0.5),
   # name = "Proportion of winged aphids on focal plant",
   # sec.axis = sec_axis(
  #    ~ . * 32,
   #   name = "Number of lady beetles",
  #    breaks = seq(0, 16, by = 4))) +
  #theme(axis.line.y.right = element_line(colour = "forestgreen"),
   # axis.ticks.y.right = element_line(colour = "forestgreen"),
   # axis.text.y.right = element_text(colour = "forestgreen"),
   # axis.title.y.right = element_text(colour = "forestgreen"))

propwinged

## analysis winged

countdata$cage <- as.factor(countdata$cage)

# changed back to a binomial model due to convergence issues from adding ar1 to this model
dispersal_model_binom <- glmmTMB(
  cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment * jd_sc +
    warmingtreatment * predatortreatment * I(jd_sc^2) + ar1(count + 0 | cage),
  family = binomial(),
  data = countdata)
Anova(dispersal_model_binom, type = 2)

#### DAY BY DAY PROPORTION WINGED ANALYSIS ####

## analysis winged by julian day

#July 17th 
glmm198 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                   family = betabinomial(), data = jd198)
Anova(glmm198,type = 2)
fixef(glmm198)
boxplot(focal_proportion~warmingtreatment, data=jd198)
# weak effect of predatortreatment, strong effect of warmingtreatment (more winged in warmed)
# checking for overdispersion
testDispersion(glmm198)
# not overdispersed

#July 21st 
glmm202 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
               family = betabinomial(),  data = jd202)
Anova(glmm202,type = 2)
fixef(glmm202)
# nothing significant
# checking for overdispersion
testDispersion(glmm202)
# not overdispersed

#July 24th 
glmm205 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd205)
Anova(glmm205,type = 2)
fixef(glmm205)
# nothing significant
# checking for overdispersion
testDispersion(glmm205)
# overdispersed

#July 28th
glmm209 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd209)
Anova(glmm209,type = 2)
fixef(glmm209)
boxplot(focal_proportion~warmingtreatment, data=jd209)
boxplot(focal_proportion~predatortreatment, data=jd209)
# strong effect of warmingtreatment (higher in warmed), effect of predatortreatment (lower in pred), no interaction
# checking for overdispersion
testDispersion(glmm209)
# overdispersed

#July 31st
glmm212 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd212)
Anova(glmm212,type = 2)
fixef(glmm212)
boxplot(focal_proportion~warmingtreatment, data=jd209)
boxplot(focal_proportion~predatortreatment, data=jd209)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm212)
# overdispersed

#August 5th 
glmm217 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd217)
Anova(glmm217,type = 2)
fixef(glmm217)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm217)
# overdispersed

#August 7th 
glmm219 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd219)
Anova(glmm219,type = 2)
fixef(glmm219)
boxplot(focal_proportion~warmingtreatment*predatortreatment, data=jd219)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm219)
# overdispersed

#August 11
glmm223 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd223)
Anova(glmm223, type = 2)
fixef(glmm223)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm223)
# overdispersed

#August 14
glmm226 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd226)
Anova(glmm226,type = 2)
fixef(glmm226)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm226)
# overdispersed

# August 18
glmm230 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd230)
Anova(glmm230,type = 2)
fixef(glmm230)
# strong interactin and predator effect
# checking for overdispersion
testDispersion(glmm230)
# overdispersed

# August 21
glmm233 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd233)
Anova(glmm233,type = 2)
fixef(glmm233)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm233)
# overdispersed

# August 25
glmm237 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial, data = jd237)
Anova(glmm237,type = 2)
fixef(glmm237)
# weak effect pred, strong effect warmingtreatment
# checking for overdispersion
testDispersion(glmm237)
# overdispersed

# August 28
glmm240 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd240)
Anova(glmm240,type = 2)
fixef(glmm240)
# strong effect of everything
# checking for overdispersion
testDispersion(glmm240)
# overdispersed

# Sept 2
glmm245 <- glmmTMB(cbind(focal_winged + net_winged, focal_wingless) ~ warmingtreatment * predatortreatment,
                 family = betabinomial,data = jd245)
Anova(glmm245,type = 2)
fixef(glmm245)
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
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = dispersed_mean - se_dispersed, ymax = dispersed_mean + se_dispersed), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_dispersed, aes(x = date, y = dispersed, group = cage, linetype = predatortreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_dispersed_means, aes(x = date, y = dispersed_mean, linetype = predatortreatment), position = position_dodge(width = 0.5))

dispersedplot

# run lmm
dispersed_model <- lmer(dispersed ~ warmingtreatment * predatortreatment * jd + (1 | cage), data = countdata_dispersed)
Anova(dispersed_model, type = 2)

testDispersion(dispersed_model)

# plot dispersal to sentinel plant (success) vs to sticky card (failure)
sentineldispersedplot <- ggplot(countdata_dispersed_means, aes(x = date, y = sentinel_dispersed_mean, color = warmingtreatment, shape = predatortreatment)) +
  labs(x = "Date", y = "Proportion of dispersed aphids to sentinel plant") +
  theme_tess() +
  scale_color_manual(name = "Warming", labels = c("No", "Yes"), values = c("steelblue1", "red3")) +
  scale_shape_manual(name = "Predator", labels = c("No", "Yes"), values = c(1, 16)) +
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_errorbar(aes(ymin = sentinel_dispersed_mean - se_sentinel_dispersed, ymax = sentinel_dispersed_mean + se_sentinel_dispersed), width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(data = countdata_dispersed, aes(x = date, y = sentinel_dispersed, group = cage, linetype = predatortreatment), alpha = 1/10, position = position_dodge(width = 0.5)) +
  scale_linetype_manual(name = "Predator", labels = c("No", "Yes"), values = c("dashed", "solid")) +
  scale_x_date(
    date_labels = "%b %d"
  ) +
  geom_line(data = countdata_dispersed_means, aes(x = date, y = sentinel_dispersed_mean, linetype = predatortreatment), position = position_dodge(width = 0.5))

sentineldispersedplot

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
## using glmmTMB package which allows for quasi/betabinomials and fixes most overdispersion
glmm205 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd205)
Anova(glmm205, type = 2)
fixef(glmm205)

glmm209 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment, 
                 family = betabinomial, data = djd209)
Anova(glmm209, type = 2)
fixef(glmm209)

glmm212 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd212)
Anova(glmm212, type = 2)
fixef(glmm212)

glmm217 <- glmmTMB(cbind(sentinel_aphids, sticky_aphids)  ~ predatortreatment * warmingtreatment,
                 family = betabinomial, data = djd217)
Anova(glmm217, type = 2)
fixef(glmm217)

# total dispersal
glmm205 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd205)
Anova(glmm205, type = 2)

glmm209 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd209)
Anova(glmm209, type = 2)

glmm212 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd212)
Anova(glmm212, type = 2)

glmm217 <- lm(dispersed ~ warmingtreatment * predatortreatment, data = djd217)
Anova(glmm217, type = 2)


