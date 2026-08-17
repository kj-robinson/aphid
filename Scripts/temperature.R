#load packages

library(tidyverse)
library(lubridate)
library(cowplot)
library(lme4)

#load my ggplot theme to make plots look nice
theme_tess <- function () { 
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.text.x=element_text(size=20))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.x=element_text(size=20))+
    theme(axis.title.y=element_text(size=20))+
    theme(plot.title = element_text(hjust = 0.5,size=25))+
    theme(legend.title = element_text(size = 18))+
    theme(legend.title=element_text(size=16,face="bold"))+
    theme(legend.text=element_text(size = 16))+
    theme(strip.background = element_blank())+
    theme(strip.text = element_blank())
}

#load data
data<- read.csv("./Data/temperature.csv")
view(data)
data <- data[!is.na(data$temperature),]

#separate data and time into two columns
df <- data %>%
  mutate(datetime = str_squish(datetime)) %>% #remove extra white spaces
  separate(datetime, into = c("date", "time"), sep = " ") 
view(df)

#make a new column that specifies if a a row is in the day or night - day is 7:00 to 19:00
df <- df %>%
  mutate(
    time_parsed = hms(time),
    hour = hour(time_parsed),
    minute = minute(time_parsed),
    period = case_when(
      hour > 7 & hour < 19 ~ "day",
      hour == 7 & minute > 0 ~ "day",
      hour == 19 & minute == 0 ~ "day",
      TRUE ~ "night"
    )
  )

#make a column for date that goes 1 --> last date
df <- df %>%
  mutate(date = mdy(date),                # convert string to Date format (month-day-year)
    day_number = dense_rank(date))        # assign day numbers starting from earliest

#plot all the data with each cage as a line

cols <- c("unwarmed" = "steelblue1", "warmed" = "red3")

p<-ggplot(df, aes(x = time, y = temperature, color = warmingtreatment, group=cage)) +
  geom_line()+
  facet_wrap(~ day_number) + #make a separate panel for each day
  scale_color_manual(values = cols,
                     name = "Warming", 
                     labels=c("unwarmed", "warmed"), 
                     breaks = c("Unwarmed","Warmed")) +
  labs(x = "Time",y = "Temperature (°C)") +
  theme_tess()
p

windows();p

####FULL DAY (DAY + NIGHT) CALCULATIONS ####

#calculate means for each cage for each day (the daily average temp for each cage) 
#this gives you one value per replicate (replicate = cage)

dailymeansdn <- df %>%
  group_by(day_number, cage, warmingtreatment) %>% #put warmingtreatment here so it keeps that column
  summarise(meantemp = mean(temperature),
            n=n(),
            sd = sd(temperature),
            setemp = sd / sqrt(n))

#now use these daily averages to calculate the means and SE for each treatment (warmed vs. unwarmed)

treatmentmeansdn <- dailymeansdn %>%
  group_by(warmingtreatment) %>%
  summarise(mean = mean(meantemp),
            n=n(),
            sd = sd(meantemp),
            se = sd / sqrt(n))
view(treatmentmeansdn)

#plot this

p <- ggplot() +
  geom_point(data = treatmentmeansdn,
             aes(x =warmingtreatment, y = mean), size=3)+
  geom_errorbar(data = treatmentmeansdn,
                aes(x = warmingtreatment,
                    ymin = mean - se,
                    ymax = mean + se), width = 0, linewidth=1.2) +
  labs(x = "",y = "Temperature (°C)") +
  scale_y_continuous(limits=c(22.5,25.5))+
  scale_x_discrete(labels=c("Unwarmed", "Warmed"))+
  theme_tess()
p
  
#### JUST DAYTIME CALCULATIONS ####

dayonly<-df%>%
  filter(period=="day")

#calculate means for each cage for each day (the daily average temp for each cage) 

dailymeansd<-dayonly%>%
  group_by(day_number,cage,warmingtreatment) %>%
  summarise(meantemp = mean(temperature),
            n=n(),
            sd = sd(temperature),
            setemp = sd / sqrt(n))

#now use these daily averages to calculate the means and SE for each treatment (warmed vs. unwarmed)

treatmentmeansd<-dailymeansd%>%
  group_by(warmingtreatment) %>%
  summarise(mean = mean(meantemp),
            n=n(),
            sd = sd(meantemp),
            se = sd / sqrt(n))
view(treatmentmeansd)

#plot this

p <- ggplot() +
  geom_point(data = treatmentmeansd,
             aes(x =warmingtreatment, y = mean), size = 3)+
  geom_errorbar(data = treatmentmeansd, 
                aes(x = warmingtreatment,
                    ymin = mean - se,
                    ymax = mean + se), width = 0, linewidth=1.2) +
  labs(x = "",y = "Temperature (°C)") +
  scale_y_continuous(limits=c(29.5,33.5))+
  scale_x_discrete(labels=c("Unwarmed", "Warmed"))+
  theme_tess()

p

#### MAX AND AVG TEMPERATURES ####

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
  labs(x = "Date", y = "Temperature (°C)") +
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
    date_labels = "%b %d") +
  scale_y_continuous(limits = c(12, 34))

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
  labs(x = "Date", y = "Temperature (°C)") +
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
    date_labels = "%b %d")

maxtemp
