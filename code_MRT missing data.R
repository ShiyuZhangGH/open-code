library(clipr)
library(MRTAnalysis)
library(tidyverse)


################################################################################################################################################################

setwd("C:/Users/zsy/Documents/MARS/Jamie_code/MARS-Repository-Primary-Aim-master")
source("paths.R")
source("prepare-analytic-datasets/functions-for-analytical-dataset-creation.R")


################################################################################
# Load datasets
################################################################################
# This dataset contains data from only those sequences which began
# within the 10-day MRT period.
dat_matched_to_decision_points <- readRDS(file = file.path(path_manipulated_data, "dat_matched_to_decision_points.rds"))
# The following participants were excluded from all analyses:
# - Pilot participants
# - Among participants who have any data from sequences which began within the
#   10-day MRT period, those participants who did not complete
#   at least 3 EMA between the 2nd and 9th day inclusive.
mars_ids_excluded_from_all_analytic_datasets <- readRDS(file = file.path(path_manipulated_data, "mars_ids_excluded_from_all_analytic_datasets.rds"))

################################################################################
# First, drop data from participants which we will not be using in any further
# analysis.
################################################################################
dat_mars_analysis <- dat_matched_to_decision_points %>% filter(!(mars_id %in% mars_ids_excluded_from_all_analytic_datasets))

################################################################################
# Next, create the primary proximal outcome using EMA responses
################################################################################
dat_mars_analysis[["Y"]] <- construct_primary_proximal_outcome(cleaned_data_frame = dat_mars_analysis, q1_var_name = "Q1_response", q2_var_name = "Q2_response", q3_var_name = "Q3_response")

################################################################################
# Next, create other variables needed in data analysis
################################################################################
# Participant ID must be is numeric format (not character) when performing analysis
dat_mars_analysis <- dat_mars_analysis %>% mutate(participant_id = as.numeric(substring(mars_id, first = 6))) 
# Treatment indicators must be in numeric format (not character) when performing analysis
dat_mars_analysis <- dat_mars_analysis %>% 
  mutate(eligibility = if_else(!is.na(A), 1, 0)) %>%
  mutate(coinflip = if_else(A != "none", 1, 0),
         is_high_effort = if_else(A == "mars", 1, 0),
         is_low_effort = if_else(A == "low_effort", 1, 0))

################################################################################################################################################################

# Filter----

## person-blocks within 2-9 days
dat_mars_29<- dat_mars_analysis%>% filter(!cluster_id%in% c(1, 10))

## person-blocks eligible for treatment (not driving, not sleep mode, phone on)
dat_mars_29_elig<- dat_mars_29%>% filter(!is.na(A))


# Y (descriptive)----

## group by decision point (person-block), count # of missing

nmiss_by_dp<- dat_mars_29_elig%>%
  group_by(decision_point)%>%
  summarise(nelig= n(), nmiss= sum(is.na(Y)))

nmiss_by_dp$ncom<- nmiss_by_dp$nelig- nmiss_by_dp$nmiss
nmiss_by_dp$ninelig<- 99- nmiss_by_dp$nelig


nmiss_by_dp_long<- nmiss_by_dp%>%
  pivot_longer(cols= nmiss:ninelig,
               names_to= "response",
               values_to= "n")
nmiss_by_dp_long$response<- factor(nmiss_by_dp_long$response, 
                                   levels= c("nmiss", "ninelig", "ncom"))




dat_helper<- nmiss_by_dp_long%>%
  filter(response== "ncom")%>%
  mutate(complete_rate= n/nelig)

ggplot(dat_helper)+ # Figure 2
  geom_bar(aes(x= decision_point, y= complete_rate),
           position="stack", stat="identity")+
  geom_hline(yintercept= 1)+
  annotate(geom= "text", x= 21, y= 1.04, label= "quit day",  color= "red4")+
  annotate(geom= "segment", x= 19, xend= 24, y= 1.01, yend= 1.01, color= "red4")+
  
  scale_x_continuous(name= "decision point",
                     breaks= c(7, 13, 19, 25, 31, 37, 43, 49),
                     labels= c("7"= paste("7", "day2", sep= "\n"),
                               "13"= paste("13", "day3", sep= "\n"),
                               "19"= paste("19", "day4", sep= "\n"), 
                               "25"= paste("25", "day5", sep= "\n"),
                               "31"= paste("31", "day6", sep= "\n"),
                               "37"= paste("37", "day7", sep= "\n"),
                               "43"= paste("43", "day8", sep= "\n"),
                               "49"= paste("49", "day9", sep= "\n")))+
  scale_y_continuous(name= "proportion of completed assessments")




# Y Prediction ----

# At the participant level, do individual characteristics predict the percentage of missingness in Y?

nmiss_by_person<- dat_mars_29_elig%>%
  group_by(mars_id)%>%
  summarise(nelig= n(), nmiss= sum(is.na(Y)))


# outcome 2: percentage of missing out of eligible decision points
nmiss_by_person$per_missY_elig<- nmiss_by_person$nmiss/ nmiss_by_person$nelig
summary(nmiss_by_person$per_missY_elig)

# baseline predictors
dat_demogs <- readRDS(file = file.path(path_manipulated_data, "dat_demogs.rds"))


dat_demogs<- dat_demogs%>%
  mutate(
    female= case_when(
      gender_category== "female"~ 1,
      gender_category== "male"~ 0
    ),
    
    race_ethnicity= case_when(
      race_and_ethnicity== "not latino and black"~ "non-l black",
      race_and_ethnicity== "latino"~ "latino",
      race_and_ethnicity%in% c("not latino and white",
                               "other")~ "non-l other"
    ),
    race_ethnicity= factor(race_ethnicity,
                           levels= c("non-l other",
                                     "non-l black",
                                     "latino")),  
    
    income_cat3= case_when(
      income_category%in% c("less than or equal to USD 9,999",
                            "greater than USD 9,999 and less than or equal to USD 19,999",
                            "greater than USD 19,999 and less than or equal to USD 29,999")~ 1,
      income_category%in% c("greater than USD 29,999 and less than or equal to USD 39,999",
                            "greater than USD 39,999 and less than or equal to USD 49,999",
                            "greater than USD 49,999 and less than or equal to USD 59,999",
                            "greater than USD 59,999 and less than or equal to USD 69,999")~ 2,
      income_category%in% c("greater than USD 69,999 and less than or equal to USD 79,999",
                            "greater than USD 79,999 and less than or equal to USD 89,999",
                            "greater than USD 89,999 and less than or equal to USD 99,999",
                            "greater than or equal to USD 100,000")~ 3
    ),
    
    marital= case_when(
      partner_status_category== "married"~ "married",
      partner_status_category== "single"~ "single",
      partner_status_category== "living with significant other"~ "cohabit",
      partner_status_category%in% c("widowed", "divorced", "separated")~ "w/d/s"
    ),
    marital= factor(marital, 
                    levels= c("single", "married", "cohabit", "w/d/s"))
  )


# merging
nmiss_by_person<- nmiss_by_person%>%
  left_join(dat_demogs%>% select(mars_id, age, female, race_ethnicity, income_cat3, marital, baseline_tobacco_history),
            by= "mars_id")




# modeling - outcome 2
summary(nmiss_by_person$per_missY_elig)

lm1<- lm(per_missY_elig~  age, data= nmiss_by_person)
summary(lm1)

lm1<- lm(per_missY_elig~  female, data= nmiss_by_person)
summary(lm1)

lm1<- lm(per_missY_elig~  race_ethnicity, data= nmiss_by_person)
summary(lm1)
summary(lm1)$coefficients%>% data.frame()%>% write_clip()

lm1<- lm(per_missY_elig~  factor(income_cat3), data= nmiss_by_person)
summary(lm1)

lm1<- lm(per_missY_elig~  marital, data= nmiss_by_person)
summary(lm1)

lm1<- lm(per_missY_elig~  baseline_tobacco_history, data= nmiss_by_person)
summary(lm1)
summary(lm1)$coefficients%>% data.frame()

#=================================================================================
#=================================================================================
# Figure 3

# This is the data frame that is the output of the very last step of our data pipeline
dat_primary_aim <- readRDS(file = file.path(path_manipulated_data, "dat_primary_aim.rds")) 


# How the analytical sample size changes 

# 99 participants * 10 days * 6 decision point/day

n<- 99*10*6

# only day 2-9
n_29<- 99*(10-2)*6

# only eligible
n_29_elig<- dat_primary_aim%>%
  filter(decision_point>= 7 & decision_point<= 54)%>%
  filter(eligibility== 1)%>%
  nrow()

# valid Y
n_29_elig_Y<- dat_primary_aim%>%
  filter(decision_point>= 7 & decision_point<= 54)%>%
  filter(eligibility== 1)%>%
  filter(!is.na(Y))%>%
  nrow()

# Responded to EMA at time (t-1), i.e., may be included in moderation analysis
dat_primary_aim<- dat_primary_aim%>%
  group_by(mars_id)%>%
  mutate(
    status_survey_ema_lag1= lag(status_survey_ema_collapsed),
    substance_is_any_lag1 = lag(substance_is_any) # an example moderator from EMA
  )%>%
  ungroup()%>% data.frame()%>%
  
  mutate(
    ema_respondent_lag1= case_when(
      status_survey_ema_lag1%in% c("fully_completed", "partially_completed")~ 1,
      TRUE~ 0
    )
  )

n_29_elig_Y_preEMA<- dat_primary_aim%>%
  filter(decision_point>= 7 & decision_point<= 54)%>%
  filter(eligibility== 1)%>%
  filter(!is.na(Y))%>%
  filter(ema_respondent_lag1==1)%>%
  nrow()



# ggplot
n
n_29
n_29_elig 
n_29_elig_Y
n_29_elig_Y_preEMA

d<- data.frame(
  helper= c("4", "3", "2", "1"),
  
  filter= c("day 2-9",
            "exclude ineligible decision points",
            "exclude missing in Y",
            "excluding missing in time (t-1)"),
  n= c(n_29, 
       n_29_elig, 
       n_29_elig_Y, 
       n_29_elig_Y_preEMA)
)


ggplot(d, aes(x= helper, y= n, fill= helper))+
  geom_bar( stat="identity")+
  coord_flip()+
  scale_fill_manual(values= c("4"= "gray90",
                              "3"= "gray70",
                              "2"= "gray35",
                              "1"= "gray20"))+
  #geom_hline(yintercept= 4752)+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        
        axis.line = element_line(colour = "white"),
        legend.position = "none")


