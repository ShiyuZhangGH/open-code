library(tidyverse)

# ========= SECTION: PROXIMAL OUTCOME =======

#READ: eligible decision points 
dat_mars_29_elig<- read_rds("d1.rds")
dim(dat_mars_29_elig)

# Y (descriptive)----

## group by decision point (person-block), count # of missing

nmiss_by_dp<- dat_mars_29_elig%>%
  group_by(decision_point)%>%
  summarise(nelig= n(), nmiss= sum(Y_miss==1))

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





# Prediction ----

# At the participant level, do individual characteristics predict the percentage of missingness in Y?

nmiss_by_person<- dat_mars_29_elig%>%
  group_by(mars_id)%>%
  summarise(nelig= n(), nmiss= sum(Y_miss==1), 
            syn_race_ethnicity= first(syn_race_ethnicity), syn_tobacco_history= first(syn_tobacco_history))


# outcome: percentage of missing out of eligible decision points
nmiss_by_person$per_missY_elig<- nmiss_by_person$nmiss/ nmiss_by_person$nelig


# modeling outcome

lm1<- lm(per_missY_elig~  syn_race_ethnicity, data= nmiss_by_person)
summary(lm1)


lm1<- lm(per_missY_elig~  factor(syn_tobacco_history), data= nmiss_by_person)
summary(lm1)


# ========= SECTION: COVARIATES =======

# READ: eligible decision points
dat_for_analysis_elig<- read_rds("d2.rds") 
dim(dat_for_analysis_elig)

dat_for_analysis_elig$miss_covariate_sum<- dat_for_analysis_elig$miss_hour_coinflip_local+
  dat_for_analysis_elig$miss_days_between_v1_and_coinflip_local+
  dat_for_analysis_elig$miss_any_response_2qs+
  dat_for_analysis_elig$miss_any_recent_eligible_dp+
  dat_for_analysis_elig$miss_engagement_most_recent_eligible+
      dat_for_analysis_elig$miss_age+
  dat_for_analysis_elig$miss_is_female+
  dat_for_analysis_elig$miss_is_latino+
  dat_for_analysis_elig$miss_is_not_latino_and_black+
  dat_for_analysis_elig$miss_is_not_latino_and_other+
  dat_for_analysis_elig$miss_baseline_tobacco_history+
  dat_for_analysis_elig$miss_has_partner+
  dat_for_analysis_elig$miss_income_val

dat_for_analysis_elig%>%
  filter(Y_miss==0)%>% #among decision points that have a non-missing Y
  select(miss_covariate_sum)%>% #do they have missing values in any of the covariates? 
  table()


# ========= SECTION: ELIGIBILITY =======

# READ: all possible decision points between day 2-9
dat_analysis<- read_rds("d3.rds")
dim(dat_analysis)

# The eligbility=1 could be because participants were truly eligible at those moments,
# or because their eligibility information was missing
table(dat_analysis$sys_info_conditions,dat_analysis$availability, exclude= F)



# ========= SECTION: EMBEDDED TAILORING VARIABLE THAT DO NOT RESTRICT RANDOMIZATIONS =======

# READ: decision points that are randomized to the low-effort condition
dat_mars_29_lowEff<- read_rds( "d4.rds") 
dim(dat_mars_29_lowEff)

table(dat_mars_29_lowEff$miss_cig_available, dat_mars_29_lowEff$miss_negative_affect, exclude= FALSE)



# ========= SECTION: CANDIDATE TAILORING VARIABLE ============

# READ: decision points that have a non-missing Y value
dat_Y<- read_rds("d5.rds")
dim(dat_Y)

table(dat_Y$miss_ashamed_lag1)
table(dat_Y$miss_guilty_lag1)
table(dat_Y$miss_happy_lag1)
table(dat_Y$miss_stressor_is_any_lag1)
table(dat_Y$miss_src_scored_lag1)
table(dat_Y$miss_substance_is_any_lag1)
table(dat_Y$miss_substance_is_any_nicotine_lag1)
table(dat_Y$miss_substance_is_alcohol_lag1)
table(dat_Y$miss_substance_is_marijuana_or_cannabis_lag1)
table(dat_Y$miss_cigarette_counts_lag1)
table(dat_Y$miss_alcohol_counts_lag1)




