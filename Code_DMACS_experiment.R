library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(survey)
library(PracTools)
library(clipr) #copy table to clipboard
library(naniar) #replace na
library(ggplot2)
library(plyr) #rounding to 0.01
library(tidycensus)
library(tigris)
library(sf)
library(mice)


# ====Reading datasets====
datRes<- read_dta("W12 Data as of June 11_Detroiters Only.dta") #Revised dataset after race bug

# + ACS point of refernce ====
detroit_tract_ID #this is a variable storing Detroit census tract ID
# Download ACS data as benchmark
acs19Var<- load_variables(2019, "acs1", cache= TRUE)
#downloading Wayne tract shape from census
census_tract<- tracts(state= "MI", county = "Wayne") 
#subset to Detroit tract
census_tract_detroit<- census_tract[census_tract$GEOID%in% detroit_tract_ID,]

acs19<- get_acs(geography = "tract",
                state= "MI",
                county = "Wayne",
                year = 2019,
                variables = c(totPop= "B01001_001",
                              totHouseholder= "B07013_001",
                              totPop_Mless5= "B01001_003",
                              totPop_M5_9= "B01001_004",
                              totPop_M10_14= "B01001_005",
                              totPop_M15_17= "B01001_006",
                              totPop_Fless5= "B01001_027",
                              totPop_F5_9= "B01001_028",
                              totPop_F10_14= "B01001_029",
                              totPop_F15_17= "B01001_030",
                              
                              lanDenominator= "B06007_001", #pop Under 5
                              englishOnly= "B06007_002",
                              spanish= "B06007_003",
                              otherLan= "B06007_006",
                              ownerHouseholder= "B07013_002",
                              renterHouseholder= "B07013_003",
                              
                              computerHHDenominator= "B28001_001",
                              anyComputerHH= "B28001_002",
                              desklaptopHH= "B28001_003",
                              smartphoneOnlyHH= "B28001_006",
                              otherComputerOnlyHH= "B28001_010",
                              
                              computerPerDenominator_18_64= "B28005_008",
                              computerPer_18_64= "B28005_009",
                              internetPer1_18_64= "B28005_010",
                              internetPer2_18_64= "B28005_011",
                              noInternetPer_18_64= "B28005_012",
                              noComputerPer_18_64= "B28005_013",
                              
                              computerPerDenominator_over65= "B28005_014",
                              computerPer_over65= "B28005_015",
                              internetPer1_over65= "B28005_016",
                              internetPer2_over65= "B28005_017",
                              noInternetPer_over65= "B28005_018",
                              noComputerPer_over65= "B28005_019",
                              
                              maritalDenominator= "B12002_001",
                              maritalDenominatorM= "B12002_002",
                              neverMarrM= "B12002_003",
                              nowMarrM= "B12002_018",
                              widowM= "B12002_065",
                              divorceM= "B12002_080",
                              
                              neverMarrF= "B12002_096",
                              nowMarrF= "B12002_111",
                              widowF= "B12002_158",
                              divorceF= "B12002_173",
                              
                              insureDenominator_19_34= "B27010_018",
                              noInsure_19_34= "B27010_033",
                              
                              insureDenominator_35_64= "B27010_034",
                              noInsure_35_64= "B27010_050",
                              
                              insureDenominator_over65= "B27010_051",
                              noInsure_over65= "B27010_066"
                ))

acs19_wide<- pivot_wider(acs19, names_from= "variable", 
                         values_from= "estimate", id_cols= "GEOID")
acs19_wide<- acs19_wide[acs19_wide$GEOID%in%detroit_tract_ID,] #subsetting Detroit


#householders lived in owner-occupied units 
ref_owner<- sum(acs19_wide$ownerHouseholder)/sum(acs19_wide$totHouseholder) #Proportion

#population who speak other language 
totPop_over5= sum(acs19_wide$totPop)- 
  (sum(acs19_wide$totPop_Mless5)+sum(acs19_wide$totPop_Fless5))
ref_language<- (sum(acs19_wide$spanish)+ sum(acs19_wide$otherLan))/totPop_over5 #Proportion

#Computer 
haveComputer= sum(acs19_wide$computerPer_18_64)+ sum(acs19_wide$computerPer_over65)
denominator= sum(acs19_wide$computerPerDenominator_18_64)+ sum(acs19_wide$computerPerDenominator_over65)
haveComputer/denominator
#problem with computer is that ACS asks about any type of computer including smartphone
#but people typically dont consider smartphone as computer
#Here is rough estimate to exclude smartphone. 
#The computer type variable is at the household level
#Of hh that reports having computing device, % that have desktop or laptop
desklaptopPercent= sum(acs19_wide$desklaptopHH)/ sum(acs19_wide$anyComputerHH)
#Of hh that reports having computing device, % that have desktop, laptop or tablet
nominator= sum(acs19_wide$anyComputerHH)- sum(acs19_wide$smartphoneOnlyHH)
desklapTabletPercent= nominator/sum(acs19_wide$anyComputerHH)
#adjusted have computer
haveComputer/denominator*desklaptopPercent #Proportion
ref_computer<- haveComputer/denominator*desklapTabletPercent #Proportion


#Internet
nominator= sum(acs19_wide$computerPer_18_64)- sum(acs19_wide$noInternetPer_18_64)+ 
  sum(acs19_wide$computerPer_over65)- sum(acs19_wide$noInternetPer_over65)
denomator= sum(acs19_wide$computerPerDenominator_18_64)+ sum(acs19_wide$computerPerDenominator_over65)
ref_internet<- nominator/denomator #Proportion

#divorced 
totPop_over15<- sum(acs19_wide$totPop)- 
  (sum(acs19_wide$totPop_Mless5)+sum(acs19_wide$totPop_Fless5)+
     sum(acs19_wide$totPop_M5_9)+sum(acs19_wide$totPop_F5_9)+
     sum(acs19_wide$totPop_M10_14)+sum(acs19_wide$totPop_F10_14))

ref_never_marry<- (sum(acs19_wide$neverMarrF)+ sum(acs19_wide$neverMarrM))/sum(acs19_wide$maritalDenominator)
ref_marry<- (sum(acs19_wide$nowMarrM)+ sum(acs19_wide$nowMarrF))/sum(acs19_wide$maritalDenominator)
ref_divorce<- (sum(acs19_wide$divorceM)+ sum(acs19_wide$divorceF))/sum(acs19_wide$maritalDenominator)

#insurance
nominator= sum(acs19_wide$noInsure_19_34)+ 
  sum(acs19_wide$noInsure_35_64)+ 
  sum(acs19_wide$noInsure_over65)
denominator= sum(acs19_wide$insureDenominator_19_34)+ 
  sum(acs19_wide$insureDenominator_35_64)+ 
  sum(acs19_wide$insureDenominator_over65)
ref_insured<- nominator/denominator




# ====Managing the response dataset====

# + Demographic features ====
datRes<- replace_with_na_at(datRes,
                            .vars= c("ageinyears_d12", "educat6", 
                                     "income_d12", "income_1_d12",
                                     "income_2_d12", "income_3_d12"),
                            condition= ~.x== -99)

datRes<- datRes%>%
  mutate(
    gender= dplyr::case_when(
      gender_d12%in% c(1, 4)~ "man", #men and trans men
      gender_d12%in% c(2, 3)~ "woman" #women and trans women
      #everyone else coded as NA, including fluid identity
    ))

datRes<- datRes%>%
  mutate(
    gender_age= dplyr::case_when(
      gender=="man" & agecat4== 1 ~ "man<35",
      gender=="man" & agecat4== 2 ~ "man35-54",
      gender=="man" & agecat4== 3 ~ "man55-64",
      gender=="man" & agecat4== 4 ~ "man>65",
      gender=="woman" & agecat4== 1 ~ "woman<35",
      gender=="woman" & agecat4== 2 ~ "woman35-54",
      gender=="woman" & agecat4== 3 ~ "woman55-64",
      gender=="woman" & agecat4== 4 ~ "woman>65",
      
    ))
datRes$gender_age<- factor(datRes$gender_age,
                           levels= c("man<35", "man35-54", 
                                     "man55-64", "man>65",
                                     "woman<35", "woman35-54", 
                                     "woman55-64", "woman>65"))

# ***** Imputing demographic features====
# There are missing values in the demographic features, this creates problem
# for weighting. 

# Re-creating demographic variables These variables will be replaced with imputed values
datRes$gender_impute<- datRes$gender%>%
  dplyr::recode("man"= 1, "woman"= 0, .missing= NA_real_)
datRes$gender_impute<- factor(datRes$gender_impute)
datRes$age_impute<- datRes$ageinyears_d12
datRes$racecat5_impute<- factor(datRes$racecat5)
datRes$educat6_impute<- factor(datRes$educat6)
datRes$incomecat5_impute<- factor(datRes$incomecat5)
datRes$income_1_enrich<- factor(datRes$income_1_d12)

impute_vari<- c("gender_impute", "age_impute", 
                "racecat5_impute", "educat6_impute", 
                "incomecat5_impute")
help_vari<- c("income_1_enrich",
              "percent_nonhisp_black", #ACS
              "percent_hisp", 
              "percent_lessthan10000",
              "percent_10kto29k",
              "percent_30kto49k",
              "percent_50kto99k")
#imputation
tempData <- mice(datRes[c(impute_vari, help_vari)], 
                 m=1,
                 maxit=20,
                 meth= c("logreg", "pmm", 
                         "polyreg", "polyreg", "polyreg",
                         "logreg",
                         "pmm", "pmm", "pmm", 
                         "pmm", "pmm", "pmm"),seed=11)
#assinging the imputed variables back to the dataset
datRes[c(impute_vari, help_vari)]<- complete(tempData)

#create a imputed agecat4 variable
datRes<- datRes%>% 
  mutate(agecat4_impute= dplyr::case_when(
    age_impute< 35~ 1,
    age_impute>= 35 & age_impute< 55~ 2,
    age_impute>= 55 & age_impute< 65~ 3,
    age_impute>= 65~ 4
  ))

#Then, create an imputed (and numeric) version of the gender-age variable
table(datRes$agecat4_impute, exclude = F)
datRes<- datRes%>%
  mutate(
    gender_age_impute= dplyr::case_when(
      gender_impute==1 & agecat4_impute== 1 ~ 1, #"man<35",
      gender_impute==1 & agecat4_impute== 2 ~ 2, #"man35-54",
      gender_impute==1 & agecat4_impute== 3 ~ 3, #"man55-64",
      gender_impute==1 & agecat4_impute== 4 ~ 4, #"man>65",
      gender_impute==0 & agecat4_impute== 1 ~ 5, #"woman<35",
      gender_impute==0 & agecat4_impute== 2 ~ 6, #"woman35-54",
      gender_impute==0 & agecat4_impute== 3 ~ 7, #"woman55-64",
      gender_impute==0 & agecat4_impute== 4 ~ 8 #"woman>65",
    ))

datRes<- datRes%>%
  mutate(educat4_impute= dplyr::case_when(
    educat6_impute==1~ 1,
    educat6_impute==2~ 2,
    educat6_impute%in% c(3, 4)~ 3,
    educat6_impute%in% c(5, 6)~ 4
  ))

datRes$gender_age_impute<- as.factor(datRes$gender_age_impute)
datRes$educat6_impute<- as.factor(datRes$educat6_impute)
datRes$educat4_impute<- as.factor(datRes$educat4_impute)
datRes$racecat5_impute<- as.factor(datRes$racecat5_impute)
datRes$incomecat5_impute<- as.factor(datRes$incomecat5_impute)



# + Key survey variables ====
#setting -99 to NA
vari<- c("housing_d12", "home_owner_d12",
         "nb_reputation_d12", "nb_satis_d12", 
         "vac_covid_d12", "trust_doctor_d12", "trust_faith_leader_d12",
         "trust_friends_d12", "trust_acquaint_d12", "trust_news_d12",
         "trust_socialmedia_d12", "trust_usgovt_d12", "trust_taskforce_d12",
         "language_d12",
         "walk_safety_d12", 
         "insured_d12",
         "home_computer_d12", "home_internet_d12",
         "marital_d12")

datRes<- replace_with_na_at(datRes,
                            .vars= vari,
                            condition= ~.x== -99)

datRes<- replace_with_na_at(datRes,
                            .vars= c("language_d12", "insured_d12", 
                                     "home_computer_d12", "home_internet_d12"),
                            condition= ~.x== 2)

datRes<- replace_with_na_at(datRes,
                            .vars= c("walk_safety_d12"),
                            condition= ~.x== 4)

datRes<- replace_with_na_at(datRes,
                            .vars= c("home_owner_d12"),
                            condition= ~.x== -88)

table(datRes$home_internet_d12, exclude= F)

# Some of the key variables need to be turned into binary indicator
datRes<- datRes%>%
  dplyr::mutate(
    owner= dplyr::case_when(
      is.na(housing_d12)~ housing_d12,
      housing_d12%in% c(1, 2)~ 1,
      TRUE~ 0
    ),
    distrust_doctor= dplyr::case_when(
      is.na(trust_doctor_d12)~ trust_doctor_d12,
      trust_doctor_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_faith_leader= dplyr::case_when(
      is.na(trust_faith_leader_d12)~ trust_faith_leader_d12,
      trust_faith_leader_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_friends= dplyr::case_when(
      is.na(trust_friends_d12)~ trust_friends_d12,
      trust_friends_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_acquaint= dplyr::case_when(
      is.na(trust_acquaint_d12)~ trust_acquaint_d12,
      trust_acquaint_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_news= dplyr::case_when(
      is.na(trust_news_d12)~ trust_news_d12,
      trust_news_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_socialmedia= dplyr::case_when(
      is.na(trust_socialmedia_d12)~ trust_socialmedia_d12,
      trust_socialmedia_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_usgovt= dplyr::case_when(
      is.na(trust_usgovt_d12)~ trust_usgovt_d12,
      trust_usgovt_d12== 1~ 1,
      TRUE~ 0
    ),
    distrust_taskforce= dplyr::case_when(
      is.na(trust_taskforce_d12)~ trust_taskforce_d12,
      trust_taskforce_d12== 1~ 1,
      TRUE~ 0
    ),
    walk_unsafe= dplyr::case_when(
      is.na(walk_safety_d12)~ walk_safety_d12,
      walk_safety_d12== 1~ 1,
      TRUE~ 0
    ),
    divorce= dplyr::case_when(
      is.na(marital_d12)~ marital_d12,
      marital_d12== 3~ 1,
      TRUE~ 0
    )
    
  )

#+ Split up experimental and control data ====
datRes_exp= datRes[datRes$adapteddesigntreatment!= "Control",]
datRes_con= datRes[datRes$adapteddesigntreatment== "Control",]


# ==== Analysis  ====
# + Demographic distribution comparison ====
#sample sizes of panel and refreshment sample
datRes%>%
  filter(!is.na(gender_age))%>%
  select(Adaptive.Design)%>%
  table()
datRes%>%
  filter(!is.na(educat4))%>%
  select(Adaptive.Design)%>%
  table()
datRes%>%
  filter(!is.na(racecat5))%>%
  select(Adaptive.Design)%>%
  table()
datRes%>%
  filter(!is.na(incomecat5_enrich))%>%
  select(Adaptive.Design)%>%
  table()

#Gender and age distribution 
# --> separate: experiment vs control + new vs established
t1<- table(datRes$gender_age, datRes$experiment_Panel)
prop.table(t1, 2)%>%
  clipr::write_clip(.)
# --> separate: experiment vs control
t2<- table(datRes$gender_age, datRes$experiment)
prop.table(t2, 2)%>%
  clipr::write_clip(.)

#Education distribution 
# --> separate: experiment vs control + new vs established
t1<- table(datRes$educat4, datRes$experiment_Panel)
prop.table(t1, 2)%>%
  clipr::write_clip(.)
# --> separate: experiment vs control
t2<- table(datRes$educat4, datRes$experiment)
prop.table(t2, 2)%>%
  clipr::write_clip(.)


#Race distribution 
# --> separate: experiment vs control + new vs established
t1<- table(datRes$racecat5, datRes$experiment_Panel)
prop.table(t1, 2)%>%
  clipr::write_clip(.)
# --> separate: experiment vs control
t2<- table(datRes$racecat5, datRes$experiment)
prop.table(t2, 2)%>%
  clipr::write_clip(.)


#Income distribution 
# --> separate: experiment vs control + new vs established
t1<- table(datRes$incomecat5_enrich, datRes$experiment_Panel)
prop.table(t1, 2)%>%
  clipr::write_clip(.)
# --> separate: experiment vs control
t2<- table(datRes$incomecat5_enrich, datRes$experiment)
prop.table(t2, 2)%>%
  clipr::write_clip(.)



# ==== Analysis w/ geographic post-stratification (experimental and control separately)====

# Population size by strata (target of poststratification)
postsurvey_strata<- c("Campau / Banglatown", 
                      "East Warren / Cadieux", 
                      "Grand River / Northwest", 
                      "Gratiot / 7-Mile",  
                      "Islandview / Greater Villages", 
                      "Jefferson / Chalmers", 
                      "Livernois / McNichols",
                      "Russell Woods / Nardin Park",
                      "Warrendale / Cody-Rouge",
                      "Southwest / Vernor highHisp", 
                      "Southwest / Vernor lowHisp",
                      "theRest highHisp", 
                      "theRest lowHisp")
strata_pop<- c(7886, 
               17667, 
               29329, 
               13437, 
               10580, 
               5875,
               31935, 
               4401,
               25391,
               11067,
               9106,
               8851,
               331864)
N= sum(strata_pop)
N.ps<- data.frame(postsurvey_strata, strata_pop)

n_exp<- dim(datRes_exp)[1] #sample size
n_con<- dim(datRes_con)[1]

d_exp<- rep(N/n_exp, n_exp) #SRS weights
d_con<- rep(N/n_con, n_con)

f_exp<- rep(n_exp/N, n_exp) #finite population correlation
f_con<- rep(n_con/N, n_con)

#specifying SRS survey design
srs.design_exp<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                           data= datRes_exp, 
                           weight= ~d_exp, fpc= ~f_exp)
srs.design_con<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                           data= datRes_con, 
                           weight= ~d_con, fpc= ~f_con)

#specifying post-stratification survey design
ps.design_exp<- postStratify(design= srs.design_exp,
                             strata= ~postsurvey_strata, 
                             population= N.ps)
ps.design_con<- postStratify(design= srs.design_con,
                             strata= ~postsurvey_strata, 
                             population= N.ps)

#extracting post-stratification weights
weightPs_exp<- weights(ps.design_exp)
weightPs_con<- weights(ps.design_con)


# + Demographic distribution comparison ====
#Gender and age distribution 
# --> separate: experiment vs control
svytable(~gender_age, ps.design_exp)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)
svytable(~gender_age, ps.design_con)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)


#Education distribution 
# --> separate: experiment vs control
svytable(~educat4, ps.design_exp)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)
svytable(~educat4, ps.design_con)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)


#Race distribution 
# --> separate: experiment vs control
svytable(~racecat5, ps.design_exp)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)
svytable(~racecat5, ps.design_con)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)


#Income distribution 
# --> separate: experiment vs control
svytable(~incomecat5_enrich, ps.design_exp)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)
svytable(~incomecat5_enrich, ps.design_con)%>%
  prop.table()%>%
  data.frame()%>%
  clipr::write_clip(.)



# ==== Analysis w/ raking (experimental and control separately)====
# specifying raking distributions
ACS2019_total_pop = 670052
ACS2019_18over_pop = 503934
ACS2019_25over_pop = 441647
ACS2019_hh=267139
#gender_age sum to ACS2019_18over_pop
N.gender_age<- c("gender_age1"= 84118, 
                 "gender_age2"= 73963, 
                 "gender_age3"= 37824, 
                 "gender_age4"= 39890,
                 "gender_age5"= 84975, 
                 "gender_age6"= 83769, 
                 "gender_age7"= 42187, 
                 "gender_age8"= 57208)
P.gender_age<- N.gender_age/ACS2019_18over_pop
#educ sum to ACS2019_25over_pop
N.educ<- c("educat1"= 72271, 
           "educat2"= 154232,
           "educat3"= 141557,
           "educat4"= 73587)


P.educ<- N.educ/ACS2019_25over_pop
#race sum to ACS2019_total_pop
N.race<- c("racecat1"= 70901, 
           "racecat2"= 518305, 
           "racecat3"= 8846, #Important: 3 and 4 reversed after debuging race
           "racecat4"= 16711, 
           "racecat5"= 55289)
P.race<- N.race/ACS2019_total_pop
#income sum to ACS2019_hh
N.income<-  c("incomecat1"= 45045, 
              "incomecat2"= 73065, 
              "incomecat3"= 57748, 
              "incomecat4"= 63124, 
              "incomecat5"= 28157)
P.income<- N.income/ACS2019_hh

#combining the categories into one long vector (required by raking)
pop.P<- c('(Intercept)'= 1, 
          P.gender_age[-1], 
          P.educ[-1],
          P.race[-1],
          P.income[-1])
pop.P_exp<- pop.P*n_exp #scale to the size of the experimental data
pop.P_con<- pop.P*n_con #scale to the size of the control data

#survey design variable - raking
rake.design_exp<- calibrate(design= ps.design_exp, 
                            formula= ~gender_age_impute+ 
                              educat4_impute+
                              racecat5_impute+
                              incomecat5_impute,
                            calfun= "raking", 
                            population= pop.P_exp)

#survey design variable - raking
rake.design_con<- calibrate(design= ps.design_con, 
                            formula= ~gender_age_impute+ 
                              educat4_impute+
                              racecat5_impute+
                              incomecat5_impute,
                            calfun= "raking", 
                            population= pop.P_con)

#trim weights
trim.design_exp<- trimWeights(design= rake.design_exp, lower= 0.2, upper= 6)
trim.design_con<- trimWeights(design= rake.design_con, lower= 0.2, upper= 6)

weights_trim_exp<- weights(trim.design_exp)
weights_trim_con<- weights(trim.design_con)


#====Bootstrapping====
# + bootstrapping experimental sample====
#I decide to save bootstrap samples and perform analysis on them separately

#Wrap bootstrapping in a false if to mute this chunk of code 
if (0==1){
  setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/exp")
  sink("quiet")
  for (i in 1:5000){
    set.seed(i)
    sam_bs<- datRes_exp%>%
      dplyr::group_by(Region)%>%
      sample_frac(1, replace= T) #bootstrap sample stratified by Region
    
    sam_bs$sample<- i 
    
    #Imputing demographic features for each of the bootstrap sample
    if (1==1){
      sam_bs$gender_impute<- sam_bs$gender%>%
        dplyr::recode("man"= 1, "woman"= 0, .missing= NA_real_)
      sam_bs$gender_impute<- factor(sam_bs$gender_impute)
      sam_bs$age_impute<- sam_bs$ageinyears_d12
      sam_bs$racecat5_impute<- factor(sam_bs$racecat5)
      sam_bs$educat6_impute<- factor(sam_bs$educat6)
      sam_bs$incomecat5_impute<- factor(sam_bs$incomecat5)
      sam_bs$income_1_enrich<- factor(sam_bs$income_1_d12)
      
      impute_vari<- c("gender_impute", "age_impute", 
                      "racecat5_impute", "educat6_impute", 
                      "incomecat5_impute")
      help_vari<- c("income_1_enrich",
                    "percent_nonhisp_black", #ACS
                    "percent_hisp", 
                    "percent_lessthan10000",
                    "percent_10kto29k",
                    "percent_30kto49k",
                    "percent_50kto99k")
      #imputation
      tempData <- mice(sam_bs[c(impute_vari, help_vari)], 
                       m=1,
                       maxit=20,
                       meth= c("logreg", "pmm", 
                               "polyreg", "polyreg", "polyreg",
                               "logreg",
                               "pmm", "pmm", "pmm", 
                               "pmm", "pmm", "pmm"),
                       seed=i)

      sam_bs[c(impute_vari, help_vari)]<- complete(tempData)
      
      sam_bs<- sam_bs%>% 
        mutate(agecat4_impute= dplyr::case_when(
          age_impute< 35~ 1,
          age_impute>= 35 & age_impute< 55~ 2,
          age_impute>= 55 & age_impute< 65~ 3,
          age_impute>= 65~ 4
        ))
  
      sam_bs<- sam_bs%>%
        mutate(
          gender_age_impute= dplyr::case_when(
            gender_impute==1 & agecat4_impute== 1 ~ 1, #"man<35",
            gender_impute==1 & agecat4_impute== 2 ~ 2, #"man35-54",
            gender_impute==1 & agecat4_impute== 3 ~ 3, #"man55-64",
            gender_impute==1 & agecat4_impute== 4 ~ 4, #"man>65",
            gender_impute==0 & agecat4_impute== 1 ~ 5, #"woman<35",
            gender_impute==0 & agecat4_impute== 2 ~ 6, #"woman35-54",
            gender_impute==0 & agecat4_impute== 3 ~ 7, #"woman55-64",
            gender_impute==0 & agecat4_impute== 4 ~ 8 #"woman>65",
          ))
      
      sam_bs<- sam_bs%>%
        mutate(educat4_impute= dplyr::case_when(
          educat6_impute==1~ 1,
          educat6_impute==2~ 2,
          educat6_impute%in% c(3, 4)~ 3,
          educat6_impute%in% c(5, 6)~ 4
        ))
    }
    
    name<- paste0("Samples_exp", i, ".csv")
    write.csv(sam_bs, name)
  }
  sink()
}


#+ bootstrapping control samples ====
#Wrap bootstrapping in a false if to mute this chunk of code 
if (0==1){
  setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/con")
  sink("quiet")
  for (i in 1:5000){
    set.seed(i)
    sam_bs<- datRes_con%>%
      dplyr::group_by(Region)%>%
      sample_frac(1, replace= T) #bootstrap sample stratified by Region
    
    sam_bs$sample<- i #numbering the # of iteration
    
    #Imputing demographic features for each of the bootstrap sample
    if (1==1){
      # Re-creating demographic variables These variables will be replaced with imputed values
      sam_bs$gender_impute<- sam_bs$gender%>%
        dplyr::recode("man"= 1, "woman"= 0, .missing= NA_real_)
      sam_bs$gender_impute<- factor(sam_bs$gender_impute)
      sam_bs$age_impute<- sam_bs$ageinyears_d12
      sam_bs$racecat5_impute<- factor(sam_bs$racecat5)
      sam_bs$educat6_impute<- factor(sam_bs$educat6)
      sam_bs$incomecat5_impute<- factor(sam_bs$incomecat5)
      sam_bs$income_1_enrich<- factor(sam_bs$income_1_d12)
      
      impute_vari<- c("gender_impute", "age_impute", 
                      "racecat5_impute", "educat6_impute", 
                      "incomecat5_impute")
      help_vari<- c("income_1_enrich",
                    "percent_nonhisp_black", #ACS
                    "percent_hisp", 
                    "percent_lessthan10000",
                    "percent_10kto29k",
                    "percent_30kto49k",
                    "percent_50kto99k")
      #imputation
      tempData <- mice(sam_bs[c(impute_vari, help_vari)], 
                       m=1,
                       maxit=20,
                       meth= c("logreg", "pmm", 
                               "polyreg", "polyreg", "polyreg",
                               "logreg",
                               "pmm", "pmm", "pmm", 
                               "pmm", "pmm", "pmm"),
                       seed=i)
   
      sam_bs[c(impute_vari, help_vari)]<- complete(tempData)
      
      sam_bs<- sam_bs%>% 
        mutate(agecat4_impute= dplyr::case_when(
          age_impute< 35~ 1,
          age_impute>= 35 & age_impute< 55~ 2,
          age_impute>= 55 & age_impute< 65~ 3,
          age_impute>= 65~ 4
        ))
      
      sam_bs<- sam_bs%>%
        mutate(
          gender_age_impute= dplyr::case_when(
            gender_impute==1 & agecat4_impute== 1 ~ 1, #"man<35",
            gender_impute==1 & agecat4_impute== 2 ~ 2, #"man35-54",
            gender_impute==1 & agecat4_impute== 3 ~ 3, #"man55-64",
            gender_impute==1 & agecat4_impute== 4 ~ 4, #"man>65",
            gender_impute==0 & agecat4_impute== 1 ~ 5, #"woman<35",
            gender_impute==0 & agecat4_impute== 2 ~ 6, #"woman35-54",
            gender_impute==0 & agecat4_impute== 3 ~ 7, #"woman55-64",
            gender_impute==0 & agecat4_impute== 4 ~ 8 #"woman>65",
          ))
      
      sam_bs<- sam_bs%>%
        mutate(educat4_impute= dplyr::case_when(
          educat6_impute==1~ 1,
          educat6_impute==2~ 2,
          educat6_impute%in% c(3, 4)~ 3,
          educat6_impute%in% c(5, 6)~ 4
        ))
    }
    
    name<- paste0("Samples_con", i, ".csv")
    write.csv(sam_bs, name)
  }
  sink()
  
}


#+ Analyzing experimental bootstrap samples====
key_vari<- c("owner", 
             "language_d12",
             "home_computer_d12", "home_internet_d12",
             "divorce", "insured_d12")

means_bs_PS_exp<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_PS_exp)<- key_vari 
SEs_bs_PS_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_PS_exp)<- key_vari 
deffs_bs_PS_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_PS_exp)<- key_vari 

means_bs_RK_exp<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_RK_exp)<- key_vari 
SEs_bs_RK_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_RK_exp)<- key_vari 
deffs_bs_RK_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_RK_exp)<- key_vari 

means_bs_TM_exp<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_TM_exp)<- key_vari 
SEs_bs_TM_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_TM_exp)<- key_vari 
deffs_bs_TM_exp<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_TM_exp)<- key_vari 

n_exp<- nrow(datRes_exp)
setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/exp")
for (i in 1:5000){
  print(i)
  set.seed(i)
  name<- paste0("Samples_exp", i, ".csv")
  sam_bs<- read.csv(name)
  
  sam_bs$gender_age_impute<- as.factor(sam_bs$gender_age_impute)
  sam_bs$educat6_impute<- as.factor(sam_bs$educat6_impute)
  sam_bs$educat4_impute<- as.factor(sam_bs$educat4_impute)
  sam_bs$racecat5_impute<- as.factor(sam_bs$racecat5_impute)
  sam_bs$incomecat5_impute<- as.factor(sam_bs$incomecat5_impute)
  
  d_exp<- rep(N/n_exp, n_exp) 
  f_exp<- rep(n_exp/N, n_exp) 
  
  #specifying SRS survey design
  srs.design_exp<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                             data= sam_bs,
                             weight= ~d_exp, fpc= ~f_exp)
  #specifying post-stratification survey design
  ps.design_exp<- postStratify(design= srs.design_exp,
                               strata= ~postsurvey_strata, 
                               population= N.ps) 
  #specifying raking survey design
  rake.design_exp<- calibrate(design= ps.design_exp, 
                              formula= ~gender_age_impute+ 
                                educat4_impute+
                                racecat5_impute+
                                incomecat5_impute,
                              calfun= "raking", 
                              population= pop.P_exp)
  #trim raking weights
  trim.design_exp<- trimWeights(design= rake.design_exp, lower= 0.2, upper= 6)
  
  for (kv in key_vari){
    #post-stratified
    k_mean<- svymean(~ eval(parse(text= kv)), ps.design_exp, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_PS_exp[i, kv]<- m
    SEs_bs_PS_exp[i, kv]<- se
    deffs_bs_PS_exp[i, kv]<- deff
    
    #raking
    k_mean<- svymean(~ eval(parse(text= kv)), rake.design_exp, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_RK_exp[i, kv]<- m
    SEs_bs_RK_exp[i, kv]<- se
    deffs_bs_RK_exp[i, kv]<- deff
    
    #trimmed weights
    k_mean<- svymean(~ eval(parse(text= kv)), trim.design_exp, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_TM_exp[i, kv]<- m
    SEs_bs_TM_exp[i, kv]<- se
    deffs_bs_TM_exp[i, kv]<- deff
  }
}


# + Analyzing control bootstrap sample====
means_bs_PS_con<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_PS_con)<- key_vari 
SEs_bs_PS_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_PS_con)<- key_vari 
deffs_bs_PS_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_PS_con)<- key_vari 

means_bs_RK_con<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_RK_con)<- key_vari 
SEs_bs_RK_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_RK_con)<- key_vari 
deffs_bs_RK_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_RK_con)<- key_vari 

means_bs_TM_con<- data.frame(matrix(ncol= length(key_vari))) 
colnames(means_bs_TM_con)<- key_vari 
SEs_bs_TM_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(SEs_bs_TM_con)<- key_vari 
deffs_bs_TM_con<- data.frame(matrix(ncol= length(key_vari)))
colnames(deffs_bs_TM_con)<- key_vari 

n_con<- dim(datRes_con)[1]
setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/con")
for (i in 1:5000){
  print(i)
  set.seed(i)
  name<- paste0("Samples_con", i, ".csv")
  sam_bs<- read.csv(name)
  
  sam_bs$gender_age_impute<- as.factor(sam_bs$gender_age_impute)
  sam_bs$educat6_impute<- as.factor(sam_bs$educat6_impute)
  sam_bs$educat4_impute<- as.factor(sam_bs$educat4_impute)
  sam_bs$racecat5_impute<- as.factor(sam_bs$racecat5_impute)
  sam_bs$incomecat5_impute<- as.factor(sam_bs$incomecat5_impute)
  
  d_con<- rep(N/n_con, n_con) 
  f_con<- rep(n_con/N, n_con) #
  
  #specifying SRS survey design
  srs.design_con<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                             data= sam_bs, #bootstrap sample here 
                             weight= ~d_con, fpc= ~f_con)
  #specifying post-stratification survey design
  ps.design_con<- postStratify(design= srs.design_con,
                               strata= ~postsurvey_strata, 
                               population= N.ps) #N.ps is specified in the block above
  #specifying raking survey design
  rake.design_con<- calibrate(design= ps.design_con, 
                              formula= ~gender_age_impute+ 
                                educat4_impute+
                                racecat5_impute+
                                incomecat5_impute,
                              calfun= "raking", 
                              population= pop.P_con)
  #trim raking weights
  trim.design_con<- trimWeights(design= rake.design_con, lower= 0.2, upper= 6)
  
  for (kv in key_vari){
    #post-stratified
    k_mean<- svymean(~ eval(parse(text= kv)), ps.design_con, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_PS_con[i, kv]<- m
    SEs_bs_PS_con[i, kv]<- se
    deffs_bs_PS_con[i, kv]<- deff
    
    #raking
    k_mean<- svymean(~ eval(parse(text= kv)), rake.design_con, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_RK_con[i, kv]<- m
    SEs_bs_RK_con[i, kv]<- se
    deffs_bs_RK_con[i, kv]<- deff
    
    #raking
    k_mean<- svymean(~ eval(parse(text= kv)), trim.design_con, na.rm= T, deff= "replace")
    m<- data.frame(k_mean)[, 1]
    se<- data.frame(k_mean)[, 2]
    deff<- data.frame(k_mean)[, 3]
    
    means_bs_TM_con[i, kv]<- m
    SEs_bs_TM_con[i, kv]<- se
    deffs_bs_TM_con[i, kv]<- deff
  }
}


# + plotting ====
#plotting (post-stratified & raking; 4 distributions per plot)
key_vari<- c("owner", 
             "language_d12",
             "home_computer_d12", "home_internet_d12",
             "divorce", "insured_d12")


kv<- "insured_d12"
if (1==1){
  
  meanline_TM_exp<- mean(means_bs_TM_exp[,kv])
  CI_TM_exp<- quantile(means_bs_TM_exp[,kv], c(0.025, 0.975))
  meanline_TM_con<- mean(means_bs_TM_con[,kv])
  CI_TM_con<- quantile(means_bs_TM_con[,kv], c(0.025, 0.975))
  #Fot setting the x-axis
  x_low<- min(c(means_bs_TM_exp[, kv], means_bs_TM_con[, kv]))
  x_low<- round_any(x_low, 0.01, f= floor)
  x_high<- max(c(means_bs_TM_exp[, kv], means_bs_TM_con[, kv]))
  x_high<- round_any(x_high, 0.01, f= ceiling)
  
  ggplot()+
    geom_histogram(aes(x= eval(parse(text= kv))), 
                   data= means_bs_TM_exp, #results of experimental data with trimmed weights
                   bins = 50, fill= "royalblue4", alpha = 0.4)+
    geom_vline(xintercept= meanline_TM_exp, 
               color= "royalblue4")+
    geom_errorbarh(aes(y= 400, xmin= CI_TM_exp[[1]], xmax= CI_TM_exp[[2]]), 
                   color= "royalblue4", height = 30)+
    
    
    geom_histogram(aes(x= eval(parse(text= kv))), 
                   data= means_bs_TM_con, #results of control data with trimmied weights
                   bins = 50, fill= "darkorange", alpha = 0.4)+
    geom_vline(xintercept= meanline_TM_con, 
               color= "darkorange")+
    geom_errorbarh(aes(y= 200, xmin= CI_TM_con[[1]], xmax= CI_TM_con[[2]]), 
                   color= "darkorange", height = 30)+
    
    annotate(geom= "text",
             x= x_high-0.032, y=450, label="adaptive (experimental)", 
             color= "royalblue4", size= 3.5)+
    annotate(geom= "text",
             x= x_high-0.12, y=250, label="homogeneous (control)", 
             color= "darkorange", size= 3.5)+
    
    geom_vline(xintercept= 1- ref_insured, color= "black")+
    
    scale_x_continuous(limits= c(x_low, x_high), oob = scales::oob_keep)+
    theme(axis.title.x=element_blank(),
          plot.title = element_text(size= 11))
}


var_ratio_TM<-(SEs_bs_TM_exp[,kv])**2/(SEs_bs_TM_con[,kv])**2
sum(var_ratio_TM<3/7)/5000
ggplot()+
  geom_histogram(aes(var_ratio_TM), bins= 50,
                 fill= "black", alpha = 0.5)+
  annotate(geom= "text",
           x= 0.6, y=300, label="ratio of variances", 
           color= "black", size= 3.5)+
  geom_vline(xintercept= 3/7)+
  theme(axis.title.x=element_blank())


deffs_bs_TM_exp[,kv]%>%
  mean()%>%
  round(2)
deffs_bs_TM_con[,kv]%>%
  mean()%>%
  round(2)


# ==== Multivariate analysis====

#+ weighted regression====

formula1<- vac_covid_d12~  factor(gender_impute)+ factor(racecat5_impute)+ 
  factor(educat4_impute)+ factor(incomecat5_impute)+age_impute+
  distrust_doctor+ distrust_usgovt+
  distrust_faith_leader+ distrust_acquaint+ distrust_socialmedia

formula2<- nb_satis_d12~ nb_reputation_d12+ walk_unsafe+ 
  owner+ home_computer_d12

# + Analyzing experimental bootstrap samples====

n_exp<- nrow(datRes_exp)
setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/exp")
for (i in 1:5000){
  print(i)
  set.seed(i)
  name<- paste0("Samples_exp", i, ".csv")
  sam_bs<- read.csv(name)
  
  #turn these variables into factors for raking
  sam_bs$gender_age_impute<- as.factor(sam_bs$gender_age_impute)
  sam_bs$educat6_impute<- as.factor(sam_bs$educat6_impute)
  sam_bs$educat4_impute<- as.factor(sam_bs$educat4_impute)
  sam_bs$racecat5_impute<- as.factor(sam_bs$racecat5_impute)
  sam_bs$incomecat5_impute<- as.factor(sam_bs$incomecat5_impute)
  
  d_exp<- rep(N/n_exp, n_exp) #SRS weights #N is specified in the block above
  f_exp<- rep(n_exp/N, n_exp) #finite population correlation
  
  #specifying SRS survey design
  srs.design_exp<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                             data= sam_bs, #bootstrap sample here 
                             weight= ~d_exp, fpc= ~f_exp)
  #specifying post-stratification survey design
  ps.design_exp<- postStratify(design= srs.design_exp,
                               strata= ~postsurvey_strata, 
                               population= N.ps) #N.ps is specified in the block above
  #specifying raking survey design
  pop.P_exp<- pop.P*n_exp #scale to the size of the experimental data
  rake.design_exp<- calibrate(design= ps.design_exp, 
                              formula= ~gender_age_impute+ 
                                educat4_impute+
                                racecat5_impute+
                                incomecat5_impute,
                              calfun= "raking", 
                              population= pop.P_exp)
  #trim raking weights
  trim.design_exp<- trimWeights(design= rake.design_exp, lower= 0.2, upper= 6)
  
  #weighted regression
  glm_TR<- svyglm(formula1, trim.design_exp) #family= quasibinomial, 
  
  #storing results
  res_TR<- data.frame(summary(glm_TR)$coefficients)
  
  if (i==1){
    #TR coefficients
    coefficient_TR_exp<- data.frame("c1"= res_TR$Estimate)
    rownames(coefficient_TR_exp)<- rownames(res_TR)
    #TR SE
    SE_TR_exp<- data.frame("c1"= res_TR$Std..Error)
    rownames(SE_TR_exp)<- rownames(res_TR)
    #TR t value
    t_TR_exp<- data.frame("c1"= res_TR$t.value)
    rownames(t_TR_exp)<- rownames(res_TR)
    #TR p value
    p_TR_exp<- data.frame("c1"= res_TR$Pr...t..)
    rownames(p_TR_exp)<- rownames(res_TR)
    
  } else{
    cName<- paste0("c", i)
    
    coefficient_TR_exp[cName]<- res_TR$Estimate
    SE_TR_exp[cName]<- res_TR$Std..Error
    t_TR_exp[cName]<- res_TR$t.value
    p_TR_exp[cName]<- res_TR$Pr...t..
  }
  
}



# + Analyzing control bootstrap samples====
n_con<- nrow(datRes_con)
setwd("/Volumes/Samsung_T5/DMACS/BootstrapSamples/con")
for (i in 1:5000){
  print(i)
  set.seed(i)
  name<- paste0("Samples_con", i, ".csv")
  sam_bs<- read.csv(name)
  
  #turn these variables into factors for raking
  sam_bs$gender_age_impute<- as.factor(sam_bs$gender_age_impute)
  sam_bs$educat6_impute<- as.factor(sam_bs$educat6_impute)
  sam_bs$educat4_impute<- as.factor(sam_bs$educat4_impute)
  sam_bs$racecat5_impute<- as.factor(sam_bs$racecat5_impute)
  sam_bs$incomecat5_impute<- as.factor(sam_bs$incomecat5_impute)
  
  d_con<- rep(N/n_con, n_con) #SRS weights #N is specified in the block above
  f_con<- rep(n_con/N, n_con) #finite population correlation
  
  #specifying SRS survey design
  srs.design_con<- svydesign(ids= ~0, strata= ~postsurvey_strata,
                             data= sam_bs, #bootstrap sample here 
                             weight= ~d_con, fpc= ~f_con)
  #specifying post-stratification survey design
  ps.design_con<- postStratify(design= srs.design_con,
                               strata= ~postsurvey_strata, 
                               population= N.ps) #N.ps is specified in the block above
  #specifying raking survey design
  pop.P_con<- pop.P*n_con #scale to the size of the control data
  rake.design_con<- calibrate(design= ps.design_con, 
                              formula= ~gender_age_impute+ 
                                educat4_impute+
                                racecat5_impute+
                                incomecat5_impute,
                              calfun= "raking", 
                              population= pop.P_con)
  #trim raking weights
  trim.design_con<- trimWeights(design= rake.design_con, lower= 0.2, upper= 6)
  
  #weighted regression
  glm_TR<- svyglm(formula1,  trim.design_con) #family= quasibinomial,  
  
  #storing results
  res_TR<- data.frame(summary(glm_TR)$coefficients)
  
  if (i==1){
    #TR coefficients
    coefficient_TR_con<- data.frame("c1"= res_TR$Estimate)
    rownames(coefficient_TR_con)<- rownames(res_TR)
    #TR SE
    SE_TR_con<- data.frame("c1"= res_TR$Std..Error)
    rownames(SE_TR_con)<- rownames(res_TR)
    #TR t value
    t_TR_con<- data.frame("c1"= res_TR$t.value)
    rownames(t_TR_con)<- rownames(res_TR)
    #TR p value
    p_TR_con<- data.frame("c1"= res_TR$Pr...t..)
    rownames(p_TR_con)<- rownames(res_TR)
    
  } else{
    cName<- paste0("c", i)

    coefficient_TR_con[cName]<- res_TR$Estimate
    SE_TR_con[cName]<- res_TR$Std..Error
    t_TR_con[cName]<- res_TR$t.value
    p_TR_con[cName]<- res_TR$Pr...t..
  }
  
}


# + Storing results temporarily====
setwd("/Volumes/Samsung_T5/DMACS/temporaryResults/vac_covid_d12")
# write.csv(coefficient_TR_exp, "exp/coefficient_TR_exp.csv")
# write.csv(SE_TR_exp, "exp/SE_TR_exp.csv")
# write.csv(t_TR_exp, "exp/t_TR_exp.csv")
# write.csv(p_TR_exp, "exp/p_TR_exp.csv")

# write.csv(coefficient_TR_con, "con/coefficient_TR_con.csv")
# write.csv(SE_TR_con, "con/SE_TR_con.csv")
# write.csv(t_TR_con, "con/t_TR_con.csv")
# write.csv(p_TR_con, "con/p_TR_con.csv")



coefficient_TR_exp<- read.csv("exp/coefficient_TR_exp.csv")
SE_TR_exp<- read.csv("exp/SE_TR_exp.csv")
t_TR_exp<- read.csv("exp/t_TR_exp.csv")
p_TR_exp<- read.csv("exp/p_TR_exp.csv")


coefficient_TR_con<- read.csv( "con/coefficient_TR_con.csv")
SE_TR_con<- read.csv("con/SE_TR_con.csv")
t_TR_con<- read.csv("con/t_TR_con.csv")
p_TR_con<- read.csv("con/p_TR_con.csv")


rownames(coefficient_TR_exp)<- coefficient_TR_exp$X
rownames(SE_TR_exp)<- SE_TR_exp$X
rownames(t_TR_exp)<- t_TR_exp$X
rownames(p_TR_exp)<- p_TR_exp$X

rownames(coefficient_TR_con)<- coefficient_TR_con$X
rownames(SE_TR_con)<- SE_TR_con$X
rownames(t_TR_con)<- t_TR_con$X
rownames(p_TR_con)<- p_TR_con$X


coefficient_TR_exp$X<- NULL
SE_TR_exp$X<- NULL
t_TR_exp$X<- NULL
p_TR_exp$X<- NULL

coefficient_TR_con$X<- NULL
SE_TR_con$X<- NULL
t_TR_con$X<- NULL
p_TR_con$X<- NULL




# + Mean coefficients ====
a<- coefficient_PS_exp%>%
  rowMeans()
data.frame(a)%>%
  write_clip()

a<- coefficient_TR_exp%>%
  rowMeans()
data.frame(a)%>%
  write_clip()

a<- coefficient_PS_con%>%
  rowMeans(na.rm= T)
data.frame(a)%>%
  write_clip()

a<- coefficient_TR_con%>%
  rowMeans(na.rm= T)
data.frame(a)%>%
  write_clip()


# + Adjusting SE and significance====

#The experimental group 
a<- data.frame(p_PS_exp<0.05)%>%
  rowSums()/5000
sigPct_PS_exp<- data.frame(a) #% of significance
write_clip(sigPct_PS_exp) 


a<- data.frame(p_TR_exp<0.05)%>%
  rowSums()/5000
sigPct_TR_exp<- data.frame(a) #% of significance
write_clip(sigPct_TR_exp)

a<- SE_TR_exp%>%
  rowMeans() #Ageraged SE of the experimental group
data.frame(a)%>%
  write_clip() 

#The control group has a smaller sample size --> adjust
sampSizeRatio<- 7/3
sampSizeRatio_sqrt<- sqrt(sampSizeRatio)
df0= 609 

coefficient_PS_con #post-stratified control group
SE_PS_con_Adjusted<- SE_PS_con/sampSizeRatio_sqrt
t_PS_con_Adjusted<- coefficient_PS_con/SE_PS_con_Adjusted
p_PS_con_Adjusted<- lapply(t_PS_con_Adjusted, function(x) 2*pt(-abs(x), df= df0))%>% 
  data.frame() #adjusted p value

a<- data.frame(p_PS_con_Adjusted<0.05)%>%
  rowSums()/5000
sigPct_PS_con<- data.frame(a) #% of significance
write_clip(sigPct_PS_con)

coefficient_TR_con #raked control group
SE_TR_con_Adjusted<- SE_TR_con/sampSizeRatio_sqrt
t_TR_con_Adjusted<- coefficient_TR_con/SE_TR_con_Adjusted
p_TR_con_Adjusted<- lapply(t_TR_con_Adjusted, function(x) 2*pt(-abs(x), df= df0))%>% 
  data.frame()

a<- data.frame(p_TR_con_Adjusted<0.05)%>% 
  rowSums()/5000
sigPct_TR_con<- data.frame(a) #% of significance
write_clip(sigPct_TR_con)

a<- SE_TR_con_Adjusted%>%
  rowMeans() #Ageraged SE of the control group
data.frame(a)%>%
  write_clip() 


#====Revision 1: Mode (unweighted)====
table(datRes$surv_mode_d12, exclude= F)
datRes<- datRes%>%
  dplyr::mutate(telephone= case_when(
    surv_mode_d12== 3~ 1,
    TRUE~ 0
  ))
table(datRes$surv_mode_d12, datRes$telephone, exclude= F)

## ====+ Telephone/web vs. experimental/control?====
t<- table(datRes$telephone, datRes$experiment)
t
t%>% prop.table(1)
chisq.test(t)

## ====+ Telephone/web on demographics?====
t<- table(datRes$telephone, datRes$gender_enrich)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

t<- table(datRes$telephone, datRes$agecat4)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

t<- table(datRes$telephone, datRes$gender_age)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

t<- table(datRes$telephone, datRes$educat4)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

t<- table(datRes$telephone, datRes$racecat5)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

t<- table(datRes$telephone, datRes$incomecat5_enrich)
t%>% t()%>% write_clip()
t%>% prop.table(2)%>% t()%>% write_clip()
chisq.test(t)

## ====+ Telephone/web on substantive variables?====

binary<- c("owner", 
           "language_d12",
           "home_computer_d12", "home_internet_d12",
           "divorce", "insured_d12",
           "inRcd_full",
           "distrust_doctor", "distrust_usgovt", "distrust_faith_leader", 
           "distrust_acquaint", "distrust_socialmedia", 
           "walk_unsafe")

continuous<- c("vac_covid_d12", "nb_satis_d12", "nb_reputation_d12")

# binary variables
t<- datRes%>% dplyr::select(owner, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(language_d12, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(home_computer_d12, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(home_internet_d12, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(divorce, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(insured_d12, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes_enrich%>% dplyr::select(inRcd_full, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(distrust_doctor, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(distrust_usgovt, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(distrust_faith_leader, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(distrust_acquaint, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(distrust_socialmedia, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)

t<- datRes%>% dplyr::select(walk_unsafe, telephone)%>% table()
t%>% prop.table(2)%>% write_clip()
chisq.test(t)


t.test(vac_covid_d12~ telephone, data= datRes_enrich)
t.test(nb_satis_d12~ telephone, data= datRes_enrich)
t.test(nb_reputation_d12~ telephone, data= datRes_enrich)

table(datRes_enrich$vac_covid_d12, exclude= F)
table(datRes_enrich$nb_satis_d12, exclude= F)
table(datRes_enrich$nb_reputation_d12, exclude= F)

