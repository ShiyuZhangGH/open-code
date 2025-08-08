library(MASS)
library(dplyr)
library(mice)
library(ggplot2)
library(ggridges)
library(clipr)



rm(list=ls())


z<- qt(0.975, 999)
z

z_MI<- qt(0.975, 19) # in MI, use m-1 as the degree of freedom
z_MI

set.seed(9)



table_T<- data.frame() # storing coefficient of the Treatment
table_X1<- data.frame() # storing coefficient of X1 (without missing)
table_X2<- data.frame() # storing coefficient of X2 (with missing)

table_T_CIcoverage<- data.frame() #storing the CI coverage (only T)


for (i in 1:5000){
  
  
  ##### DATA GENERATION #####
  n<- 1000 
  
  # X: covariates
  mu<- c(0, 0)
  covariance<- matrix(c(1, 0.3,
                        0.3, 1),
                      2, 2)
  X<- mvrnorm(n, mu, covariance)
  colnames(X)<- c("X1", "X2")
  dat<- data.frame(X)
  
  
  # Treatment
  T<- rbinom(n= n, size= 1, prob= 0.5) 
  dat<- cbind(dat, data.frame("T"= T))
  
  
  
  # Missing mechanisms
  # X1 has no missing; X2 contains missingness
  
  ## MAR ----
  M2_logit<- 0+ 3* dat$X1 
  dat$M2_MAR_prob<- exp(M2_logit)/(1+ exp(M2_logit))
  dat$M2_MAR<- as.numeric(dat$M2_MAR_prob>= runif(n)) # missing indicator under MAR
  

  dat<- dat%>% mutate(
    X2_observed_MAR= case_when(
      M2_MAR== 1~ NA_real_,
      M2_MAR== 0~ X2
    )
  )
  
  # Potential outcome
  # Y(1)
  dat$Y1<- 1+ dat$X1+ dat$X2+ (2+2*dat$X1+ 2*dat$X2)*1+ rnorm(n, mean= 0, sd= 1)
  # Y(0)
  dat$Y0<- 1+ dat$X1+ dat$X2+ (2+2*dat$X1+ 2*dat$X2)*0+ rnorm(n, mean= 0, sd= 1)
  # The experimental outcome, either treated or not treated
  dat<- dat%>% mutate(
    Y_experiment= case_when(
      T==1~ Y1,
      T==0~ Y0
    )
  )
  
  
  ##### DATASET PREPARATION #####
  
  # for single imputation ----
  dat<- dat%>% mutate(
    X2_imputeMean_MAR= case_when(
      is.na(X2_observed_MAR)~ mean(dat$X2_observed_MAR, na.rm= TRUE),
      TRUE~ X2_observed_MAR
    )
  )
  
  # for the missing indicator method ----
  dat<- dat%>% mutate(
    X2_impute0_MAR= case_when(
      is.na(X2_observed_MAR)~ 0,
      TRUE~ X2_observed_MAR
    )
  )
  
  # for multiple imputation with baseline variables ----
  imputationVar<- c("X1", "X2_observed_MAR") 
  
  imp.dat<- mice(data= dat[, imputationVar],
                 m= 20, maxit= 15, print= FALSE)
  

  dat_rest<- dat[, !colnames(dat)%in% imputationVar]
  
  dat_multiple_MAR_designStage<- complete(imp.dat, "long")%>% 
    group_by(.imp)%>% 
    mutate(cbind(dat_rest))%>% 
    ungroup()
  
  
  # for multiple imputation overall ----
  imputationVar<- c("X1", "X2_observed_MAR", "Y_experiment", "T") 
  
  imp.dat<- mice(data= dat[, imputationVar],
                 m= 20, maxit= 15, print= FALSE)
  
  dat_rest<- dat[, !colnames(dat)%in% imputationVar]
  
  dat_multiple_MAR_overall<- complete(imp.dat, "long")%>% 
    group_by(.imp)%>% 
    mutate(cbind(dat_rest))%>% 
    ungroup()
  
  
  # for multiple imputation by arm ----
  imputationVar<- c("X1", "X2_observed_MAR", "Y_experiment") 
  
  imp.dat_T1<- mice(data= dat[dat$T==1, imputationVar], #treated
                    m= 20, maxit= 15, print= FALSE)
  
  imp.dat_T0<- mice(data= dat[dat$T==0, imputationVar], #not treated
                    m= 20, maxit= 15, print= FALSE)
  
  dat_rest_T1<- dat[dat$T==1, !colnames(dat)%in% imputationVar]
  dat_rest_T0<- dat[dat$T==0, !colnames(dat)%in% imputationVar]
  
  aa<- complete(imp.dat_T1, "long")%>% 
    group_by(.imp)%>% 
    mutate(cbind(dat_rest_T1))%>% 
    ungroup()
  
  bb<- complete(imp.dat_T0, "long")%>% 
    group_by(.imp)%>% 
    mutate(cbind(dat_rest_T0))%>% 
    ungroup()
  
  dat_multiple_MAR_bygroup<- rbind(aa, bb)
  
  
  ##### Modeling ##### 
  
  # 0. no missing ----
  # centering
  dat<- dat%>% mutate(
    X1_predi= X1- mean(X1),
    X2_predi= X2- mean(X2)
  )
  
  # model 1
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "noMiss", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "noMiss", "mainEffect", "X1", summary(lm)$coefficients["X1_predi", "Estimate"], summary(lm)$coefficients["X1_predi", "Std. Error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "noMiss", "mainEffect", "X2", summary(lm)$coefficients["X2_predi", "Estimate"], summary(lm)$coefficients["X2_predi", "Std. Error"])
  table_X2<- rbind(table_X2, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "noMiss",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  
  
  # 1. CCA ----

  # centering
  dat<- dat%>% mutate(
    X1_predi= X1- mean(X1),
    X2_predi= X2_observed_MAR- mean(X2_observed_MAR, na.rm= TRUE) 
  )
  
  # model 1
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "completeCase_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "completeCase_MAR", "mainEffect", "X1", summary(lm)$coefficients["X1_predi", "Estimate"], summary(lm)$coefficients["X1_predi", "Std. Error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "completeCase_MAR", "mainEffect", "X2", summary(lm)$coefficients["X2_predi", "Estimate"], summary(lm)$coefficients["X2_predi", "Std. Error"])
  table_X2<- rbind(table_X2, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "completeCase_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)

  
  
  # 2. single imputation (impute mean)----
  
  # centering
  dat<- dat%>% mutate(
    X1_predi= X1- mean(X1), 
    X2_predi= X2_imputeMean_MAR- mean(X2_imputeMean_MAR, na.rm= TRUE), 
  )
  
  # model 1
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "singleImputeMean_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "singleImputeMean_MAR", "mainEffect", "X1", summary(lm)$coefficients["X1_predi", "Estimate"], summary(lm)$coefficients["X1_predi", "Std. Error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "singleImputeMean_MAR", "mainEffect", "X2", summary(lm)$coefficients["X2_predi", "Estimate"], summary(lm)$coefficients["X2_predi", "Std. Error"])
  table_X2<- rbind(table_X2, row)
  
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "singleImputeMean_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  

  
  
  # 3. missing indicator ----

  # centering
  dat<- dat%>% mutate(
    X1_predi= X1- mean(X1), 
    X2_predi= X2_impute0_MAR- mean(X2_impute0_MAR, na.rm= TRUE), 
    M2_predi= M2_MAR- mean(M2_MAR)
  )
  
  # model 1
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + M2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "missIndicator_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "missIndicator_MAR", "mainEffect", "X1", summary(lm)$coefficients["X1_predi", "Estimate"], summary(lm)$coefficients["X1_predi", "Std. Error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "missIndicator_MAR", "mainEffect", "X2", summary(lm)$coefficients["X2_predi", "Estimate"], summary(lm)$coefficients["X2_predi", "Std. Error"])
  table_X2<- rbind(table_X2, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "missIndicator_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  
  
  # 5. MI (baseline only/design-stage) ----
  
  # model 1
  lm<- dat_multiple_MAR_designStage%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_designStage_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "MI_designStage_MAR", "mainEffect", "X1", summary(lm)[summary(lm)$term== "X1_predi", "estimate"], summary(lm)[summary(lm)$term== "X1_predi", "std.error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "MI_designStage_MAR", "mainEffect", "X2", summary(lm)[summary(lm)$term== "X2_predi", "estimate"], summary(lm)[summary(lm)$term== "X2_predi", "std.error"])
  table_X2<- rbind(table_X2, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_designStage_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  
  
  # 5. MI (overall) ----
  
  # model 1
  lm<- dat_multiple_MAR_overall%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_overall_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "MI_overall_MAR", "mainEffect", "X1", summary(lm)[summary(lm)$term== "X1_predi", "estimate"], summary(lm)[summary(lm)$term== "X1_predi", "std.error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "MI_overall_MAR", "mainEffect", "X2", summary(lm)[summary(lm)$term== "X2_predi", "estimate"], summary(lm)[summary(lm)$term== "X2_predi", "std.error"])
  table_X2<- rbind(table_X2, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_overall_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  
  
  
  # 5. MI (by arm) ----
  
  # model 1
  lm<- dat_multiple_MAR_bygroup%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_bygroup_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T<- rbind(table_T, row)
  
  row<- c(i, "MI_bygroup_MAR", "mainEffect", "X1", summary(lm)[summary(lm)$term== "X1_predi", "estimate"], summary(lm)[summary(lm)$term== "X1_predi", "std.error"])
  table_X1<- rbind(table_X1, row)
  
  row<- c(i, "MI_bygroup_MAR", "mainEffect", "X2", summary(lm)[summary(lm)$term== "X2_predi", "estimate"], summary(lm)[summary(lm)$term== "X2_predi", "std.error"])
  table_X2<- rbind(table_X2, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_bygroup_MAR",  T_CIcoverage)
  table_T_CIcoverage<- rbind(table_T_CIcoverage, row)
  
  
  
}



##### RESULT TABLE #####

mse <- function(estimates, estimand){
  n<- length(estimates)
  sum((estimates- estimand)^2)/n
}



## Treatment effect ----
# table formatting
colnames(table_T)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_T$estimate<- as.numeric(table_T$estimate)
table_T$SE<- as.numeric(table_T$SE)

table_T$missingMethod_Mechanism<- factor(table_T$missingMethod_Mechanism, 
                                         levels= c("noMiss", 
                                                   "completeCase_MAR",
                                                   "singleImputeMean_MAR", 
                                                   "missIndicator_MAR", 
                                                   "MI_designStage_MAR", 
                                                   "MI_overall_MAR", 
                                                   "MI_bygroup_MAR", 
                                                   
                                                   "completeCase_MNAR", 
                                                   "singleImputeMean_MNAR",
                                                   "missIndicator_MNAR",
                                                   "MI_designStage_MNAR",
                                                   "MI_overall_MNAR",
                                                   "MI_bygroup_MNAR"))

# results 
table_T%>%
  group_by(missingMethod_Mechanism, model)%>%
  summarise(mean(estimate), sd(estimate), mean(SE), mse(estimate, 2))%>%
  data.frame()

## CI coverage of T ----

colnames(table_T_CIcoverage)<- c("rep", "missingMethod_Mechanism", "yes_no_coverage")
table_T_CIcoverage$yes_no_coverage<- as.numeric(table_T_CIcoverage$yes_no_coverage)

table_T_CIcoverage$missingMethod_Mechanism<- factor(table_T_CIcoverage$missingMethod_Mechanism, 
                                         levels= c("noMiss", 
                                                   "completeCase_MAR",
                                                   "singleImputeMean_MAR", 
                                                   "missIndicator_MAR", 
                                                   "MI_designStage_MAR", 
                                                   "MI_overall_MAR", 
                                                   "MI_bygroup_MAR", 
                                                   
                                                   "completeCase_MNAR", 
                                                   "singleImputeMean_MNAR",
                                                   "missIndicator_MNAR",
                                                   "MI_designStage_MNAR",
                                                   "MI_overall_MNAR",
                                                   "MI_bygroup_MNAR"))


table_T_CIcoverage%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(yes_no_coverage))%>%
  data.frame()




## Coefficient of X1 ----

colnames(table_X1)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_X1$estimate<- as.numeric(table_X1$estimate)
table_X1$SE<- as.numeric(table_X1$SE)

table_X1$missingMethod_Mechanism<- factor(table_X1$missingMethod_Mechanism, 
                                          levels= c("noMiss", 
                                                    "completeCase_MAR",
                                                    "singleImputeMean_MAR", 
                                                    "missIndicator_MAR", 
                                                    "MI_designStage_MAR", 
                                                    "MI_overall_MAR", 
                                                    "MI_bygroup_MAR", 
                                                    
                                                    "completeCase_MNAR", 
                                                    "singleImputeMean_MNAR",
                                                    "missIndicator_MNAR",
                                                    "MI_designStage_MNAR",
                                                    "MI_overall_MNAR",
                                                    "MI_bygroup_MNAR"))

# results 
table_X1%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(estimate), sd(estimate), mean(SE))%>%
  data.frame()



## Coefficient of X2 ----

colnames(table_X2)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_X2$estimate<- as.numeric(table_X2$estimate)
table_X2$SE<- as.numeric(table_X2$SE)

table_X2$missingMethod_Mechanism<- factor(table_X2$missingMethod_Mechanism, 
                                          levels= c("noMiss", 
                                                    "completeCase_MAR",
                                                    "singleImputeMean_MAR", 
                                                    "missIndicator_MAR", 
                                                    "MI_designStage_MAR", 
                                                    "MI_overall_MAR", 
                                                    "MI_bygroup_MAR", 
                                                    
                                                    "completeCase_MNAR", 
                                                    "singleImputeMean_MNAR",
                                                    "missIndicator_MNAR",
                                                    "MI_designStage_MNAR",
                                                    "MI_overall_MNAR",
                                                    "MI_bygroup_MNAR"))

# results 
table_X2%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(estimate), sd(estimate), mean(SE))%>%
  data.frame()



