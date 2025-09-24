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



table_T_primary<- data.frame() # storing coefficient of the Treatment
table_T_CIcoverage_primary<- data.frame() #storing the CI coverage (only T)



table_T_secondary<- data.frame() # storing coefficient of the Treatment
table_X1T_secondary<- data.frame() # storing coefficient of the interaction term
table_X2T_secondary<- data.frame() # storing coefficient of the interaction term

table_T_CIcoverage_secondary<- data.frame() #storing the CI coverage (T)
table_X1T_CIcoverage_secondary<- data.frame() #storing the CI coverage (X1T)
table_X2T_CIcoverage_secondary<- data.frame() #storing the CI coverage (X2T)


for (i in 1:5000){
  
  print(i)
  
  
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
  M2_logit<- 0+ 3* dat$X1 # -1+ 3* dat$X1
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
  
  # 00. unadjusted ----
  
  lm<- lm(Y_experiment~ T, data= dat)
  summary(lm)
  row<- c(i, "unadjusted", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T_primary<- rbind(table_T_primary, row)
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "unadjusted",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
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
  table_T_primary<- rbind(table_T_primary, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "noMiss",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  # model 2
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T*X1_predi+ T*X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "noMiss", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  
  row<- c(i, "noMiss", "interaction", "X1T", summary(lm)$coefficients["X1_predi:T", "Estimate"], summary(lm)$coefficients["X1_predi:T", "Std. Error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "noMiss", "interaction", "X2T", summary(lm)$coefficients["X2_predi:T", "Estimate"], summary(lm)$coefficients["X2_predi:T", "Std. Error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "noMiss",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)$coefficients["X1_predi:T", "Estimate"]
  X1T_SE<- summary(lm)$coefficients["X1_predi:T", "Std. Error"]
  X1T_CI_low<- X1T_estimate- z*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "noMiss",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)$coefficients["X2_predi:T", "Estimate"]
  X2T_SE<- summary(lm)$coefficients["X2_predi:T", "Std. Error"]
  X2T_CI_low<- X2T_estimate- z*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "noMiss",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  
  
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
  table_T_primary<- rbind(table_T_primary, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "completeCase_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  # model 2
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T*X1_predi+ T*X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "completeCase_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  row<- c(i, "completeCase_MAR", "interaction", "X1T", summary(lm)$coefficients["X1_predi:T", "Estimate"], summary(lm)$coefficients["X1_predi:T", "Std. Error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "completeCase_MAR", "interaction", "X2T", summary(lm)$coefficients["X2_predi:T", "Estimate"], summary(lm)$coefficients["X2_predi:T", "Std. Error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "completeCase_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)$coefficients["X1_predi:T", "Estimate"]
  X1T_SE<- summary(lm)$coefficients["X1_predi:T", "Std. Error"]
  X1T_CI_low<- X1T_estimate- z*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "completeCase_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)$coefficients["X2_predi:T", "Estimate"]
  X2T_SE<- summary(lm)$coefficients["X2_predi:T", "Std. Error"]
  X2T_CI_low<- X2T_estimate- z*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "completeCase_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  

  
  
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
  table_T_primary<- rbind(table_T_primary, row)

  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "singleImputeMean_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
  # model 2
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + T*X1_predi+ T*X2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "singleImputeMean_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  
  row<- c(i, "singleImputeMean_MAR", "interaction", "X1T", summary(lm)$coefficients["X1_predi:T", "Estimate"], summary(lm)$coefficients["X1_predi:T", "Std. Error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "singleImputeMean_MAR", "interaction", "X2T", summary(lm)$coefficients["X2_predi:T", "Estimate"], summary(lm)$coefficients["X2_predi:T", "Std. Error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "singleImputeMean_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)$coefficients["X1_predi:T", "Estimate"]
  X1T_SE<- summary(lm)$coefficients["X1_predi:T", "Std. Error"]
  X1T_CI_low<- X1T_estimate- z*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "singleImputeMean_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)$coefficients["X2_predi:T", "Estimate"]
  X2T_SE<- summary(lm)$coefficients["X2_predi:T", "Std. Error"]
  X2T_CI_low<- X2T_estimate- z*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "singleImputeMean_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  

  
  
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
  table_T_primary<- rbind(table_T_primary, row)
  

  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "missIndicator_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
  # model 2
  lm<- lm(Y_experiment~ X1_predi+ X2_predi
          + M2_predi
          + T*X1_predi+ + T*X2_predi+ T*M2_predi
          + T, data= dat)
  summary(lm)
  
  row<- c(i, "missIndicator_MAR", "mainEffect", "T", summary(lm)$coefficients["T", "Estimate"], summary(lm)$coefficients["T", "Std. Error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  row<- c(i, "missIndicator_MAR", "interaction", "X1T", summary(lm)$coefficients["X1_predi:T", "Estimate"], summary(lm)$coefficients["X1_predi:T", "Std. Error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "missIndicator_MAR", "interaction", "X2T", summary(lm)$coefficients["X2_predi:T", "Estimate"], summary(lm)$coefficients["X2_predi:T", "Std. Error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  
  T_estimate<- summary(lm)$coefficients["T", "Estimate"]
  T_SE<- summary(lm)$coefficients["T", "Std. Error"]
  T_CI_low<- T_estimate- z*T_SE
  T_CI_hig<- T_estimate+ z*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "missIndicator_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)$coefficients["X1_predi:T", "Estimate"]
  X1T_SE<- summary(lm)$coefficients["X1_predi:T", "Std. Error"]
  X1T_CI_low<- X1T_estimate- z*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "missIndicator_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)$coefficients["X2_predi:T", "Estimate"]
  X2T_SE<- summary(lm)$coefficients["X2_predi:T", "Std. Error"]
  X2T_CI_low<- X2T_estimate- z*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "missIndicator_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  
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
  table_T_primary<- rbind(table_T_primary, row)
  

  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_designStage_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
  # model 2
  lm<- dat_multiple_MAR_designStage%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T*X1_predi+ T*X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_designStage_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  row<- c(i, "MI_designStage_MAR", "interaction", "X1T", summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "MI_designStage_MAR", "interaction", "X2T", summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_designStage_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"]
  X1T_SE<- summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"]
  X1T_CI_low<- X1T_estimate- z_MI*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z_MI*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "MI_designStage_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"]
  X2T_SE<- summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"]
  X2T_CI_low<- X2T_estimate- z_MI*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z_MI*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "MI_designStage_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  
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
  table_T_primary<- rbind(table_T_primary, row)

  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_overall_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
  
  # model 2
  lm<- dat_multiple_MAR_overall%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T*X1_predi+ T*X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_overall_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  
  row<- c(i, "MI_overall_MAR", "interaction", "X1T", summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "MI_overall_MAR", "interaction", "X2T", summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_overall_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"]
  X1T_SE<- summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"]
  X1T_CI_low<- X1T_estimate- z_MI*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z_MI*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "MI_overall_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"]
  X2T_SE<- summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"]
  X2T_CI_low<- X2T_estimate- z_MI*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z_MI*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "MI_overall_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  
  
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
  table_T_primary<- rbind(table_T_primary, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_bygroup_MAR",  T_CIcoverage)
  table_T_CIcoverage_primary<- rbind(table_T_CIcoverage_primary, row)
  
  
  
  
  # model 2
  lm<- dat_multiple_MAR_bygroup%>%
    group_by(.imp)%>%
    mutate(X1_predi= X1- mean(X1),
           X2_predi= X2_observed_MAR- mean(X2_observed_MAR))%>%
    do(model = lm(Y_experiment~ X1_predi+ X2_predi
                  + T*X1_predi+ T*X2_predi
                  + T,
                  data = .)) %>%
    as.list() %>%
    .[[-1]] %>%
    pool()
  
  row<- c(i, "MI_bygroup_MAR", "mainEffect", "T", summary(lm)[summary(lm)$term== "T", "estimate"], summary(lm)[summary(lm)$term== "T", "std.error"])
  table_T_secondary<- rbind(table_T_secondary, row)
  
  row<- c(i, "MI_bygroup_MAR", "interaction", "X1T", summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"])
  table_X1T_secondary<- rbind(table_X1T_secondary, row)
  
  row<- c(i, "MI_bygroup_MAR", "interaction", "X2T", summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"], summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"])
  table_X2T_secondary<- rbind(table_X2T_secondary, row)
  
  
  
  T_estimate<- summary(lm)[summary(lm)$term== "T", "estimate"]
  T_SE<- summary(lm)[summary(lm)$term== "T", "std.error"]
  T_CI_low<- T_estimate- z_MI*T_SE
  T_CI_hig<- T_estimate+ z_MI*T_SE
  T_CIcoverage<- as.numeric(T_CI_low<= 2 & T_CI_hig>= 2)
  row<- c(i, "MI_bygroup_MAR",  T_CIcoverage)
  table_T_CIcoverage_secondary<- rbind(table_T_CIcoverage_secondary, row)
  
  
  X1T_estimate<- summary(lm)[summary(lm)$term== "X1_predi:T", "estimate"]
  X1T_SE<- summary(lm)[summary(lm)$term== "X1_predi:T", "std.error"]
  X1T_CI_low<- X1T_estimate- z_MI*X1T_SE
  X1T_CI_hig<- X1T_estimate+ z_MI*X1T_SE
  X1T_CIcoverage<- as.numeric(X1T_CI_low<= 2 & X1T_CI_hig>= 2)
  row<- c(i, "MI_bygroup_MAR",  X1T_CIcoverage)
  table_X1T_CIcoverage_secondary<- rbind(table_X1T_CIcoverage_secondary, row)
  
  
  X2T_estimate<- summary(lm)[summary(lm)$term== "X2_predi:T", "estimate"]
  X2T_SE<- summary(lm)[summary(lm)$term== "X2_predi:T", "std.error"]
  X2T_CI_low<- X2T_estimate- z_MI*X2T_SE
  X2T_CI_hig<- X2T_estimate+ z_MI*X2T_SE
  X2T_CIcoverage<- as.numeric(X2T_CI_low<= 2 & X2T_CI_hig>= 2)
  row<- c(i, "MI_bygroup_MAR",  X2T_CIcoverage)
  table_X2T_CIcoverage_secondary<- rbind(table_X2T_CIcoverage_secondary, row)
  
  
  
  
}



##### RESULT TABLE #####

mse <- function(estimates, estimand){
  n<- length(estimates)
  sum((estimates- estimand)^2)/n
}




### Model 1 ----

## Treatment effect ----
# table formatting
colnames(table_T_primary)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_T_primary$estimate<- as.numeric(table_T_primary$estimate)
table_T_primary$SE<- as.numeric(table_T_primary$SE)

table_T_primary$missingMethod_Mechanism<- factor(table_T_primary$missingMethod_Mechanism, 
                                         levels= c("unadjusted",
                                                   "noMiss", 
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
table_T_primary%>%
  group_by(missingMethod_Mechanism, model)%>%
  summarise(mean(estimate), sd(estimate), mean(SE), mse(estimate, 2))%>%
  data.frame()%>% write_clip()

## CI coverage of T ----

colnames(table_T_CIcoverage_primary)<- c("rep", "missingMethod_Mechanism", "yes_no_coverage")
table_T_CIcoverage_primary$yes_no_coverage<- as.numeric(table_T_CIcoverage_primary$yes_no_coverage)

table_T_CIcoverage_primary$missingMethod_Mechanism<- factor(table_T_CIcoverage_primary$missingMethod_Mechanism, 
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


table_T_CIcoverage_primary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(yes_no_coverage))%>%
  data.frame()%>% write_clip()




### Model 2 ----

## Treatment effect ----
# table formatting
colnames(table_T_secondary)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_T_secondary$estimate<- as.numeric(table_T_secondary$estimate)
table_T_secondary$SE<- as.numeric(table_T_secondary$SE)

table_T_secondary$missingMethod_Mechanism<- factor(table_T_secondary$missingMethod_Mechanism, 
                                                   levels= c("unadjusted",
                                                             "noMiss", 
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
table_T_secondary%>%
  group_by(missingMethod_Mechanism, model)%>%
  summarise(mean(estimate), sd(estimate), mean(SE), mse(estimate, 2))%>%
  data.frame()%>% write_clip()


## CI coverage of T ----

colnames(table_T_CIcoverage_secondary)<- c("rep", "missingMethod_Mechanism", "yes_no_coverage")
table_T_CIcoverage_secondary$yes_no_coverage<- as.numeric(table_T_CIcoverage_secondary$yes_no_coverage)

table_T_CIcoverage_secondary$missingMethod_Mechanism<- factor(table_T_CIcoverage_secondary$missingMethod_Mechanism, 
                                                              levels= c("unadjusted",
                                                                        "noMiss", 
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


table_T_CIcoverage_secondary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(yes_no_coverage))%>%
  data.frame()





## Coefficient of X1*T ----

colnames(table_X1T_secondary)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_X1T_secondary$estimate<- as.numeric(table_X1T_secondary$estimate)
table_X1T_secondary$SE<- as.numeric(table_X1T_secondary$SE)

table_X1T_secondary$missingMethod_Mechanism<- factor(table_X1T_secondary$missingMethod_Mechanism, 
                                                     levels= c("unadjusted",
                                                               "noMiss", 
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
table_X1T_secondary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(estimate), sd(estimate), mean(SE), mse(estimate, 2))%>%
  data.frame()


## CI coverage of X1T ----

colnames(table_X1T_CIcoverage_secondary)<- c("rep", "missingMethod_Mechanism", "yes_no_coverage")
table_X1T_CIcoverage_secondary$yes_no_coverage<- as.numeric(table_X1T_CIcoverage_secondary$yes_no_coverage)

table_X1T_CIcoverage_secondary$missingMethod_Mechanism<- factor(table_X1T_CIcoverage_secondary$missingMethod_Mechanism, 
                                                                levels= c("unadjusted",
                                                                          "noMiss", 
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


table_X1T_CIcoverage_secondary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(yes_no_coverage))%>%
  data.frame()



## Coefficient of X2*T ----

colnames(table_X2T_secondary)<- c("rep", "missingMethod_Mechanism", "model", "estimand", "estimate", "SE")
table_X2T_secondary$estimate<- as.numeric(table_X2T_secondary$estimate)
table_X2T_secondary$SE<- as.numeric(table_X2T_secondary$SE)

table_X2T_secondary$missingMethod_Mechanism<- factor(table_X2T_secondary$missingMethod_Mechanism, 
                                                     levels= c("unadjusted",
                                                               "noMiss", 
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
table_X2T_secondary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(estimate), sd(estimate), mean(SE), mse(estimate, 2))%>%
  data.frame()


## CI coverage of X2T ----

colnames(table_X2T_CIcoverage_secondary)<- c("rep", "missingMethod_Mechanism", "yes_no_coverage")
table_X2T_CIcoverage_secondary$yes_no_coverage<- as.numeric(table_X2T_CIcoverage_secondary$yes_no_coverage)

table_X2T_CIcoverage_secondary$missingMethod_Mechanism<- factor(table_X2T_CIcoverage_secondary$missingMethod_Mechanism, 
                                                                levels= c("unadjusted",
                                                                          "noMiss", 
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


table_X2T_CIcoverage_secondary%>%
  group_by(missingMethod_Mechanism)%>%
  summarise(mean(yes_no_coverage))%>%
  data.frame()


