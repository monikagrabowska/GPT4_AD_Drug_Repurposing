library(glue)
library(lubridate)
library(tidyverse)
library(vroom)
library(MatchIt)
library(optmatch)
library(boot)
library(lmtest)
library(sandwich)
library(survival)
library(survminer)
library(marginaleffects)
library(finalfit)

set.seed(54321)

# directory containing data files
files <- list.files(path=".", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in seq_along(files)) {
  x <- files[[i]]
  print(x)

  data <- vroom::vroom(x,
                       .name = janitor::make_clean_names, 
                       col_types = cols())

  # standardize column names for drugs with comorbidities
  if (ncol(data) == 16) {
    names(data)[8] <- "drug_status" 
    names(data)[12] <- "comorbidity_status" 
 
    # total number drug-exposed
    total_drug <- nrow(data %>% filter(drug_status == 1))

    if (!(total_drug == 0)) { 
      # matching - sex, race, ehr length after 65, comorbidities  
      m.out <- matchit(drug_status ~ gender_concept_id + race_concept_id + ehr_length + comorbidity_status,
                       data = data,
                       method = "nearest",
                       ratio = 2, caliper = 0.1, verbose = TRUE)
      print(summary(m.out))

      # estimating effects after matching 
      md <- match.data(m.out)
    }
  }
  # standardize column names for drugs without comorbidities
  else {
    names(data)[8] <- "drug_status"
 
    # total number drug-exposed
    total_drug <- nrow(data %>% filter(drug_status == 1))

    if (!(total_drug == 0)) {
      # matching - sex, race, ehr length after 65    
      m.out <- matchit(drug_status ~ gender_concept_id + race_concept_id + ehr_length,
                       data = data,
                       method = "nearest",
                       ratio = 2, caliper = 0.1, verbose = TRUE)
      print(summary(m.out))

      # estimating effects after matching 
      md <- match.data(m.out)  
    }
  }

  total_drug <- nrow(data %>% filter(drug_status == 1))
  if (!(total_drug == 0)) {
    matched_drug_ad <- nrow(md %>% filter(drug_status == 1 & ad == 1)) # number ad/exposed
    matched_drug_noad <- nrow(md %>% filter(drug_status == 1 & ad == 0)) # number no ad/exposed
    matched_nodrug_ad <- nrow(md %>% filter(drug_status == 0 & ad == 1)) # number ad/unexposed
    matched_nodrug_noad <- nrow(md %>% filter(drug_status == 0 & ad == 0)) # number no ad/unexposed

    # cox regression
    fit.cox <- coxph(Surv(ehr_length, ad) ~ drug_status, data = md, robust = TRUE, 
            weights = weights, cluster = subclass)

    print(summary(fit.cox))
  }
}
