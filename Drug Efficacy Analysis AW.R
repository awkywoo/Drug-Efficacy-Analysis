#Load packages
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(tableone)
library(MatchIt)
install.packages("naniar")
library(naniar)


#Create data frame
dfDemo <- read_csv("Patient_characteristics.csv")
dfDuration <- read_csv("Event_duration.csv")

#Merge data frame by id and treatment
dfDrug <- merge(x = dfDemo, y = dfDuration, by = c("patient_id", "treatment_variable"))

#Dichotomize age as predictive variable for later
dfDrug <- dfDrug %>% mutate(age_group = ifelse(age >= 65, "65 & above", "below 65"))
dfDrug$age_group <- factor(dfDrug$age_group)
dfDrug$age_group <- recode_factor(dfDrug$age_group,
                                  `below 65` = 1,
                                  `65 & above` = 2)

#Recode variables from chr to factor
for (col in 5:27) {
  dfDrug[, colnames(dfDrug)[col]] <- recode_factor(dfDrug[, colnames(dfDrug)[col]],
                                                   `Yes` = 1,
                                                   `No` = 2)
}

dfDrug$treatment_variable <- factor(dfDrug$treatment_variable)
dfDrug$treatment_variable <- recode_factor(dfDrug$treatment_variable,
                                  `Drug_A` = 1,
                                  `Drug_B` = 2)

#Inspect data frame
str(dfDrug)
sapply(lapply(dfDrug, unique), length)
gg_miss_var(dfDrug, show_pct = TRUE)

summary(dfDrug)
table1 <- CreateTableOne(vars = c('age_group', 'sex', 'Bleeding_event',
                                  'other_drugs_1','other_drugs_2','other_drugs_3','other_drugs_4',
                                  'other_drugs_5','other_drugs_6','other_drugs_7','other_drugs_8',
                                  'diagnosis_1','diagnosis_2','diagnosis_3','diagnosis_4',
                                  'diagnosis_5','diagnosis_6','diagnosis_7','diagnosis_8',
                                  'diagnosis_9','diagnosis_10','diagnosis_11','diagnosis_12',
                                  'diagnosis_13','diagnosis_14','diagnosis_15',
                                  'lab_1','lab_2','lab_3','lab_4',
                                  'lab_5','lab_6','lab_7','lab_8',
                                  'Diag_Score_1', 'Diag_Score_2'), 
                         data = dfDrug, 
                         factorVars = c('age_group', 'sex', 'Bleeding_event',
                                        'other_drugs_1','other_drugs_2','other_drugs_3','other_drugs_4',
                                        'other_drugs_5','other_drugs_6','other_drugs_7','other_drugs_8',
                                        'diagnosis_1','diagnosis_2','diagnosis_3','diagnosis_4',
                                        'diagnosis_5','diagnosis_6','diagnosis_7','diagnosis_8',
                                        'diagnosis_9','diagnosis_10','diagnosis_11','diagnosis_12',
                                        'diagnosis_13','diagnosis_14','diagnosis_15', 'Diag_Score_1', 'Diag_Score_2'),
                         strata = 'treatment_variable',
                         smd = TRUE)
table1 <- print(table1, 
                showAllLevels = TRUE,
                smd = TRUE)
table1
# Propensity score matching
# logit model
m_ps <- glm(treatment_variable ~ age_group + sex
            +other_drugs_2 +other_drugs_4
            +other_drugs_5 +other_drugs_6+other_drugs_7
            +diagnosis_1+diagnosis_2+diagnosis_3
            +diagnosis_4+diagnosis_7+diagnosis_9
            +diagnosis_10+diagnosis_11+diagnosis_12
            +diagnosis_13+diagnosis_14+diagnosis_15
            +lab_3+lab_5+lab_6+lab_7+lab_8
            +Diag_Score_1+Diag_Score_2,
            family = binomial(), data = dfDrug)
summary(m_ps)

#propensity score calculation
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     treatment_variable = m_ps$model$treatment_variable)
head(prs_df)

#region of common support for treatment status
labs <- paste("treatment given:", c("Drug_A", "Drug_B"))
prs_df %>%
  mutate(treatment = ifelse(treatment_variable == 1, labs[1], labs[2])) %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~treatment_variable) +
  xlab("Probability of being given either treatment") +
  theme_bw()

#matching using matchit, covariates significant from logit model
match.it <- matchit(treatment_variable ~ age_group + sex
                    +other_drugs_5+diagnosis_2+diagnosis_3
                    +diagnosis_11+diagnosis_13+Diag_Score_2, 
                    data = dfDrug, method="nearest", distance = "logit", ratio=1)
a <- summary(match.it)

knitr::kable(a$nn, digits = 2, align = 'c', 
      caption = 'Matched sample sizes')

#Create new data frame for matched data
df.match <- match.data(match.it)[1:ncol(dfDrug)]

#examine covariate balance in matched sample
table2 <- CreateTableOne(vars = c('age_group', 'sex',
                                  'other_drugs_5','diagnosis_2','diagnosis_3',
                                  'diagnosis_11','diagnosis_13','Diag_Score_2'), 
                         data = df.match, 
                         factorVars = c('age_group', 'sex',
                                        'other_drugs_5','diagnosis_2','diagnosis_3',
                                        'diagnosis_11','diagnosis_13','Diag_Score_2'),
                         strata = 'treatment_variable')
table2

#treatment effect without matched data
lm_treat <- lm(Bleeding_event ~ treatment_variable, data = dfDrug)
summary(lm_treat)

#treatment effect with matched data
lm_treat1 <- lm(Bleeding_event ~ treatment_variable, data = df.match)
summary(lm_treat1)

# Fit survival data using the Kaplan-Meier method

#Survival with treatment
surv_object <- Surv(time = df.match$duration_in_years, event = df.match$Bleeding_event)
surv_object 

#Overall survival
survfit1 <- survfit(surv_object ~ 1, data = df.match)
survfit1 # 17036 patients, 3164 events

# Survival curve by treatment
fit1 <- survfit(surv_object ~ treatment_variable, data = df.match)
summary(fit1)

# Plot survival curve
ggsurvplot(fit1, data = df.match, 
           pval = TRUE,
           conf.int = TRUE, 
           legend.labs=c("Drug A", "Drug B"),
           risk.table = TRUE,
           surv.median.line = "hv",
           ggtheme = theme_minimal())

#apply univariate coxph function to multiple co-variates
covariates <- c('age_group', 'sex',
                'other_drugs_5','diagnosis_2','diagnosis_3',
                'diagnosis_11','diagnosis_13','Diag_Score_2')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('surv_object~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df.match)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         return(exp(cbind(coef(x),confint(x))))
                       })
univ_results
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#The output shows the regression beta coefficients, 
#the effect sizes (given as hazard ratios) and 
#statistical significance for each of the variables 
#in relation to overall survival. 
#Each factor is assessed through separate univariate Cox regressions.

## Fit a Cox proportional hazards model

fit.coxph <- coxph(surv_object ~treatment_variable+age_group + sex
                  +other_drugs_5+diagnosis_2+diagnosis_3
                  +diagnosis_11+diagnosis_13+Diag_Score_2,
                  data = df.match)
summary(fit.coxph)

#Assess model fit
pscl::pR2(fit.coxph)["McFadden"]

#Assess multicollinearity
car::vif(fit.coxph)

ggforest(fit.coxph, data = df.match)

# Visualize estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(fit.coxph, data = df.match), palette = "#2E9FDF",
           ggtheme = theme_minimal())

# Create the new data, covariates fixed to average/ lowest level  
data = df.match
treat_df <- with(df.match,
               data.frame(treatment_variable = c(1, 2),
                          age_group = c(1,1),
                          sex = c(1,1),
                          other_drugs_5 = c(1,1),
                          diagnosis_2 = c(1,1),
                          diagnosis_3 = c(1,1),
                          diagnosis_11 = c(1,1),
                          diagnosis_13 = c(1,1),
                          Diag_Score_2 = c(0,0)
               )
      )
str(treat_df)

#Survival curves
treat_df[sapply(treat_df, is.numeric)] <- lapply(treat_df[sapply(treat_df, is.numeric)], as.factor)

fit2 <- survfit(fit.coxph, newdata = treat_df)
ggsurvplot(fit2, data = treat_df, 
           conf.int = TRUE, 
           legend.labs=c("Drug A", "Drug B"),
           ggtheme = theme_minimal())

