library(wakefield)  # library for generating random dataset
library(survival)
library(ggplot2)
library(survminer)  # library for plotting hazard and survival

set.seed(12)

n = 1150  # sample size

################################################################################
################################################################################
############################ generate covariates ###############################
################################################################################
################################################################################
# age
ages = age(n, x = 20:65, prob = NULL)  # age range from 20 to 65
# gender
n_male = rbinom(1, n, 0.2)
n_female = n - n_male
gender = c(rep("M", n_male), rep("F", n_female))
# height
height_m = height_cm(n_male, 185, 7)
height_f = height_cm(n_female, 170, 6)
height = c(height_m, height_f) / 100   # height in meters
# surgery types. SG: sleeve gastrectomy (SG). RYGB: Roux-en-Y gastric bypass
n_SG = rbinom(1, n, 0.5)
n_RYGB = n - n_SG
surgery_type = sample(c(rep("SG", n_SG), rep("RYGB", n_RYGB)))
# BMI
BMI = rnorm(n, 45, 5)
while (sum(BMI <= 30) > 0) {  # everyone in the dataset is obese
  not_obese = BMI <= 30
  BMI[not_obese] = rnorm(sum(not_obese), 45, 5)  
}
# protein intake, in grams
protein = rnorm(n, 50, 10)
while (sum(protein <= 0) > 0) {  # everyone in the data set eat at least some protein
  protein_leq_0 = protein <= 0
  protein[protein_leq_0] = rnorm(sum(protein_leq_0), 50, 10)
}
# calorie intake except protein, in kcal
calorie_intake = rnorm(n, 1000 - protein * 4, 100)
while (sum(calorie_intake <= 200)) {  # minimum calorie intake other than protein
  calorie_leq_200 = calorie_intake <= 200
  calorie_intake[calorie_leq_200] = rnorm(sum(calorie_leq_200), 1000 - protein[calorie_leq_200] * 4, 100)
}

df = data.frame(
  ages = ages,
  gender = ifelse(gender == "M", 1, 0),
  surgery = ifelse(surgery_type == "RYGB", 1, 0),
  BMI = BMI,
  protein_intake = protein,
  calorie_intake = calorie_intake
)


################################################################################
############################# generate outcome #################################
################################################################################
# regression coefficient, these are completely arbitrary
beta = c(ages = 0.05, 
         gender = 0.5, 
         surgery = 0, 
         BMI = 0.25, 
         protein_intake = -0.01, 
         calorie_intake = -0.02)

noise = rnorm(n, 0, 0.5) # add noise so that the problem is not linearly separable
                         # (for logistic regression to converge)

hazard = function(baseline_hazard, df, beta, noise) {
  baseline_hazard * exp(as.matrix(df) %*% beta) + noise
}

baseline_hazards = c(0.5, 1.5, 4)  # set hazard rate so that the period has ~ 10% new events
baseline_hazards_rare = c(0.1, 0.2, 0.3)  # set hazard rate so the the total events is about 3% 
                                          # (rare enough for the binomial to approximate Poisson)

hazard_T1 = hazard(baseline_hazards[1], df, beta, noise)
hazard_T3 = hazard(baseline_hazards[2], df, beta, noise)
hazard_T6 = hazard(baseline_hazards[3], df, beta, noise)
# for rare events (3%), uncomment the following code
# hazard_T1 = hazard(baseline_hazards_rare[1], df, beta, noise)
# hazard_T3 = hazard(baseline_hazards_rare[2], df, beta, noise)
# hazard_T6 = hazard(baseline_hazards_rare[3], df, beta, noise)


# the more the following results are similar, 
# the more the cox&poisson estimates are closer to the original coefficients
sum(trunc(hazard_T1))
sum(trunc(hazard_T1) > 0)
# the more the following results are similar, 
# the more the cox&poisson estimates are closer to the original coefficients
sum(trunc(hazard_T3))
sum(trunc(hazard_T3) > 0)
# the more the following results are similar, 
# the more the cox&poisson estimates are closer to the original coefficients
sum(trunc(hazard_T6))
sum(trunc(hazard_T6) > 0)


trunc_hazard = function(hazard_T1, hazard_T3, hazard_T6) {
  # event indicator
  event_T1 = ifelse(trunc(hazard_T1) > 0, 1, 0)
  event_T3 = ifelse(trunc(hazard_T3) > 0, 1, 0)
  event_T6 = ifelse(trunc(hazard_T6) > 0, 1, 0)
  # check. All should be 0.
  print(sum((event_T1 == 1) & (event_T3 != 1)))
  print(sum((event_T1 == 1) & (event_T6 != 1)))
  print(sum((event_T3 == 1) & (event_T6 != 1)))
  # check. Since hazard_T1 < hazard_T2 < hazard_T3, mean values should increase
  print(mean(event_T1))
  print(mean(event_T3))
  print(mean(event_T6))
  
  list(event_T1 = event_T1, event_T3 = event_T3, event_T6 = event_T6)
}

events = trunc_hazard(hazard_T1, hazard_T3, hazard_T6)
event_T1 = events$event_T1
event_T3 = events$event_T3
event_T6 = events$event_T6


################################################################################
############################# data for modeling ################################
################################################################################

# time
create_time = function(event_T1, event_T3, event_T6) {
  time = numeric(n)
  time[event_T1 == 1] = 1
  time[(event_T3 == 1) & (event_T1 != 1)] = 3
  time[(event_T6 == 1) & (event_T3 != 1) & (event_T1) != 1] = 6
  time[(event_T6 == 0)] = 6
  time
}
time = create_time(event_T1, event_T3, event_T6)
# status
status = event_T6

# make the data set for Poisson and logistic regression
# split the dataset at times when any event happened, and set nonevent to 0 (censored)
make_dataset_Poisson = function(df_time_event, event_T1, event_T3, event_T6) {
  df_pois_censor_T1 = df_time_event[event_T1 == 0, ]
  df_pois_censor_T1$time = 1
  df_pois_censor_T1$status = 0
  df_pois_censor_T3 = df_time_event[event_T3 == 0 & event_T1 == 0, ]
  df_pois_censor_T3$time = 3
  df_pois_censor_T3$status = 0
  df_pois = rbind(df_time_event, df_pois_censor_T1, df_pois_censor_T3)
  df_pois$time = as.factor(df_pois$time)
  df_pois
}
df_pois = make_dataset_Poisson(cbind(df, time, status), event_T1, event_T3, event_T6)

mean(df_pois$status)  ## overall event rate
################################################################################
################################################################################
################################ modeling ######################################
################################################################################
################################################################################
# Cox model
cox = coxph(Surv(time, status) ~ ., df, ties = 'breslow')
# Poisson model
pois = glm(status ~ 0 + ., family = poisson, data = df_pois)
# logistic regression
lr = glm(status ~ 0 + ., family = binomial, data = df_pois)

summary(cox)
summary(pois)
summary(lr)

################################################################################
####################### Comparing cumhaz and survival ##########################
################################################################################
## Cumulative hazard: baseline
baseline = basehaz(cox, centered = F)
cumhaz_Cox = unique(baseline$hazard)
cumhaz_Poisson = unname(cumsum(exp(pois$coefficients[-1:-6])))
abs(cumhaz_Cox - cumhaz_Poisson) < 1e-6

## Cumulative hazard: centered
# Poisson model
df_centered = colMeans(df)
df_centered["gender"] = 0
df_centered["surgery"] = 0
cumhaz_Poi_c = cumsum(exp(pois$coefficients[-1:-6]) * as.vector(exp(df_centered %*% pois$coefficients[1:6])))

H = data.frame(
  time = c(0, 1, 3, 6),
  cumhaz = c(0, cumhaz_Poi_c)
)
ggplot(H, aes(x = time, y = cumhaz)) +
  geom_step(color = "blue") +    # Stair plot using geom_step
  labs(title = "Cumulative Hazard - Poisson Model",
       x = "Time",
       y = "Cumulative Hazard") +
  ylim(c(0, 0.32)) +
  theme_minimal()
# Cox model
surv_fit = survfit(cox)
ggsurvplot(surv_fit, 
           data = df,
           conf.int = TRUE,          # Add confidence intervals
           fun = "cumhaz",
           palette = "red",          # Customize color
           ggtheme = theme_minimal(), # Minimal theme
           title = "Cumulative Hazard - Cox Model",  # Title for the plot
           xlab = "Time",             # X-axis label
           ylab = "Cumulative Hazard",
           legend = "none")  # Y-axis label

## Survival
# Poisson
survival_Poisson = exp(-cumhaz_Poi_c)
S = data.frame(
  time = c(0, 1, 3, 6),
  survival = c(1, survival_Poisson)
)
ggplot(S, aes(x = time, y = survival)) +
  geom_step(color = "blue") +    # Stair plot using geom_step
  labs(title = "Survival Curve - Poisson Model",
       x = "Time",
       y = "Survival Probability") +
  ylim(c(0, 1)) +
  theme_minimal()
# Cox
ggsurvplot(surv_fit, 
           data = df,
           conf.int = TRUE,          # Add confidence intervals
           palette = "red",          # Customize color
           ggtheme = theme_minimal(), # Minimal theme
           title = "Survival Curve - Cox Model",  # Title for the plot
           xlab = "Time",             # X-axis label
           ylab = "Survival Probability", # Y-axis label
           legend = "none") 





################################################################################
################################################################################
############################ Stratified Model ##################################
################################################################################
################################################################################
# note that for convenience, no new dataset is created
cox = coxph(Surv(time, status) ~ ages + 
                                 strata(gender) + 
                                 BMI + 
                                 protein_intake + 
                                 calorie_intake,
            data = df, 
            ties = 'breslow')
pois = glm(status ~ 0 +
                    ages + 
                    time +
                    gender:time +
                    BMI + 
                    protein_intake + 
                    calorie_intake, 
           family = poisson, 
           data = df_pois)

summary(cox)
summary(pois)
################################################################################
################################################################################
######################### time varying coefficient #############################
################################################################################
################################################################################

## create data
# hazard with changing BMI
BMI_T1 = BMI - rnorm(n, 5, 1)
BMI_T3 = BMI_T1 - rnorm(n, 4, 1)
df_copy = df
baseline_hazards = c(0.2, 5, 20)
hazard_T1 = hazard(baseline_hazards[1], df_copy, beta, noise)
df_copy$BMI= BMI_T1
hazard_T3 = hazard(baseline_hazards[2], df_copy, beta, noise)
df_copy$BMI = BMI_T3
hazard_T6 = hazard(baseline_hazards[3], df_copy, beta, noise)
# events by time
events = trunc_hazard(hazard_T1, hazard_T3, hazard_T6)
event_T1 = events$event_T1
event_T3 = events$event_T3
event_T6 = events$event_T6
# time
time = create_time(event_T1, event_T3, event_T6)
# status
status = event_T6


## create data for Poisson model
df = cbind(df, time, status, BMI_T1, BMI_T3)
df$BMI_T0 = df$BMI
df_pois = make_dataset_Poisson(df, event_T1, event_T3, event_T6)
df_pois[df_pois$time == 1, "BMI"] = df_pois[df_pois$time == 1, "BMI_T0"]
df_pois[df_pois$time == 3, "BMI"] = df_pois[df_pois$time == 3, "BMI_T1"]
df_pois[df_pois$time == 6, "BMI"] = df_pois[df_pois$time == 6, "BMI_T3"]

## create data for cox model
start_time = numeric(nrow(df_pois))
start_time[df_pois$time == 1] = 0
start_time[df_pois$time == 3] = 1
start_time[df_pois$time == 6] = 3

end_time = numeric(nrow(df_pois))
end_time[df_pois$time == 1] = 1
end_time[df_pois$time == 3] = 3
end_time[df_pois$time == 6] = 6

## cox model
cox = coxph(Surv(start_time, end_time, status) ~ BMI, data = df_pois, ties = 'breslow')
## Poisson model
pois = glm(status ~ 0 + time + BMI, family = poisson, data = df_pois)

summary(cox)
summary(pois)
