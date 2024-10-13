library(wakefield)  # library for generating random dataset
set.seed(12)

n = 1150  # sample size


#########################generate covariates:
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
# protein intake
protein = rnorm(n, 50, 10)
while (sum(protein <= 0) > 0) {  # everyone in the data set eat at least some protein
  protein_leq_0 = protein <= 0
  protein[protein_leq_0] = rnorm(sum(protein_leq_0), 50, 10)
}
# calorie intake
calorie_intake = rnorm(n, 1000 - protein * 4, 100)
while (sum(calorie_intake <= 200)) {  # minimum calorie intake other than protein
  calorie_leq_200 = calorie_intake <= 200
  calorie_intake[calorie_leq_200] = rnorm(sum(calorie_leq_200), 1000 - protein[calorie_leq_200] * 4, 100)
}
total_calorie_intake = protein * 4 + calorie_intake

df = data.frame(
  ages = ages,
  gender = ifelse(gender == "M", 1, 0),
  surgery = ifelse(surgery_type == "RYGB", 1, 0),
  BMI = BMI,
  protein_intake = protein,
  total_calorie_intake = total_calorie_intake,
  protein_calorie = protein * total_calorie_intake
)


#########################################generate outcome
# regression coefficient
beta = c(ages = 0.5, 
         gender = 3, 
         surgery = 0.01, 
         BMI = 1.7, 
         protein_intake = -0.12, 
         total_calorie_intake = -0.08,
         protein_calorie = -0.0004)

# Generate event (muscle mass loss)
baseline_hazard_T1 = 20
hazard_T1 = baseline_hazard_T1 * exp(as.matrix(df) %*% beta)

baseline_hazard_T3 = 50
hazard_T3 = baseline_hazard_T3 * exp(as.matrix(df) %*% beta)

baseline_hazard_T6 = 80
hazard_T6 = baseline_hazard_T6 * exp(as.matrix(df) %*% beta)

# Time-to-event data with event indicator
event_T1 = ifelse((1 - exp(-hazard_T1 * 1)) > 0.5, 1, 0)
event_T3 = ifelse((1 - exp(-hazard_T3 * 1)) > 0.5, 1, 0)
event_T6 = ifelse((1 - exp(-hazard_T6 * 1)) > 0.5, 1, 0)

# check
sum((event_T1 == 1) & (event_T3 != 1))
sum((event_T1 == 1) & (event_T6 != 1))
sum((event_T3 == 1) & (event_T6 != 1))
# check
mean(event_T1)
mean(event_T3)
mean(event_T6)






############################################data for modeling

# time
time = numeric(n)
time[event_T1 == 1] = 1
time[(event_T3 == 1) & (event_T1 != 1)] = 3
time[(event_T6 == 1) & (event_T3 != 1) & (event_T1) != 1] = 6
time[(event_T6 == 0)] = 6
# status
status = event_T6

# censor = 0 for each interval
df_pois = cbind(df, time, status) 
df_pois_censor_T1 = df_pois[event_T1 == 0, ]
df_pois_censor_T1$time = 1
df_pois_censor_T1$status = 0
df_pois_censor_T3 = df_pois[event_T3 == 0 & event_T1 == 0, ]
df_pois_censor_T3$time = 3
df_pois_censor_T3$status = 0
df_pois = rbind(df_pois, df_pois_censor_T1, df_pois_censor_T3)
df_pois$time = as.factor(df_pois$time)
#########################################################
#################### modeling ###########################
#########################################################
library(survival)
## Cox model
cox = coxph(Surv(time, status) ~ ., df, ties = 'breslow')
## Poisson model
pois = glm(status ~ 0 + ., family = poisson, data = df_pois)

summary(cox)
summary(pois)

## Cumulative hazard
baseline <- basehaz(cox, centered = FALSE)
cumhaz_Cox = unique(baseline$hazard)
cumhaz_Poisson = unname(cumsum(exp(pois$coefficients[-1:-7])))

## Survival probability
survival_Cox = exp(-cumhaz_Cox)
survival_Poisson = exp(-cumhaz_Poisson)
survival_Cox
survival_Poisson
