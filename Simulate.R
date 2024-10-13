library(wakefield)  # library for generating random dataset
set.seed(12)

n = 1150  # sample size

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

## T0
# weight
weight_T0 = BMI * height^2
# Fat Free Mass
fat_pct_T0 = rnorm(n, 50, 5) / 100
FM_T0 = weight_T0 * fat_pct_T0
FFM_T0 = weight_T0 - FM_T0
# data
df = data.frame(
  ages = ages,
  gender = ifelse(gender == "M", 1, 0),
  surgery = ifelse(surgery_type == "RYGB", 1, 0),
  BMI = BMI,
  protein_intake = protein,
  total_calorie_intake = total_calorie_intake,
  protein_calorie = protein * total_calorie_intake,
  time = NA
)
# regression coefficient
beta_weight_pct_loss = c(age = 0.02, 
                         gender = 0.1, 
                         surgery = 0.01, 
                         BMI = 0.1, 
                         protein_intake = -0.05, 
                         total_calorie_intake = -0.01,
                         protein_calorie = 0,
                         time = 15,
                         timesq = -0.01)
beta_fat_pct_loss = c(age = 0.016 - 0.002, 
                      gender = 0.08 + 0.01, 
                      surgery = 0.008, 
                      BMI = 0.08 + 0.05, 
                      protein_intake = -0.04 - 0.01, 
                      total_calorie_intake = -0.008,
                      protein_calorie = 0 - 0.000005,
                      time = 10,
                      timesq = -0.2)
# df$time = NULL
# df$timesq = NULL
# range((as.matrix(df) %*% beta_weight_pct_loss))
# range(as.matrix(df) %*% beta_fat_pct_loss)
for (t in 1:3) {
  df$time = t
  df$timesq = df$time^2
  weight_pct_loss = (as.matrix(df) %*% beta_weight_pct_loss) - 2 + rnorm(n, 0, 0.05)
  print(all(weight_pct_loss > 0))
  print(median(weight_pct_loss))
  
  fat_pct_loss = (as.matrix(df) %*% beta_fat_pct_loss) - 1 + rnorm(n, 0, 0.05) 
  print(all(fat_pct_loss > 0))
  print(median(fat_pct_loss))
  
  print(all(weight_pct_loss - fat_pct_loss > 0))
  print(median(fat_pct_loss / weight_pct_loss))
  print(" ")
}


EFFML = function(df, t, 
                 beta_weight_pct_loss, beta_fat_pct_loss, 
                 weight_T0, FFM_T0, 
                 weight_previous, fat_previous) {
  df$time = t
  df$timesq = df$time^2
  d = as.matrix(df)
  # weight
  weight_pct_loss = (as.matrix(df) %*% beta_weight_pct_loss) - 2 
  weight_next = weight_previous * (1 - weight_pct_loss / 100)
  # Fat Free Mass
  fat_pct_loss = (as.matrix(df) %*% beta_fat_pct_loss) - 1 
  fat_next = fat_previous * (1 - fat_pct_loss / 100)
  FFM_next = weight_next - fat_next
  # FFML
  FFML_next = (FFM_next - FFM_T0) / (weight_next - weight_T0)
  # excessive FFML
  EFFML = abs(FFML_next) > 0.35
  list(weight_next, fat_next, FFML_next) 
}

l1 = EFFML(df, 1, 
      beta_weight_pct_loss, beta_fat_pct_loss, 
      weight_T0, FFM_T0, 
      weight_T0, FM_T0)

weight_T1 = l1[[1]]
fat_T1 = l1[[2]]
mean(l1[[3]])

l3 = EFFML(df, 2, 
      beta_weight_pct_loss, beta_fat_pct_loss, 
      weight_T0, FFM_T0, 
      weight_T1, fat_T1)

weight_T3 = l3[[1]]
fat_T3 = l3[[2]]
mean(l3[[3]])

l6 = EFFML(df, 3, 
      beta_weight_pct_loss, beta_fat_pct_loss, 
      weight_T0, FFM_T0, 
      weight_T3, fat_T3)

mean(l6[[3]])
## T3
# time difference is 2
df$time = 2
df$timesq = df$time^2
d = as.matrix(df)
# weight
weight_pct_loss = (as.matrix(df) %*% beta_weight_pct_loss) + 10 + rnorm(n, 0, 0.1)
weight_T3 = weight_T1 * (1 - weight_pct_loss / 100)
# Fat Free Mass
fat_pct_loss = (as.matrix(df) %*% beta_fat_pct_loss) + 6 + rnorm(n, 0, 0.1)
fat_pct_T3 = fat_pct_T1 - fat_pct_loss / 100
FFM_T3 = (1 - fat_pct_T3) * weight_T3
# FFML
FFML_T3 = (FFM_T3 - FFM_T0) / (weight_T3 - weight_T0)
# excessive FFML
EFFML_T3 = abs(FFML_T3) > 0.25
mean(EFFML_T3)

## T6
# time difference is 3
df$time = 3
df$timesq = df$time^2
d = as.matrix(df)
# weight
weight_pct_loss = (as.matrix(df) %*% beta_weight_pct_loss) + 8 + rnorm(n, 0, 0.1)
weight_T6 = weight_T3 * (1 - weight_pct_loss / 100)
# Fat Free Mass
fat_pct_loss = (as.matrix(df) %*% beta_fat_pct_loss) + 4 + rnorm(n, 0, 0.1)
fat_pct_T6 = fat_pct_T3 - fat_pct_loss / 100
FFM_T6 = (1 - fat_pct_T6) * weight_T6
# FFML
FFML_T6 = (FFM_T6 - FFM_T0) / (weight_T6 - weight_T0)
# excessive FFML
EFFML_T6 = abs(FFML_T6) > 0.25
mean(EFFML_T6)


### putting them together


## generate outcome with linear regression




#########################################################
#################### modeling ###########################
#########################################################
library(survival)
## Cox model
kfit1 <- coxph(Surv(time, status) ~ age + sex, kidney, ties='breslow')
## Poisson model
utime <- sort(unique(kidney$time[kidney$status==1])) # unique deaths
kdata3 <- survSplit(Surv(time, status) ~., data=kidney, cut=utime, episode="interval")
kdata3 <- subset(kdata3, time == c(utime,0)[interval]) # remove partials
kfit3 <- glm(status ~ age + sex + factor(interval) - 1, family=poisson, data=kdata3)

## Cumulative hazard
baseline <- basehaz(cox_, centered = FALSE)
cumhaz_Cox = unique(baseline$hazard)
cumhaz_Poisson = unname(cumsum(exp(kfit3$coefficients[-1:-2])))

## Survival probability
survival_Cox = exp(-cumhaz_Cox)
survival_Poisson = exp(-cumhaz_Poisson)
