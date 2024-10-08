library(wakefield)  # library for generating random dataset

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
protein = rnorm(n, 50, 20)
while (sum(protein <= 0) > 0) {  # everyone in the data set eat at least some protein
  protein_leq_0 = protein <= 0
  protein[protein_leq_0] = rnorm(sum(protein_leq_0), 50, 20)
}
# # calorie intake
# calorie_intake = rnorm(n, 1000 - protein * 4, 300)
# while (sum(calorie_intake <= 200)) {
#   calorie_leq_200 = calorie_intake <= 200
#   calorie_intake[calorie_leq_200] = rnorm(sum(calorie_leq_200), 1000 - protein[calorie_leq_200] * 4, 300)
# }
# total_calorie_intake = protein * 4 + calorie_intake
# 
# ## T0
# # weight
# weight_T0 = BMI_T0 * height^2
# # Fat Free Mass
# fat_pct_T0 = rnorm(n, 50, 5) / 100
# FFM_T0 = (1 - fat_pct_T0) * weight_T0
# 
# ## T1
# # weight
# weight_T2 = weight_T0 - 
# # Fat Free Mass
# fat_pct_T0 = rnorm(n, 50, 5) / 100
# FFM_T0 = (1 - fat_pct_T0) * weight_T0
# # FFML
# 
# ## T3
# # weight
# 
# # Fat Free Mass
# 
# # FFML
# 
# ## T6
# # weight
# 
# # Fat Free Mass
# 
# # FFML


### putting them together
df = data.frame(
  ages = ages,
  gender = ifelse(gender == "M", 1, 0),
  surgery = ifelse(surgery_type == "RYGB", 1, 0),
  BMI = BMI,
  protein_intake = protein
)

## generate binary outcome with logistic regression
beta = c(0.06, -2, 2.5, 0.1, -0.2)
z = -2 + as.matrix(df) %*% beta  # intercept: 2 for 50-50, -2 for rare events
p = 1 / (1 + exp(-z)) 
hist(p)
EFFML = rbinom(n, 1, p)      
hist(EFFML)
df$EFFML = EFFML

## generate outcome with linear regression


## generate outcome with cox model



#########################################################
########### cox regression ##############################
#########################################################
library(survival)
time = rep(6, n)

summary(glm(EFFML ~ ., data = df, family = "binomial"))

summary(glm(EFFML ~ ., data = df, family = "poisson"))
summary(coxph(Surv(time, df$EFFML) ~ ., data = df, ties = "breslow"))

