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

# generate surgery types. SG: sleeve gastrectomy (SG). RYGB: Roux-en-Y gastric bypass
n_SG = rbinom(1, n, 0.5)
n_RYGB = n - n_SG
surgery_type = sample(c(rep("SG", n_SG), rep("RYGB", n_RYGB)))

# BMI
BMI_T0 = rnorm(n, 45, 5)
not_obese = BMI_T0 <= 30
BMI_T0[not_obese] = rnorm(sum(not_obese), 45, 1)  # everyone in the dataset is obese

# weight
weight_T0 = BMI_T0 * height^2

# # Fat Free Mass
# fat_pct_T0 = rnorm(n, 50, 5) / 100
# FFM_T0 = (1 - fat_pct_T0) * weight_T0

### T6
# protein intake
protein = rnorm(n, 50, 20)
protein_less_than_0 = protein < 0
protein[protein_less_than_0] = rnorm(sum(protein_less_than_0), 50, 1)

## calorie intake
# total_calorie_intake = protein * 4 + rnorm(n, 1000 - protein * 4, 300)

# post operative fat free mass

# post operative weight

# FFML

### putting them together
df = data.frame(
  ages = ages,
  gender = ifelse(gender == "M", 1, 0),
  surgery = ifelse(surgery_type == "RYGB", 1, 0),
  BMI = BMI_T0,
  protein_intake = protein
)

## binary outcome
beta = c(2, -2, 0.1, 0.5, -1)
z = -50 + as.matrix(df) %*% beta
p = 1 / (1 + exp(-z)) 
hist(p)
EFFML = rbinom(n, 1, p)      
df$EFFML = EFFML

glm(EFFML ~ ., data = df, family = "binomial")
