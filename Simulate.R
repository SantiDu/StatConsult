library(wakefield)  # library for generating random dataset

n = 1150  # sample size

# generate random age
ages = age(n, x = 20:65, prob = NULL)  # age range from 20 to 65

# generate random gender
n_male = rbinom(1, n, 0.2)
n_female = n - n_male
gender = c(rep("M", n_male), rep("F", n_female))

# generate random height
height_m = height_cm(n_male, 185, 7)
height_f = height_cm(n_female, 170, 6)
height = c(height_m, height_f) / 100   # height in meters

# generate BMI that is greater than 30
BMI_preop = rnorm(n, 43, 6)
not_obese = BMI_preop <= 30
BMI_preop[not_obese] = rnorm(sum(not_obese), 43, 1)  # everyone in the dataset is obese

# generate surgery types. SG: sleeve gastrectomy (SG). RYGB: Roux-en-Y gastric bypass
n_SG = rbinom(1, n, 0.5)
n_RYGB = n - n_SG
surgery_type = c(rep("SG", n_SG), rep("RYGB", n_RYGB))

# calculate pre operative weight
weight_preop = BMI_preop * height^2

# pre operative Fat Free Mass
fat_pct_preop = rnorm(n, 50, 5) / 100
FFM_preop = (1 - fat_pct_preop) * weight_preop

# post operative protein intake

# (post operative MVPA moderate-to-vigorous physical activity)

# post operative fat free mass

# post operative weight

# FFML
