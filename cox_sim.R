# Load necessary libraries
set.seed(42)
library(survival)

# Simulate data
n <- 1000  # number of individuals
time_follow_up <- runif(n, 1, 10)  # follow-up times (1-10 years)
protein_intake <- sample(c(0, 1), n, replace = TRUE)  # protein intake (binary: 0=low, 1=high)

# Generate event (muscle mass loss)
baseline_hazard <- 0.02
hazard_ratio_protein <- 0.5  # reduces hazard by 50% for high protein intake
hazard <- baseline_hazard * exp(-protein_intake * log(hazard_ratio_protein))

# Time-to-event data with event indicator
event <- rbinom(n, 1, 1 - exp(-hazard * time_follow_up))

# Create a data frame
data <- data.frame(
  time_follow_up = time_follow_up,
  protein_intake = protein_intake,
  event = event
)

# Cox Proportional Hazards Regression (includes time)
cox_model <- coxph(Surv(time_follow_up, event) ~ protein_intake, data = data)

# Poisson regression to model rate of muscle mass loss
# Using the log of follow-up time as an offset
poisson_model <- glm(event ~ protein_intake + offset(log(time_follow_up)), family = poisson, data = data)

# Summary of Poisson regression
summary(poisson_model)

# Summary of Cox regression
summary(cox_model)





