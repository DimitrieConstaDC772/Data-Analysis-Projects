# Problem Set 1

# Load Libraries.
library(AER)
library(readxl)
library(xts)
library(zoo)
library(dynlm)
library(stargazer)
library(urca)

# 1.2 Time Series Plots
# Load the UK Output per worker data:
UK_output <- read_excel(path = "UKOutputperworker.xls")
Une_rateee <- read_excel(path = "NON-UKnational-unemployment.xls")
# Set the datetime format of the data frame:
UK_output$Date <- as.yearqtr(UK_output$Date, format = "%Y Q%q")
Une_rateee$Date <- as.yearqtr(Une_rateee$Date, format = "%Y Q%q")
# Convert the data to xts to make it compatible with time series analysis tools:
UK_output_xts <- xts(UK_output$Output, UK_output$Date)
Une_rate <- xts(Une_rateee$Rate, Une_rateee$Date)[-1]
# Transform the series in a growth rate:
UKGR_Output <- xts(100 * log(UK_output_xts/lag(UK_output_xts)))[-1]
# Plot the data to understand the patterns:
# Plot the output:
plot(as.zoo(UK_output_xts), col = "orange",
     lwd = 2,
     ylab = "Quarterly Output per worker",
     xlab = "Date",
     main = "U.K. Quarterly Output per worker - Whole Economy")

# Large outliers may violate the 3rd assumption of fourth moments.
# Smooth out shocks in the year 2020:

print(UKGR_Output["2019::2021"])

# Calculate the average growth rate for 1 year prior and after the shocks:

average <- mean(c(
  coredata(window(UKGR_Output, start = as.yearqtr("2019 Q2"), end = as.yearqtr("2020 Q1"))),
  coredata(window(UKGR_Output, start = as.yearqtr("2020 Q4"), end = as.yearqtr("2021 Q4")))
))

average

# Now we want to replace the 2 quarters with this average:
UKGR_Output[as.yearqtr("2020 Q2")] <- average
UKGR_Output[as.yearqtr("2020 Q3")] <- average

# Plot the Output per worker growth rate and the second order lag:
UKGR_Output_lag1 <- lag(UKGR_Output, 2)
par(mar=c(4.3, 4.3, 7, 2), xpd=TRUE)
plot(as.zoo(UKGR_Output), col = "blue",
     lwd = 4,
     ylab = "Quarterly GR Output/worker",
     xlab = "Date",
     mgp = c(3, 1, 0),
     mar= c(1,1,1,1), 
     width=5,
     height=5,
     main = "U.K. Quarterly Output Growth Rate per worker after smoothing - Whole Economy")

lines(as.zoo(UKGR_Output_lag1),
      type = "l", lwd = 2,
      col = "red")
legend("top",
       c("GR Output per worker", "GR Output per worker, lag 2"),
       col=c("blue", "red"),
       lty=c(1,1),
       inset=c(0,-0.25),
       cex=0.8)

# 2 Autoregression Analysis.

# 2.1 Estimate an autoregression model:
# Run an AR(1) of the Output per worker:
UK_output_AR1 <- lm(UKGR_Output ~ lag(UKGR_Output,1))
coeftest(UK_output_AR1)

# Get stargazer table:
stargazer(UK_output_AR1,
          digits = 3,
          header = F,
          type = "html",
          title = "AR(1) Regression of Output Growth",
          out = "AR(1)Regression")

# Run a Bayes Information Criterion of lag 4:
# Define a function that computes the BIC for each lag:
BIC <- function(model) {
  ssr <- sum(model$residuals^2) # Sum of squared residuals.
  t <- length(model$residuals) # The length of residuals needed, this is the same as the number of observations.
  p <- length(model$coef) - 1 # Subtract 1 from the list of coefficients to exclude the intercept.
  return(
    round(c("p" = p,
            "BIC" = log(ssr/t) + ((p+1)*log(t)/t)),
          4)
  )
}

# Loop through the function for different lags:
for (p in 1:4) {
  print(BIC(lm(UKGR_Output ~ lag(UKGR_Output,1:p))))
  }

# Estimate the model with the lowest BIC:
# Estimate AR(2) as well for comparison:
UK_output_AR2 <- lm(UKGR_Output ~ lag(UKGR_Output,1:2))
coeftest(UK_output_AR2)

# Get stargazer for both models:
stargazer(UK_output_AR1, UK_output_AR2,
          digits = 3,
          header = F,
          type = "html",
          title = "AR(1) and AR(2) Regression of Output Growth",
          out = "AR(1,2)Regression")

# Test for Key Time Series Assumptions:
# Conduct an Augmented Dickey Fuller test for AR(2) to test for a unit root process:
UK_output_AR1_URtest <- ur.df(UKGR_Output, lags = 2, type = "drift") # Include drift.
ADF_AR1_valuetest <- UK_output_AR1_URtest@teststat[1]
ADF_AR1_critvalues <- UK_output_AR1_URtest@cval[1,]

print(c("ADF test statistic:", round(ADF_AR1_valuetest,2)))
print("ADF Critical values:")
print(ADF_AR1_critvalues)

# Add a trend to the model as indicated by economic theory to see if it improves the target estimate:
UK_output_AR1_URtest_T <- ur.df(UKGR_Output, lags = 2, type = "trend")
ADFT_AR1_valuetest <- UK_output_AR1_URtest_T@teststat[1]
ADFT_AR1_critvalues <- UK_output_AR1_URtest_T@cval[1,]

print(c("ADF with trend test statistic:", round(ADFT_AR1_valuetest,2)))
print("ADF with trend Critical values:")
print(ADFT_AR1_critvalues)


# QLR test:
num_periods <- length(UKGR_Output)
tau_zero = round(0.25*num_periods,digits=0)
tau_one = round((1-0.05)*num_periods,digits=0)
n_tests <- tau_one - tau_zero + 1
tau <- seq(tau_zero, tau_one) (Exter, line)

# Set a binary variable and run a Chow test:
D <- 1*(time(UKGR_Output) > time(UKGR_Output)[tau[1]]) 
Chow_test = lm(UKGR_Output ~ lag(UKGR_Output,1) + D + (D*lag(UKGR_Output,1)))
coeftest(Chow_test)

# Array of chow test stats and run a for loop to compute a chow test for each tau:
chow_test_stat <- array(n_tests)
for (i in 1:n_tests) {
  D <- 1*(time(UKGR_Output) > time(UKGR_Output)[tau[i]])
  chow_test_AR1 = lm(UKGR_Output ~ lag(UKGR_Output,1) + D + (D*lag(UKGR_Output,1)))
  chow_test_hyp = linearHypothesis(chow_test_AR1, c("D=0", "lag(UKGR_Output, 1):D=0"), test="F", white.adjust = FALSE)
  chow_test_stat[i] = chow_test_hyp$F[2]
}
data.frame("Level" = tau, "F-stat" = chow_test_stat)

# Get the max F-stat:
QLR_test = max(chow_test_stat)
cat("QLR Test Statistic: ", QLR_test)
tau_est <- tau[which.max(chow_test_stat)]
UKGR_Output[tau_est]

D <- 1*(time(UKGR_Output) > time(UKGR_Output)[tau_est])
chow_test_AR = lm(UKGR_Output ~ lag(UKGR_Output,1) + D + (D*lag(UKGR_Output,1)) + lag(UKGR_Output,2))
coeftest(chow_test_AR)

# Finally, we can see what this looks like by plotting our data with a line at
# break period to see if the results look like what we see on the plot.
plot(as.zoo(UKGR_Output), col = "blue",
     lwd = 1,
     ylab = "Quarterly Output Growth rate per worker",
     xlab = "Date",
     main = "U.K. Quarterly Output Growth rate brake point")
abline(v=as.yearqtr("2020 Q1"), col = "red", lwd = 2)


# Get coefficient table for UR test with drift and trend:
UKGR_URtest <- lm(diff(UKGR_Output) ~ lag(UKGR_Output,1) + lag(diff(UKGR_Output),1:2)) 
coeftest(UKGR_URtest)
UKGR_URtestT <- lm(diff(UKGR_Output) ~ time(UKGR_Output) + lag(UKGR_Output,1) + lag(diff(UKGR_Output),1:2))
coeftest(UKGR_URtestT)


# Show the coefficient table:
stargazer(UK_output_AR1, UK_output_AR2, UKGR_URtest, chow_test_AR,
          digits = 3,
          header = F,
          type = "html",
          title = "AR (1), AR (2), AR (2) with drift, and AR (2) with brake Coefficient Table",
          out = "AR1 Complete Coeff table")

stargazer(UKGR_URtestT,
          digits = 3,
          header = F,
          type = "html",
          title = "ADF test with trend Coefficient Table",
          out = "ADF test with trend coeftable")
# 2.2 Estimate an Autoregressive Distributed Lag Model

# Estimate an ADl(1,1) model:
ADL_OutRate <- lm(UKGR_Output ~ lag(UKGR_Output,1) + D + (D*lag(UKGR_Output,1)) + lag(Une_rate,1) + (D*lag(Une_rate,1)))
coeftest(ADL_OutRate)

# Estimate ADL(2,2):
ADL2_OutRate <- lm(UKGR_Output ~ lag(UKGR_Output,1) + lag(UKGR_Output,2) + D + (D*lag(UKGR_Output,1)) + lag(Une_rate,1) + lag(Une_rate,2) + (D*lag(Une_rate,1)))
coeftest(ADL2_OutRate)


# Use BIC to select p length.

# Run a Bayes Information Criterion of lag 4:
# Define a function that computes the BIC for each lag:
BIC <- function(model) {
  ssr <- sum(model$residuals^2) # Sum of squared residuals.
  t <- length(model$residuals) # The length of residuals needed, this is the same as the number of observations.
  p <- length(model$coef) - 1 # Subtract 1 from the list of coefficients to exclude the intercept.
  return(
    round(c("p" = p,
            "BIC" = log(ssr/t) + ((p+1)*log(t)/t)),
          4)
  )
}

# Loop through the function for different lags:
for (p in 1:4) {
  print(BIC(lm(UKGR_Output ~ lag(UKGR_Output,1:p) + D + (D*lag(UKGR_Output,1:p)) + lag(Une_rate,1:p) + (D*lag(Une_rate,1:p)))))
}

# Conduct Grange Caussality test for ADL(1,1):
linearHypothesis(ADL_OutRate, c("lag(Une_rate, 1)=0", "D:lag(Une_rate, 1)=0"), test ="F", white.adjust=FALSE)
# Conduct Grange Causality test for ADL(2,2):
linearHypothesis(ADL2_OutRate, c("lag(Une_rate, 1)=0", "lag(Une_rate, 2)=0", "D:lag(Une_rate, 1)=0"), test ="F", white.adjust=FALSE)

# ADL models coefficient table:
stargazer(ADL_OutRate, ADL2_OutRate,
          digits = 3,
          header = F,
          type = "html",
          title = "ADL (1,1) and ADL (2,2) Coefficient Table",
          out = "ADL coef table")

length(UKGR_Output)*0.25
# 109 observations, 25% is 27.25 observations. Increase this to 30 observations for enough forecasts.
# 4 quarters in a year. 30/4 = 7.5 years. Select the last 7.5 years from the sample.
length(window(UKGR_Output, start = as.yearqtr("2017 Q1"), end = as.yearqtr("2024 Q2")))

# 2.3 Out-of-sample forecast performance

# check out-of-sample forecast performance:
sample_end <- seq(2017.00, 2024.25, 0.25)
p <- length(sample_end)

forecasts <- array(c(0),dim = c(p))
true_outputs <- array(c(0),dim = c(p))
fore_errors <- array(c(0),dim = c(p))
ser <- array(c(0),dim = c(p))

for (i in 1:p){
  Sample_back = as.yearqtr(sample_end[i])
  UKGR_Output_poos = UKGR_Output[index(UKGR_Output) < Sample_back]
  Une_rate_poos = Une_rate[index(Une_rate) < Sample_back]
  # Model the data before 2017:
  ADL1_model_poos <- lm(UKGR_Output_poos ~ lag(UKGR_Output_poos,1) + lag(Une_rate_poos,1))
  # Extract information:
  ser[i] <- summary(ADL1_model_poos)$sigma
  beta0_hat = ADL1_model_poos$coefficients[1]
  beta1_hat = ADL1_model_poos$coefficients[2]
  delta1_hat = ADL1_model_poos$coefficients[3]
  true_output <-  UKGR_Output[Sample_back]
  forecast <- (beta0_hat + (beta1_hat %*% UKGR_Output[Sample_back - 0.25]) + (delta1_hat %*% Une_rate[Sample_back - 0.25]))
  fore_error <- true_output - forecast
  true_outputs[i] <- true_output
  forecasts [i] <- forecast
  fore_errors [i] <- fore_error
}


# Analysis of POOS:
ser_withins <- ser[1]
estimated_RMSFE <- sd(fore_errors)

# Print answer:
cat("Within-Sample error:", ser_withins, "\n")
cat ("Estimated Forecasted Error:", estimated_RMSFE, "\n")

# Test the significance:
t.test(fore_errors)

# Plot the difference:
true_output_plt <- xts(true_outputs, as.yearqtr(sample_end))
forecasts_plt <- xts(forecasts, as.yearqtr(sample_end))

plot(as.zoo(true_output_plt), 
     col = "purple",
     lwd = 4,
     ylab = "Percent",
     main = "AR (1,1) Pseudo Out-of-Sample Forecats of Output per worker Growth Rate")
# Add the series of pseudo-out-of-sample forecasts:
lines(as.zoo(forecasts_plt),
      lwd = 4,
      lty = 2)
# Shade area between curves (the pseudo forecast error):
polygon(x= c(time(true_output_plt), rev(time(forecasts_plt))),
        y= c(true_outputs, rev(forecasts)),
        col = "grey85",
        border = NA)
# Add a legend:
legend("bottomright",
       lty = c(1, 2, 1),
       lwd = c(2, 2, 10),
       inset=c(0,0),
       cex = 0.8,
       col = c("purple", "black", "grey85"),
       legend = c("Actual GDP growth rate",
                  "Forecasted GDP growth rate",
                  "Pseudo forecast Error"))

# Draft for the model with brake, I could not find a way to iterate through it with a brake.
# The best solution is included in this analysis.
# The issue is caused by multicollinearity. 
# Dd (the dummy variable) makes the interaction coefficients multicollinear, producing NAN. 
 
sample_end <- seq(2017.00, 2024.25, 0.25)
p <- length(sample_end)

forecasts <- array(c(0),dim = c(p))
true_outputs <- array(c(0),dim = c(p))
fore_errors <- array(c(0),dim = c(p))
ser <- array(c(0),dim = c(p))


for (i in 1:p){
  Sample_back = as.yearqtr(sample_end1[i])
  UKGR_Output_poos = UKGR_Output[index(UKGR_Output) < Sample_back]
  Une_rate_poos = Une_rate[index(Une_rate) < Sample_back]
  Dd <- 1*(index(UKGR_Output_poos) >= as.yearqtr("2020 Q1"))
  # Model the data before 2017:
  ADL1_model_poos <- lm(UKGR_Output_poos ~ lag(UKGR_Output_poos,1) + lag(Une_rate_poos,1) + Dd + (Dd*lag(UKGR_Output_poos,1)) + (Dd*lag(Une_rate_poos,1)))
  # Extract information:
  ser[i] <- summary(ADL1_model_poos)$sigma
  beta0_hat = ADL1_model_poos$coefficients[1]
  beta1_hat = ADL1_model_poos$coefficients[2]
  delta1_hat = ADL1_model_poos$coefficients[3]
  gamma0_hat = ADL1_model_poos$coefficients[4]
  gamma1_hat = ADL1_model_poos$coefficients[5]
  gamma2_hat = ADL1_model_poos$coefficients[6]
  true_output <-  UKGR_Output[Sample_back]
  forecast <- (beta0_hat + (beta1_hat %*% UKGR_Output[Sample_back - 0.25]) + 
                 (delta1_hat %*% Une_rate[Sample_back - 0.25]) + gamma0_hat + (gamma1_hat %*% (Dd * UKGR_Output[Sample_back - 0.25])) + 
                 (gamma2_hat %*% (D * Une_rate[Sample_back - 0.25])))
  fore_error <- true_output - forecast
  true_outputs[i] <- true_output
  forecasts [i] <- forecast
  fore_errors [i] <- fore_error
}


