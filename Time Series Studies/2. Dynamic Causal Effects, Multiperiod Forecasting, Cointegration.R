# Problem Set 2

# Load Libraries:
library(AER)
library(readxl)
library(xts)
library(zoo)
library(dynlm)
library(stargazer)
library(urca)

# Time Series Plots
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

# Large outliers may violate the 3rd assumption of fourth moments.
# Smooth out shocks in the year 2020:

print(UKGR_Output["2019::2021"])

# Calculate the average growth rate for 1 year prior and after the shocks:

average <- mean(c(
  coredata(window(UKGR_Output, start = as.yearqtr("2019 Q2"), end = as.yearqtr("2020 Q1"))),
  coredata(window(UKGR_Output, start = as.yearqtr("2020 Q4"), end = as.yearqtr("2021 Q4")))
))


# Now we want to replace the 2 quarters with this average:
UKGR_Output[as.yearqtr("2020 Q2")] <- average
UKGR_Output[as.yearqtr("2020 Q3")] <- average

# 1. Dynamic Causal Effects
# Distributed lag model with r(3).

# First run the model and then run the regression on residuals:
DLM1 <- lm(UKGR_Output ~ Une_rate + lag(Une_rate) + lag(Une_rate,2)) 
coeftest(DLM1)

# Get stargazer table:
stargazer(DLM1,
          digits = 3,
          header = F,
          type = "html",
          title = "DL (3) model of Output on Unemployment Rate",
          out = "DL(3)regression")


# Store residuals:
residuals1 <- DLM1$residuals
acf(DLM1$residuals, main = "Distributed Lag Residuals' Autocorrelation")
# Estimate phi 1 by regressing residuals with AR(1):
residuals1_reg <- lm(residuals1 ~ lag(residuals1))
coeftest(residuals1_reg)

# Get stargazer table:
stargazer(residuals1_reg,
          digits = 3,
          header = F,
          type = "html",
          title = "AR (1) model of Error term",
          out = "AR(1)error_regression")

# Store the estimated phi 1:
est_phi <- residuals1_reg$coefficients[2]
# Compute quasi-differencing variables using the estimated phi 1:
Y_qd <- UKGR_Output - (est_phi*lag(UKGR_Output))
X_qd <- Une_rate - (est_phi*lag(Une_rate))
X_qdlag1 <- lag(Une_rate) - (est_phi*lag(Une_rate,2)) 
# Run GLS model with the quasi-differenced variables:
GLSlm <- lm(Y_qd ~ X_qd + X_qdlag1)
coeftest(GLSlm)

# Save stargazer table:
stargazer(GLSlm,
          digits = 3,
          header = F,
          type = "html",
          title = "GLS model",
          out = "GLS_regression")

# Plot the autocorrelation in residuals:
acf(GLSlm$residuals, main = "GLS Residuals' Autocorrelation")

# 2. Multiperiod forecasting

# QLR test:
num_periods <- length(UKGR_Output)
tau_zero = round(0.25*num_periods,digits=0)
tau_one = round((1-0.05)*num_periods,digits=0)
n_tests <- tau_one - tau_zero + 1
tau <- seq(tau_zero, tau_one)

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

# Use the defined time when the break happens as worked out before:
ADL1_model <- lm(UKGR_Output ~ lag(UKGR_Output,1) + D + (D*lag(UKGR_Output,1)) + lag(Une_rate,1) + (D*lag(Une_rate,1)))
coeftest(ADL1_model)

# Get stargazer table:
stargazer(ADL1_model,
          digits = 3,
          header = F,
          type = "html",
          title = "ADL (1,1) model with break",
          out = "ADL11withbreak")


# 2.1 Iterated Multiperiod Forecast

# 10 periods beyond the data finish date.
# List of prediction dates for Output:
predi_start = as.yearqtr("2024 Q3")
predi_end = as.yearqtr("2026 Q4")
predi_period <- seq(predi_start, predi_end, 0.25)

f <- length(predi_period)

# Set-up vector for iterated forecasts for output rate and unemployment rate:
it_fores <- array(c(0), dim=c(f))
it_foresu <- array(c(0), dim=c(f))

# Estimating sample:
# Save coefficients from the ADL(1,1) model with break:
b0 = ADL1_model$coefficients[1]
b1 = ADL1_model$coefficients[2]
b2 = ADL1_model$coefficients[3]
b3 = ADL1_model$coefficients[4]
b4 = ADL1_model$coefficients[5]
b5 = ADL1_model$coefficients[6]

# Forecast using a for loop:
it_fore <- (b0 + ((b1+b4) %*% UKGR_Output[predi_start-0.25]) + b2 
            + ((b3+b5) %*% Une_rate[predi_start-0.25]))
it_fores[1] <- it_fore
it_foresu[1] <- (b3+b5) %*% Une_rate[predi_start-0.25]

for (i in 2:f) {
  it_fore <- (b0 + ((b1+b4)%*% it_fores[i-1]) + b2 
              + (b3+b5) %*% it_foresu[i-1])
  it_fores[i] <- it_fore
  it_foresu[i] <- (b3+b5) %*% it_foresu[i-1]
}

# Set up a data set:
it_fores_xts <- xts(it_fores, as.yearqtr(predi_period))

# 2.2 Direct Multiperiod Forecast

# Set vector:
di_fores <- array(c(0), dim=c(f))


for (i in 1:f) {
  
  ADL1_model <- lm(UKGR_Output ~ lag(UKGR_Output,(i)) + D + (D*lag(UKGR_Output,(i))) + lag(Une_rate,(i)) + (D*lag(Une_rate,(i))))
                                                             
  b0 = ADL1_model$coefficients[1]
  b1 = ADL1_model$coefficients[2]
  b2 = ADL1_model$coefficients[3]
  b3 = ADL1_model$coefficients[4]
  b4 = ADL1_model$coefficients[5]
  b5 = ADL1_model$coefficients[6]
  
  di_fores[i] <- (b0 + ((b1+b4) %*% UKGR_Output[predi_start-0.25]) + b2 
                          + ((b3+b5) %*% Une_rate[predi_start-0.25]))  
}
di_fores_xts <- xts(di_fores, as.yearqtr(predi_period))


# Plot the forecasts:
# Select the part of the series to plot:
UKGR_Output_fore <- tail(UKGR_Output,40)

# Plot the Output per worker:
plot(as.zoo(UKGR_Output_fore), 
     col = "purple",
     lwd = 4,
     xlim=c(2016, 2027),
     ylab = "Percent",
     main = "Multiperiod ADL(1,1) Forecasts with break of UK Output Growth Rate")
# Add the series of forecasts:
lines(as.zoo(di_fores_xts),
      lwd = 4,
      col= "green",
      lty = 2)
lines(as.zoo(it_fores_xts),
      lwd = 4,
      col= "blue",
      lty = 2)
legend("bottomleft",
       lty = c(1, 2, 1),
       lwd = c(2, 2, 10),
       cex = 0.6,
       col = c("purple", "green", "blue"),
       legend = c("Output GR", "Direct Forecast", "Iterated Forecast")) 


# Exploring the value of D:
ADL3model <- lm(UKGR_Output ~ lag(UKGR_Output,3) + D + (D*lag(UKGR_Output,3)) + lag(Une_rate,3) + (D*lag(Une_rate,3)))
coeftest(ADL3model)
ADL4model <- lm(UKGR_Output ~ lag(UKGR_Output,7) + D + (D*lag(UKGR_Output,7)) + lag(Une_rate,7) + (D*lag(Une_rate,7)))
coeftest(ADL4model)

stargazer(ADL4model,
          digits = 3,
          header = F,
          type = "html",
          title = "ADL (4,4) direct forecasting with break",
          out = "ADL44withbreak")

# 3. Cointegration

# Testing for cointegration in 2 stages:

# Plot the 2 variables and their difference:
Difference <- UKGR_Output - Une_rate

plot(as.zoo(UKGR_Output),
     plot.type = "single",
     lty = 2,
     lwd = 2,
     col = "blue",
     xlab = "Date",
     ylab = "Percent per quarter",
     ylim = c(-17, 17),
     main = "Output Growth Rate, Unemployment Rate and Their Difference")
lines(as.zoo(Une_rate),
      col = "orange",
      lwd = 2,
      xlab = "Date",
      ylab = "Percent per quarter",
      main = "Difference")
# Add the term spread series:
lines(as.zoo(Difference),
      col = "purple",
      lwd = 2,
      xlab = "Date",
      ylab = "Percent per quarter",
      main = "Difference")
legend("topright",
       lty = c(2, 1, 1),
       lwd = c(2, 2, 2),
       cex = 0.8,
       col = c("blue", "orange", "purple"),
       legend = c("Output Growth Rate", "Unemployment Rate", "Difference"))


# Test whether the trends have a unit root process using the GLS DF test:
GLS_DFOutput <- ur.ers(UKGR_Output, lag.max = 1, type = "DF-GLS", model = "trend")
summary(GLS_DFOutput)

GLS_DFUne_rate <- ur.ers(Une_rate, lag.max = 1, type = "DF-GLS", model="trend")
summary(GLS_DFUne_rate)

# Extract test statistic:
test_statistic <- GLS_DFOutput@teststat[1]
test_statistic2 <- GLS_DFUne_rate@teststat[1]

# Extract critical values as a numeric vector without dropping dimensions:
critical_values <- as.numeric(GLS_DFOutput@cval[1, ])  

# Combine into a data frame, ensuring all values are numeric:
dfgls_results <- data.frame(
  Statistic = c("Test Statistic Output per worker", 
                "Test Statistic Unemployment Rate", 
                "1% Critical Value", 
                "5% Critical Value", 
                "10% Critical Value"),
  Value = c(test_statistic, test_statistic2, critical_values[1], critical_values[2], critical_values[3]) # Explicit indexing
)

# Generate stargazer table in HTML:
stargazer(dfgls_results, summary = FALSE, rownames = FALSE, type = "html", 
          title = "DF-GLS Test summary for Output per Worker and Unemployment Rate", 
          out = "DFGLStest")

# There is some evidence that their relationship may exist, therefore theta being 1.
# The relationship is not strong, therefore theta will be estimated.

# Estimate theta:
coint_lm <- lm(UKGR_Output ~ Une_rate)
coeftest(coint_lm)
theta_est = coint_lm$coefficients[2]

# Use theta_est to generate z_est:
z_est = UKGR_Output - (theta_est*Une_rate)

# Run the regression to test for stationarity:
ADFtest_z_est <- ur.df(z_est, lags = 1, selectlags = "BIC", type = "none")
summary(ADFtest_z_est)

# Extract data:
test_statistic <- ADFtest_z_est@teststat[1]
critical_values <- as.numeric(ADFtest_z_est@cval[1, ])  

# Data frame formation:
dfgls_results <- data.frame(
  Statistic = c("Test Statistic estimated Z", 
                "1% Critical Value", 
                "5% Critical Value", 
                "10% Critical Value"),
  Value = c(test_statistic, critical_values[1], critical_values[2], critical_values[3]))

# Generate table:
stargazer(dfgls_results, summary = FALSE, rownames = FALSE, type = "html", 
          title = "GLS Test Summary", 
          out = "DFGLStestZ")