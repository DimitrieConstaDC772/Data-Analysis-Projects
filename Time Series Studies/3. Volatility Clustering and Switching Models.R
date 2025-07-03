#Problem Set 3

#Load Libraries:
library(AER)
library(readxl)
library(xts)
library(zoo)
library(dynlm)
library(stargazer)
library(urca)
library(orcutt)
library(vars)
library(fGarch)
library(quantmod)
library(TSA)
library(MSwM)
library(markovchain)

# Time Series Plots
#Load the UK Output per worker data:
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
dates <- index(UKGR_Output)


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


# Plot percentage changes in output per worker:
plot(as.zoo(UKGR_Output), 
     ylab = "Percent",
     xlab = "Date",
     main = "Quarterly Percentage Changes of Output per worker in U.K.",
     type="l", 
     col = "steelblue", 
     lwd = 0.5)
# add horizontal line
abline(0, 0)

# Fit an arma(1,1) and garch(1,1) model to the Output per worker changes:
GR_Output <- garchFit(~arma(1,1) + garch(1,1), data = UKGR_Output, trace = F)
summary(GR_Output)

error_analysis <- capture.output(summary(GR_Output))
print(error_analysis)
writeLines(error_analysis, "error_analysis.txt") #(OpenAI, 2023)


# Deviations from their mean
demean_GROutput <- UKGR_Output - GR_Output@fit$coef[1]
# The Output data is demeaned to center the quarterly changes around 0.
# Mean plus one conditional standard deviation:
mean_plusSD <- GR_Output@fit$coef[1] + GR_Output@sigma.t
# Mean minus one conditional standard deviation:
mean_minusSD <- GR_Output@fit$coef[1] - GR_Output@sigma.t
# These will give a confidence interval of the expected volatility.

# Convert to xts:
demean_GROutput_xts <- xts(demean_GROutput, dates)
mean_plusSD_xts <- xts(mean_plusSD, dates)
mean_minusSD_xts <- xts(mean_minusSD, dates)

# The plotting will not be narrowed down on a specific period because the data is quarterly.

# Plot the interval:
plot(as.zoo(demean_GROutput_xts), 
     type = "l", 
     col = "steelblue",
     ylab = "Percent", 
     xlab = "Date",
     main = "ARMA(1,1)-GARCH(1,1) \n Estimated Bands of +- One Conditional Standard Deviation",
     lwd = 0.2)
# Add horizontal line at y = 0:
abline(0, 0)
# Add GARCH confidence bands (one standard deviation) to the plot:
lines(as.zoo(mean_plusSD_xts), 
      col = "darkred", 
      lwd = 0.5)
lines(as.zoo(mean_minusSD_xts), 
      col = "darkred", 
      lwd = 0.5)
legend("topright",
       lty = c(1, 1),
       lwd = c(0.5, 0.5),
       cex = 0.6,
       col = c("steelblue", "darkred"),
       legend = c("Actual Percentage Change", "+/- Conditional Standard Deviation"))

# Switching Models

# SETAR

# Change Output to xts object:
UKGR_Output_SETAR <- tar(as.vector(UKGR_Output), p1=1, p2=1, d=1, a=0.1, b=0.9, is.constant1=FALSE, is.constant2=FALSE, print=TRUE)
UKGR_Output_thd <- UKGR_Output_SETAR$thd
UKGR_Output_R1 <- UKGR_Output_SETAR$qr1$coefficients
UKGR_Output_R2 <- UKGR_Output_SETAR$qr2$coefficients

print(UKGR_Output_thd)
print(UKGR_Output_R1)
print(UKGR_Output_R2)

# Plot model:
plot(as.zoo(UKGR_Output["2019::2022"]),
     main = "SETAR(1) Model - Switching Mean",
     ylab = "Y_t",
     xlab = "t")
abline(0.33, 0, col = "darkgreen")
abline(0, 0, col = "violetred")
abline(0.29, 0, col = "darkgreen")

# Markow swithcing model:
mod_UKGR <- lm(UKGR_Output ~ lag(UKGR_Output))
mod_MKS <- msmFit(mod_UKGR, k=2, sw=c(TRUE,TRUE,TRUE))

trans_matrix <- t(mod_MKS@transMat)
coef_est <- mod_MKS@Coef

print(trans_matrix)
print(coef_est)



