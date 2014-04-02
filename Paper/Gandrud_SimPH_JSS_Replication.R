#############
# Replication file for: simPH: An R package for showing estimates for 
# interactive and nonlinear effects from Cox proportional hazard models
# Requires R 3.0.3 or greater and simPH version 1.2 or greater
# Updated 2 April 2014
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")

#### Illustration of linear effects####
# Load hmohiv data from UCLA repository
hmohiv <- read.table(
           "http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", 
           sep = ",", header = TRUE)

# Center age at its median (35)
hmohiv$AgeMed <- hmohiv$age - 35

M1 <- coxph(Surv(time, censor) ~ AgeMed + drug, 
            method = "breslow", data = hmohiv)

# Simulate relative hazards
Sim1 <- coxsimLinear(M1, b = "AgeMed", Xj = seq(-15, 19, by = 0.2))

# Plot results
simGG(Sim1, xlab = "\nYears of Age from the Sample Median (35)",
      ylab = "Relative Hazard with Comparison\n to a 35 Year Old\n",
      alpha = 0.05, type = 'lines')

# Simulate and plot binary drug relative hazard estiamates
Sim2 <- coxsimLinear(M1, b = "drug", Xj = 0:1)

simGG(Sim2, psize = 3, xlab = "",
      ylab = "Relative Hazard\n",
      type = 'points', smoother = 'lm') + 
    scale_x_continuous(breaks = c(0, 1), 
                       labels = c('\nNo Drug Use', '\nDrug Use'))

##### Illustration of time-varying interactive effects ######
# Load Golub & Steunenberg (2007) data. The data is included with simPH.
data("GolubEUPData")

# Examine the data's format
head(GolubEUPData[, 2:5])

# Expand data into equally spaced time intervals
 GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                      Time = 'begin', Time2 = 'end', event = 'event') 

# Examine the data's new format
head(GolubEUPData[, 1:4])

# Create natural log-time interactions
Golubtvc <- function(x){
  tvc(data = GolubEUPData, b = x, tvar = "end", tfun = "log")
}
GolubEUPData$Lqmv <- Golubtvc("qmv")
GolubEUPData$Lbacklog <- Golubtvc("backlog")
GolubEUPData$Lcoop <- Golubtvc("coop")
GolubEUPData$Lcodec <- Golubtvc("codec")
GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
GolubEUPData$Lthatcher <- Golubtvc("thatcher")

# Estimate model
M2 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea + qmvpostteu +
              coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher + 
              agenda + backlog + Lqmv + Lqmvpostsea + Lcoop + Lcodec +
              Lthatcher + Lbacklog,
            data = GolubEUPData, ties = "efron")

## Create simtvc object for first difference (central interval)
Sim3.1 <- coxsimtvc(obj = M2, b = "qmv", btvc = "Lqmv",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000,
                    by = 5, ci = 0.95)

# Create simtvc object for first difference (SPIn)
Sim3.2 <- coxsimtvc(obj = M2, b = "qmv", btvc = "Lqmv",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000,
                    by = 5, ci = 0.95, spin = TRUE)

# Create first difference plots
Plot3.1 <- simGG(Sim3.1, xlab = "\nTime in Days", 
                 title = "Central Interval\n", alpha = 0.3,
                 type = "ribbons", lsize = 0.5, legend = FALSE)

Plot3.2 <- simGG(Sim3.2, ylab = "", xlab = "\nTime in Days",
                 title = "SPIn\n", alpha = 0.3,
                 type = "ribbons", lsize = 0.5, legend = FALSE)

# Combine plots
grid.arrange(Plot3.1, Plot3.2, ncol = 2)

# Create simtvc object for relative hazard
Sim4 <- coxsimtvc(obj = M2, b = "backlog", btvc = "Lbacklog",
                  qi = "Relative Hazard", Xj = seq(40, 200, 40),
                  tfun = "log", from = 1200, to = 7000, by = 100,
                  nsim = 200)

# Create relative hazard plot
simGG(Sim4, xlab = "\nTime in Days", type = "ribbons",
      leg.name = "Backlogged \n Items")

##### Illustration of spline effects ######
# Load Carpenter (2002) data. The data is included with simPH.
data("CarpenterFdaData")

# Run basic model
# From Keele (2010) replication source code. Used to create Table 7.
M3 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
              acutediz + hosp01 + pspline(hospdisc, df = 4) + 
              pspline(hhosleng, df = 4) + mandiz01 + 
              femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + 
              wpnoavg3 + pspline(condavg3, df = 4) + 
              pspline(orderent, df = 4) + pspline(stafcder, df = 4), 
              data = CarpenterFdaData)

## Simulated Fitted Values
Sim5 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = seq(1100, 1700, by = 10), 
                     Xl = seq(1099, 1699, by = 10))

# Plot simulated values
Plot5 <- simGG(Sim5, xlab = "\n Number of FDA Drug Review Staff", 
                  title = "Central Interval\n", alpha = 0.1, 
                  type = "lines")

Plot5 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
            limits = c(0, 60))


# Simulated Fitted Values: shortest probability interval
Sim6 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = seq(1100, 1700, by = 10), 
                     Xl = seq(1099, 1699, by = 10), 
                     spin = TRUE)

# Plot simulated values
Plot6 <- simGG(Sim6, xlab = "\n Number of FDA Drug Review Staff",
                title = "SPIn\n", alpha = 0.1, type = "lines")

# Place on the same scale as the central interval figure
Plot6 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
            limits = c(0, 60))