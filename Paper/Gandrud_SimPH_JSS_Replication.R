#############
# Replication file for: simPH: An R package for showing estimates for interactive and nonlinear effects from Cox proportional hazard models
# Requires R 3.0.0 or greater
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")

##### Estimate the model for illustrating time-varying effects ######
# Load Golub & Steunenberg (2007) data. The data is included with simPH.
data("GolubEUPData")

# Create natural log-time interactions
Golubtvc <- function(x){
  tvc(data = GolubEUPData, b = x, tvar = "end", tfun = "log")
}
GolubEUPData$Lcoop <- Golubtvc("coop")
GolubEUPData$Lqmv <- Golubtvc("qmv")
GolubEUPData$Lbacklog <- Golubtvc("backlog")
GolubEUPData$Lcodec <- Golubtvc("codec")
GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
GolubEUPData$Lthatcher <- Golubtvc("thatcher")

# Estimate model
M1 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea + qmvpostteu +
              coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher + 
              agenda + backlog + Lqmv + Lqmvpostsea + Lcoop + Lcodec +
              Lthatcher + Lbacklog,
            data = GolubEUPData, ties = "efron")

# Create simtvc object for first difference (central interval)
Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
                  qi = "First Difference", Xj = 1, tfun = "log", 
                  from = 80, to = 2000, by = 5)

# Create first difference plot
simGG(Sim1, xlab = "\nTime in Days", title = "Central Interval\n", 
      ribbon = TRUE, lsize = 0.5, legend = FALSE, alpha = 0.3) 

## Create simtvc object for first difference (central interval)
Sim1.1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
                  qi = "First Difference", Xj = 1,
                  tfun = "log", from = 80, to = 2000,
                  by = 15, ci = 0.95)

# Create simtvc object for first difference (SPIn)
Sim1.2 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
                  qi = "First Difference", Xj = 1,
                  tfun = "log", from = 80, to = 2000,
                  by = 15, ci = 0.95, spin = TRUE)


# Create first difference plots
Plot1.1 <- simGG(Sim1.1, xlab = "\nTime in Days", 
                title = "Central Interval\n", alpha = 0.3,
                ribbon = TRUE, lsize = 0.5, legend = FALSE) 
Plot1.2 <- simGG(Sim1.2, ylab = "", xlab = "\nTime in Days",
                 title = "SPIn\n", alpha = 0.3,
                 ribbon = TRUE, lsize = 0.5, legend = FALSE)

# Combine plots
grid.arrange(Plot1.1, Plot1.2, ncol = 2)


##### Estimate the model for illustrating spline effects ######
# Load Carpenter (2002) data. The data is included with simPH.
data("CarpenterFdaData")

# Run basic model
# From Keele (2010) replication source code. Used to create Table 7.
M2 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
              acutediz + hosp01 + pspline(hospdisc, df = 4) + 
              pspline(hhosleng, df = 4) + mandiz01 + 
              femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + 
              wpnoavg3 + pspline(condavg3, df = 4) + 
              pspline(orderent, df = 4) + pspline(stafcder, df = 4), 
              data = CarpenterFdaData)

## Simulated Fitted Values
Sim3 <- coxsimSpline(M2, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = seq(1100, 1700, by = 10), 
                     Xl = seq(1099, 1699, by = 10))

# Plot simulated values
SimPlot1 <- simGG(Sim3, xlab = "\n Number of FDA Drug Review Staff", 
        title = "Central Interval\n", alpha = 0.2)

SimPlot1 + scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60))


# Simulated Fitted Values: shortest probability interval
Sim4 <- coxsimSpline(M2, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = seq(1100, 1700, by = 10), 
                     Xl = seq(1099, 1699, by = 10), 
                     spin = TRUE)

# Plot simulated values
SimPlot2 <- simGG(Sim4, xlab = "\n Number of FDA Drug Review Staff",
                title = "SPIn\n", alpha = 0.2)

# Place on the same scale as the central interval figure
SimPlot2 + scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60))
