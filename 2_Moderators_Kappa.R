rm(list = ls())                              # Clean the Global Environment
if (is.null(dev.list()) == FALSE) dev.off()  # Clean Plots
cat ("\014")                                 # Clean the R console

df <- read.csv("moderators.csv", 
                 header = TRUE)
library(vcd)
x <- ftable(df$Student_1_Exposure, df$Student_2_Exposure)
x
categories1 <- c("IN", "ST", "LT")
categories2 <- c("IN", "ST", "LT")
dimnames(x) <- list(rater1 = categories2, rater2 = categories1)
y <-Kappa(x)
y

x <- ftable(df$Student_1_Stress, df$Student_2_Stress)
x
categories1 <- c("STS", "BO", "VT")
categories2 <- c("STS", "BO", "VT")
dimnames(x) <- list(rater1 = categories2, rater2 = categories1)
y <-Kappa(x)
y