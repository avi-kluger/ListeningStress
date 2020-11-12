rm(list = ls())                              # Clean the Global Environment
if (is.null(dev.list()) == FALSE) dev.off()  # Clean Plots
cat ("\014")                                 # Clean the R console

suppressMessages(library(metafor))
library(dplyr)
library(weights)

df      <- read.csv("SecondaryTraumaInput.csv")
dfClean <- subset(df, !is.na(r))

x <- escalc(measure="COR", ri = r, ni = N, 
            slab= Author,
            data = dfClean)
x$id <- as.factor(1:nrow(x))

x$Stress <- ifelse(is.na(x$Stress) | x$Stress == "", "Other", x$Stress)

# source("1_BuildMATable.R")
# Table1 <- buildTable ("Effect.of.exposure", line = 1, data = x)

fit <-
  rma.mv(yi,
         vi,
         random =  ~  1 |  factor(Author)/factor(Sample) / id,
         data = x,
         slab = Author)
summary(fit)
var0 <- sum(fit$sigma2)

library(ggplot2)
theme_set(theme_bw(base_size=10))

res.mv0 <- rma.mv(yi, vi,
                  random =  ~  1 | factor(Author)/factor(Sample)/id,
                  data = x,
                  control=list(optimizer="optim", optmethod="Nelder-Mead"))
print(res.mv0, digits = 2)

# Forest plot, but each line shows all the effects of a given paper
data <- 
    data.frame(
      ES = x$yi,
      SE = sqrt(x$vi),
      Type = "Study",
      Study = paste(x$Author, x$Sample, sep = ", "),
      Design = x$Design
    )

plot1 <-
  ggplot(data = data,
         aes(
           x = Study,
           y = ES,
           ymax = ES + (1.96 * SE),
           ymin = ES - (1.96 * SE),
         )) + geom_pointrange()
plot1 + coord_flip() 
plot1 <-
  ggplot(data = data[which(data$Design == "Experimental"), ],
         aes(
           x = Study,
           y = ES,
           ymax = ES + (1.96 * SE),
           ymin = ES - (1.96 * SE),
         )) + geom_pointrange()
plot1 + coord_flip() 

# Forest plot, correct plot IV and DV labels too long
data <- 
    data.frame(
      ES = x$yi,
      SE = sqrt(x$vi),
      Type = "Study",
      Study = paste(x$Author, x$IV, "with", x$DV)
    )

plot1 <-
  ggplot(data = data,
         aes(
           x = Study,
           y = ES,
           ymax = ES + (1.96 * SE),
           ymin = ES - (1.96 * SE),
         )) + geom_pointrange()
plot1 + coord_flip() 


#Calculate % true variance (I2)
W <- diag(1/x$vi)
X <- model.matrix(fit)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2<- 100 * sum(fit$sigma2) / (sum(fit$sigma2) + (fit$k-fit$p)/sum(diag(P)))
I2 <- paste0(round(I2, 1) ,"%")
print(I2)

fit <- rma.mv(yi, vi, random =  list(~ factor(id) | Author), data = x)
summary(fit)

# Tau
rd(sqrt(sum(fit$sigma2)), 2)

#Sensitivity analyses
boxplot(x$r)
boxplot(x$N)
max(x$N)
funnel(fit)

# Repeat meta-analysis after dropping an outlier in sample size
fit <-
  rma.mv(yi,
         vi,
         random =  ~  1 | Author / id,
         data = x[which(x$N < 1000), ],
         slab = Author)
summary(fit)

pdf("forestPlot")
forest(fit)
dev.off()



table(x$Exposure, useNA = "ifany")
# Test the moderating effect of exposure. Read Q only here
res.mv <- rma.mv(yi, vi, 
               mods   = ~ Exposure, 
               random =  list(~ factor(id) | Author),
               data = x,
               control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)

# Print the results for each exposure level
res.mv <- rma.mv(yi, vi, 
               mods   = ~ Exposure + 0, 
               random =  list(~ factor(id) | Author),
               data = x,
               control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)
print(res.mv, digits = 2)

# Amount of variance explained by the moderator
fit1 <- rma.mv(yi, vi, random =  ~  1 | Author/id, mods   = ~ Exposure, data = x)
summary(fit1)
var1 <- sum(fit1$sigma2)
R2 <- paste0(round((var0 - var1) / var0, 3) *100,"%")
print(R2) 

# Tau
rd(sqrt(sum(fit1$sigma2)), 2)

table(x$Stress, useNA = "ifany")
# Test the moderating effect of stress Read Q only here
res.mv <- rma.mv(yi, vi, 
                 mods   = ~ Stress, 
                 random =  list(~ factor(id) | Author),
                 data = x,
                 control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)

# Print the results for each exposure level
res.mv <- rma.mv(yi, vi, 
                 mods   = ~ Stress + 0, 
                 random =  list(~ factor(id) | Author),
                 data = x,
                 control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)
print(res.mv, digits = 2)

# Amount of variance explained by the moderator
fit1 <- rma.mv(yi, vi, random =  ~  1 | Author/id, mods   = ~ Stress, data = x)
summary(fit1)
var1 <- sum(fit1$sigma2)
R2 <- paste0(round((var0 - var1) / var0, 3) *100,"%")
print(R2) 

table(x$Stress, useNA = "ifany")

table(x$Design)
# Test the moderating effect of design Read Q only here
res.mv <- rma.mv(yi, vi, 
                 mods   = ~ Design, 
                 random =  list(~ factor(id) | Author),
                 data = x,
                 control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)

# Print the results for each exposure level
res.mv <- rma.mv(yi, vi, 
                 mods   = ~ Design + 0, 
                 random =  list(~ factor(id) | Author),
                 data = x,
                 control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(res.mv)
print(res.mv, digits = 2)

# Amount of variance explained by the moderator
fit1 <- rma.mv(yi, vi, random =  ~  1 | Author/id, mods   = ~ Design, data = x)
summary(fit1)
var1 <- sum(fit1$sigma2)
R2 <- paste0(round((var0 - var1) / var0, 3) *100,"%")
print(R2) 
