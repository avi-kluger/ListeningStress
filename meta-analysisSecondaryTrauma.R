rm(list = ls())                              # Clean the Global Environment
if (is.null(dev.list()) == FALSE) dev.off()  # Clean Plots
cat ("\014")                                 # Clean the R console

library(dplyr)

# df <- read.csv(
#   "Second Trauma Input - for meta analysis - Tzvi original coding sheet.csv", 
#   stringsAsFactors = FALSE)

df <- read.csv(
  "Second Trauma Input - for meta analysis - with Avi codes added.csv", 
  stringsAsFactors = FALSE)

dfClean <- subset(df, !is.na(r))
dfClean <- dfClean %>% filter(Effect.of.exposure == "Y")

# filter to remove non-listening samples
dfClean <- dfClean %>% filter(Sample != "Nurses")
dfClean <- dfClean %>% filter(Sample != "Doctors")
dfClean <- dfClean %>% filter(Sample != "Dentists")

dfClean$N <- as.numeric(dfClean$N)
dfClean$r <- as.numeric(dfClean$r)


suppressMessages(library(metafor))
x <- escalc(measure="COR", ri = r, ni = N, 
            slab= Author,
            data = dfClean)
x$id <- as.factor(1:nrow(x))

fit <- rma.mv(yi, vi, random =  ~  1 | Author/id, data = x,
              slab = Author)
summary(fit)
var0 <- sum(fit$sigma2)

library(ggplot2)
theme_set(theme_bw(base_size=10))

# Forest plot, but each line shows all the effects of a given paper
data <- 
    data.frame(
      ES = x$yi,
      SE = sqrt(x$vi),
      Type = "Study",
      Study = x$Author
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


pdf("forestPlot")
forest(fit)
dev.off()

# filter to keep only specific type of exposure
x$Exposure                                  <- "NA"
x$Exposure[x$Short.term.exposure == "Y"]    <- "Short"
x$Exposure[x$Long.term.exposure == "Y"]     <- "Long"
x$Exposure[x$Intensity.of.exposure == "Y"]  <- "Intense"

table(x$Exposure)

table(x$Exposure.type..AVI., x$Exposure)

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

