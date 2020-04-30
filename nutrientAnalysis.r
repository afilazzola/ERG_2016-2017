## load packages
library(tidyverse)
library(effects)
library(lme4)
library(lmerTest)

## set colours for plots
obcol <- c("#E69F00","#56B4E9") ## Orange and blue
bgcol <- c("#000000","#707070") ## black and grey
scol <- obcol[1]
ocol <- obcol[2]

## Load functions

se <- function(x, ...) sd(x)/sqrt(length(x)) ## standard error

source("functions.r")
source("rsquared.r")scol <- obcol[1]
ocol <- obcol[2]

## Load functions

se <- function(x, ...) sd(x)/sqrt(length(x)) ## standard error

source("functions.r")
source("rsquared.r")

## inverse hyperbolic sine transformation for zero laden data that fits log transformations
##Zhang, M., Fortney, J. C., Tilford, J. M., & Rost, K. M. (2000). An application of the inverse hyperbolic sine transformation—a note. Health Services and Outcomes Research Methodology, 1(2), 165-171.
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}


## Load Aridity gradient
source("continentality.r") ## arid
arid$Year <- as.factor(arid$Year)

### abiotic characteristics
nutrient <- read.csv("Data//ERG.soilnutrients.csv")

nutArid <- merge(nutrient, subset(arid, Year==2017), by=c("Site"))

nutMean <- nutArid %>% group_by(aridity, microsite) %>%  summarize(n=mean(N),p=mean(P),k=mean(K), n.se=se(N), p.se=se(P), k.se=se(K))

m1 <- lm(log(N) ~ microsite * aridity, data=nutArid)
anova(m1)
shapiro.test(m1$residuals)
r.squared(m1)

m2 <- lm(log(P) ~ microsite * poly(aridity,2), data=nutArid)
anova(m2)
shapiro.test(m2$residuals)
r.squared(m2)

m3 <- lm(log(K) ~ microsite * aridity, data=nutArid)
anova(m3)
shapiro.test(m3$residuals)
r.squared(m3)


plot1 <- ggplot(data=nutArid, aes(x=aridity, y= log(N))) + theme_Publication() + 
  ylab("Nitrogen")+
  geom_jitter(size=2, width = 0.2, alpha=1, aes(color=microsite)) + 
  geom_smooth(method="lm", lwd=2, lty=1, color="black") +
  scale_color_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

plot2 <- ggplot(data=nutArid, aes(x=aridity, y= log(P), color=microsite, fill=microsite)) + theme_Publication() +
  ylab("Phosphorus")+
  geom_jitter(size=2, width = 0.2, alpha=1) + 
  geom_smooth(method="lm", lwd=2, formula = y ~ poly(x,2)) +
  scale_color_manual(values=c(scol, ocol)) +
  scale_fill_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

plot3 <- ggplot(data=nutArid, aes(x=aridity, y= log(K))) + theme_Publication() + 
  ylab("Potassium")+
  geom_jitter(size=2, width = 0.2, alpha=1, aes(color=microsite)) + 
  geom_smooth(method="lm", lwd=2, lty=1, color="black") +
  scale_color_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

require(gridExtra)
grid.arrange(plot1, plot2, plot3, ncol=3)


## inverse hyperbolic sine transformation for zero laden data that fits log transformations
##Zhang, M., Fortney, J. C., Tilford, J. M., & Rost, K. M. (2000). An application of the inverse hyperbolic sine transformation—a note. Health Services and Outcomes Research Methodology, 1(2), 165-171.
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}


## Load Aridity gradient
source("continentality.r") ## arid
arid$Year <- as.factor(arid$Year)

### abiotic characteristics
nutrient <- read.csv("Data//ERG.soilnutrients.csv")

nutArid <- merge(nutrient, subset(arid, Year==2017), by=c("Site"))

nutMean <- nutArid %>% group_by(aridity, microsite) %>%  summarize(n=mean(N),p=mean(P),k=mean(K), n.se=se(N), p.se=se(P), k.se=se(K))

m1 <- lm(log(N) ~ microsite * aridity, data=nutArid)
anova(m1)
shapiro.test(m1$residuals)
r.squared(m1)

m2 <- lm(log(P) ~ microsite * poly(aridity,2), data=nutArid)
anova(m2)
shapiro.test(m2$residuals)
r.squared(m2)

m3 <- lm(log(K) ~ microsite * aridity, data=nutArid)
anova(m3)
shapiro.test(m3$residuals)
r.squared(m3)


plot1 <- ggplot(data=nutArid, aes(x=aridity, y= log(N))) + theme_Publication() + 
  ylab("nitrogen")+
  geom_jitter(size=2, width = 0.2, alpha=1, aes(color=microsite)) + 
  geom_smooth(method="lm", lwd=2, lty=1, color="black") +
  scale_color_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

plot2 <- ggplot(data=nutArid, aes(x=aridity, y= log(P), color=microsite)) + theme_Publication() +
  ylab("phosphorus")+
  geom_jitter(size=2, width = 0.2, alpha=1) + 
  geom_smooth(method="lm", lwd=2, formula = y ~ poly(x,2)) +
  scale_color_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

plot3 <- ggplot(data=nutArid, aes(x=aridity, y= log(K))) + theme_Publication() + 
  ylab("potassium")+
  geom_jitter(size=2, width = 0.2, alpha=1, aes(color=microsite)) + 
  geom_smooth(method="lm", lwd=2, lty=1, color="black") +
  scale_color_manual(values=c(scol, ocol)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

require(gridExtra)
grid.arrange(plot1, plot2, plot3, ncol=3)
