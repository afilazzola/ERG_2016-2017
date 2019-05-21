## load packages
library(tidyverse)
library(effects)

## Load functions

se <- function(x, ...) sd(x)/sqrt(length(x)) ## standard error

source("functions.r")

## inverse hyperbolic sine transformation for zero laden data that fits log transformations
##Zhang, M., Fortney, J. C., Tilford, J. M., & Rost, K. M. (2000). An application of the inverse hyperbolic sine transformation—a note. Health Services and Outcomes Research Methodology, 1(2), 165-171.
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}


## Load Aridity gradient
source("continentality.r") ## arid
arid$Year <- as.factor(arid$Year)

## load data
nutrient <- read.csv("Data//ERG.soilnutrients.csv")
community <- read.csv("Data//ERG.communitydata.csv")
community[is.na(community)] <- 0
census <-  read.csv("Data//ERG.phytometer.census.csv")
HOBOdata <- read.csv("Data//ERG.logger.data.csv")
sitevars <- read.csv("Data//ERG.shrub.csv")


### Positive interactions respond non-linearly to a gradient of aridity driven by winter rainfall

### Community First
comm <- community %>% group_by(Year, Site, Microsite) %>%  summarize(bio=mean(Biomass), abd=mean(Abundance), rich=mean(Richness), bio.se=se(Biomass), abd.se=se(Abundance), rich.se=se(Richness))
commArid <- merge(comm, arid, by=c("Site","Year"))
commArid[,"Year"] <- as.factor(commArid$Year)


## biomass model
m1 <- lmer(ihs(bio) ~ aridity * Microsite + (1|Year), data=commArid)
anova(m1, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m1, xlevels=12)


## Plot Biomass
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_point(data=commArid, aes(x=aridity, y=ihs(bio), color=Microsite), size=4) + theme_Publication() + ylab("biomass (IHS transformed)")+
   geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


## species richness
m2 <- lmer(rich ~ poly(aridity,2) * Microsite + (1|Year), data=commArid)
anova(m2, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m2, xlevels=12)


## Plot richness
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_point(data=commArid, aes(x=aridity, y=rich, color=Microsite), size=4) + theme_Publication() + ylab("biomass (IHS transformed)")+
  geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


## phylogenetic similarity

source("phyloAnalysis.r") ## load phylogenetic data

## join with aridity
mpd.arid <- merge(mpd.rich, commArid, by=c("Site","Year","Microsite"))


m3 <- lmer(mpd ~ aridity * Microsite + (1|Year), data=mpd.arid)
anova(m3, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m3, xlevels=12)


## Plot Biomass
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_point(data=mpd.arid, aes(x=aridity, y=mpd, color=Microsite), size=4) + theme_Publication() + ylab("biomass (IHS transformed)")+
  geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

