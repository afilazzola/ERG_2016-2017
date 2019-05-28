## load packages
library(tidyverse)
library(effects)
library(lme4)
library(lmerTest)

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

## load data
nutrient <- read.csv("Data//ERG.soilnutrients.csv")
community <- read.csv("Data//ERG.communitydata.csv")
community[is.na(community)] <- 0
census <-  read.csv("Data//ERG.phytometer.census.csv")
HOBOdata <- read.csv("Data//ERG.logger.data.csv")
sitevars <- read.csv("Data//ERG.shrub.csv")


### Positive interactions respond non-linearly to a gradient of aridity driven by winter rainfall

### Community First
commArid <- merge(community, arid, by=c("Site","Year"))
commArid[,"Year"] <- as.factor(commArid$Year)


## biomass model
m1 <- lmer(log(Biomass) ~ aridity * Microsite + (1|Year), data=subset(commArid, Biomass>0))
anova(m1, test="Chisq")
shapiro.test(residuals(m1))
r.squared(m1)

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m1, xlevels=12)


## Plot Biomass
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_jitter(data=subset(commArid, Biomass>0), aes(x=aridity, y=log(Biomass), color=Microsite), size=2, width = 0.2, alpha=0.4) + theme_Publication() + ylab("biomass")+
   geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


## species richness model
m2 <- glmer.nb(Richness ~ poly(aridity,2) * Microsite + (1|Year), data=commArid)
car::Anova(m2, type=3)


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m2, xlevels=12)


## Plot richness
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("Species richness")+
   geom_jitter(data=subset(commArid, Biomass>0), aes(x=aridity, y=Richness, color=Microsite), size=2, width = 0.2, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


## phylogenetic similarity

source("phyloAnalysis.r") ## load phylogenetic data

## join with aridity
mpd.arid <- merge(mpd.rich, arid, by=c("Site","Year"))


m3 <- lmer(mpd ~ aridity * Microsite + (1|Year), data=mpd.arid)
anova(m3, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m3, xlevels=12)


## Plot phylogenetics
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_point(data=mpd.arid, aes(x=aridity, y=mpd, color=Microsite), size=4) + theme_Publication() + ylab("biomass (IHS transformed)")+
  geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

## native vs non-native
status <- read.csv("Data//ERG.specieslist.csv")


statusComm <- community %>%  gather(Species.shorthand, abundance, 13:53) %>% merge(. , status, by="Species.shorthand") %>% 
  group_by(Year, Site, Microsite, status, ID) %>%  summarize(abd= sum(abundance))  ## plants per plot

statusComm <- merge(data.frame(statusComm), commArid, by=c("Site","Year","Microsite"))

## native only model
m4 <- glmer( abd ~ poly(aridity,2) * Microsite  + (1|Year), data=subset(statusComm, status=="native"))
car::Anova(m4, test="Chisq")

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite", "status"), m4, xlevels=12)


## abundance of natives
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit))+ theme_Publication() + ylab("plant abundance")+
  geom_point(data=statusComm, aes(x=aridity, y=ihs(abd.mean), color=Microsite), size=4) + 
  geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

## non-native only model
m5 <- glmer.nb(abd.mean ~ aridity * Microsite  + (1|Year), data=subset(statusComm, status=="non.native"))
car::Anova(m5, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity"), m5, xlevels=12)


## abundance of natives
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit))+ geom_point(data=subset(statusComm, status=="non.native"), aes(x=aridity, y=ihs(abd.mean), color=Microsite), size=4) + theme_Publication() + ylab("plant abundance")+
  geom_line() +   geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


##### analyses of phytometers

census <-  read.csv("Data//ERG.phytometer.census.csv")
end <- subset(census, Census=="end")

phyto <- merge(end, arid, by=c("Year","Site"))

## create occurrence columns
phyto[,"pha.occ"] <- ifelse(phyto$Phacelia>0, 1,0)
phyto[,"pla.occ"] <- ifelse(phyto$Plantago>0, 1,0)
phyto[,"sal.occ"] <- ifelse(phyto$Salvia>0, 1,0)

## drop panoche 2017 because too cold
phyto <- phyto %>% filter(Site!="PanocheHills" | Year!="2017") 

## Phacelia occurrence
m1.occ <- glmer(pha.occ ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=phyto, family="binomial")
car::Anova(m1.occ, test="Chisq")
## Phacelia Biomass
m1.bio <- lmer(log(Phacelia.biomass) ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=subset(phyto, pha.occ==1))
anova(m1.bio, test="Chisq")

## Plantago occurrence
m2.occ <- glmer(pla.occ ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=phyto, family="binomial")
car::Anova(m2.occ, test="Chisq")
## Plantago Biomass
m2.bio <- lmer(log(Plantago.biomass) ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=subset(phyto, pla.occ==1))
anova(m2.bio, test="Chisq")

## Salvia occurrence
m3.occ <- glmer(sal.occ ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=phyto, family="binomial")
car::Anova(m3.occ, test="Chisq")
## Plantago Biomass
m3.bio <- lmer(log(Salvia.biomass) ~ poly(aridity,2) * Microsite * Nutrient + (1|Year), data=subset(phyto, sal.occ==1))
anova(m3.bio, test="Chisq")




