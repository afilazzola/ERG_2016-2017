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

## set colours for plots
obcol <- c("#E69F00","#56B4E9") ## Orange and blue
bgcol <- c("#000000","#707070") ## black and grey
scol <- obcol[1]
ocol <- obcol[2]



### Community First
commArid <- merge(community, arid, by=c("Site","Year"))
commArid[,"Year"] <- as.factor(commArid$Year)


## biomass model
m1 <- lmer(log(Biomass) ~ aridity * Microsite + Year + (1|ID), data=subset(commArid, Biomass>0))
anova(m1, test="Chisq")
shapiro.test(residuals(m1))
r.squared(m1)

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m1, xlevels=list(aridity=1:11))


## Plot Biomass
plot1 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ 
  geom_jitter(data=subset(commArid, Biomass>0), aes(x=aridity, y=log(Biomass), color=Microsite), size=2, width = 0.2, alpha=1) + theme_Publication() + ylab("log-transformed biomass")+
   geom_line(lwd=2) +   
  # geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol)) +
  annotate("text", x=1,y=4, label="a", size=8) +
  guides(color=guide_legend(override.aes=list(fill=NA)))


## species richness model
m2 <- glmer.nb(Richness ~ poly(aridity,2) * Microsite + Year + (1|ID), data=commArid)
car::Anova(m2, type=3)


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m2, xlevels=list(aridity=1:11))


## Plot richness
plot2 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("species richness")+
   geom_jitter(data=subset(commArid, Biomass>0), aes(x=aridity, y=Richness, color=Microsite), size=2, width = 0.2, alpha=1) +
  geom_line(lwd=2) +   
  # geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol))+
  annotate("text", x=1,y=6, label="b", size=8) +
  guides(color=guide_legend(override.aes=list(fill=NA)))


## phylogenetic similarity

source("phyloAnalysis.r") ## load phylogenetic data

## join with aridity
mpd.arid <- merge(mpd.rich, arid, by=c("Site","Year"))


m3 <- lmer(mpd ~ aridity * Microsite + (1|Year), data=mpd.arid)
anova(m3, test="Chisq")


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m3, list(aridity=1:11))


## Plot phylogenetics
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite))+ geom_point(data=mpd.arid, aes(x=aridity, y=mpd, color=Microsite), size=4) + theme_Publication() + ylab("biomass (IHS transformed)")+
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

## native vs non-native
status <- read.csv("Data//ERG.specieslist.csv")


statusComm <- community %>%  gather(Species.shorthand, abundance, 13:53) 
statusLong <- merge(statusComm , status, by="Species.shorthand")
statusComm <- statusLong %>%   group_by(ID, Year, Site, Microsite, status, Rep) %>%  summarize(abd= sum(abundance)) %>%  data.frame(.) ## plants per plot

statusComm <- merge(statusComm, arid, by=c("Site","Year"))


## native only model
m4 <- glmer.nb( abd ~ aridity* Microsite  + Year + (1|ID), data=subset(statusComm, status=="native"), nAGQ=0)
car::Anova(m4, test="Chisq", type=3)

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m4, xlevels=list(aridity=1:11))


## abundance of natives
plot3 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("native plant abundance")+
  geom_jitter(data=subset(statusComm, status=="native" & abd>0), aes(x=aridity, y=abd, color=Microsite), size=2, width = 0.2, alpha=1) +
  geom_line(lwd=2) +   
  # geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol))+
scale_y_continuous(trans='log2') +annotate("text", x=.8,y=65, label="c", size=8) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

## non-native only model
m5 <- glmer.nb(abd ~ poly(aridity,2) * Microsite  + Year + (1|ID), data=subset(statusComm, status=="non.native"), nAGQ=0 )
car::Anova(m5, test="Chisq", type=2)


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity","Microsite"), m5, xlevels=list(aridity=1:11))


## abundance of non-natives
plot4 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("non-native abundance")+
  geom_jitter(data=subset(statusComm, status=="non.native" & abd>0), aes(x=aridity, y=abd, color=Microsite), size=2, width = 0.2, alpha=1) +
  geom_line(lwd=2) +   
  # geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol))+
  scale_y_continuous(trans='log2') +annotate("text", x=.8,y=260, label="d", size=8) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4)

##### analyses of phytometers

census <-  read.csv("Data//ERG.phytometer.census.csv")
end <- subset(census, Census=="end")

phyto <- merge(end, arid, by=c("Year","Site"))
phyto$Year <- as.factor(phyto$Year)
phyto$Phacelia[is.na(phyto$Phacelia)] <- 0

## Run as Hurdle models

## create occurrence columns
phyto[,"pha.occ"] <- ifelse(phyto$Phacelia>0, 1,0)
phyto[,"pla.occ"] <- ifelse(phyto$Plantago>0, 1,0)
phyto[,"sal.occ"] <- ifelse(phyto$Salvia>0, 1,0)


## drop panoche 2017 because too cold
phyto <- phyto %>% filter(Site!="PanocheHills" | Year!="2017") 

## Compare global models
occLong <- gather(phyto, species, occ, 23:25)
globalOcc <- glmer(occ ~ species:poly(aridity,2) + Microsite * species + species:Nutrient+ as.factor(Year) + (1|ID), data=occLong, family="binomial" , nAGQ=0)
car::Anova(globalOcc, test="Chisq")

bioLong <- gather(phyto, species, bio, 13:15)
globalBio <- lmer(log(bio) ~ species:poly(aridity,2) + Microsite * species + species:Nutrient + as.factor(Year) + (1|ID), data=bioLong)
anova(globalBio, test="Chisq")


## Phacelia occurrence
m1.occ <- glmer(pha.occ ~ poly(aridity,2) * Microsite * Nutrient + Year + (1|ID), data=phyto, family="binomial" , nAGQ=0)
car::Anova(m1.occ, test="Chisq")
## Phacelia Biomass
m1.bio <- lmer(log(Phacelia.biomass) ~ poly(aridity,2) * Microsite * Nutrient + as.factor(Year) + (1|ID), data=subset(phyto, pha.occ==1))
anova(m1.bio, test="Chisq")

## Plantago occurrence
m2.occ <- glmer(pla.occ ~ poly(aridity,2) * Microsite * Nutrient + Year + (1|ID), data=phyto, family="binomial", nAGQ=0)
car::Anova(m2.occ, test="Chisq")
## Plantago Biomass
m2.bio <- lmer(log(Plantago.biomass) ~ poly(aridity,2) * Microsite * Nutrient + as.factor(Year) + (1|ID), data=subset(phyto, pla.occ==1))
anova(m2.bio, test="Chisq")

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity","Microsite"), m2.occ, se=T, xlevels=data.frame(aridity=seq(0,11, by=0.2), Microsite = c(rep("shrub",56),rep("open",56)))) %>%  data.frame()

## Salvia occurrence
m3.occ <- glmer(sal.occ ~ poly(aridity,2) * Microsite * Nutrient + Year + (1|ID), data=phyto, family="binomial", nAGQ=0)
car::Anova(m3.occ, test="Chisq")
## Salvia Biomass
m3.bio <- lmer(log(Salvia.biomass) ~ poly(aridity,2) * Microsite * Nutrient +  as.factor(Year) + (1|ID), data=subset(phyto, sal.occ==1))
anova(m3.bio, test="Chisq")


## plot the phytometers 

pha.ee <- Effect(c("aridity","Microsite"), m1.occ, se=T, xlevels=data.frame(aridity=seq(0,11, by=0.2), Microsite = c(rep("shrub",56),rep("open",56)))) %>%  data.frame()
pla.ee <- Effect(c("aridity","Microsite"), m2.occ, se=T, xlevels=data.frame(aridity=seq(0,11, by=0.2), Microsite = c(rep("shrub",56),rep("open",56)))) %>%  data.frame()
sal.ee <- Effect(c("aridity","Microsite"), m3.occ, se=T, xlevels=data.frame(aridity=seq(0,11, by=0.2), Microsite = c(rep("shrub",56),rep("open",56)))) %>%  data.frame()

## phytometer occurrence
plot1 <- ggplot(data=as.data.frame(pha.ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + 
  ylab(expression(italic("P. tanacetifolia")*"occurrence"))+ xlab("") +
   geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
  
plot2 <- ggplot(data=as.data.frame(pla.ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + 
  ylab(expression(italic("P. insularis")*"occurrence"))+ xlab("") +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

plot3 <- ggplot(data=as.data.frame(sal.ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + 
  ylab(expression(italic("S. columbariae")*"occurrence"))+
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))

gridExtra::grid.arrange(plot1, plot2, plot3, nrow=3)


## Salvia occurrence
m3.occ <- glmer(sal.occ ~ poly(aridity,2) * Microsite * Nutrient + as.factor(Year) + (1|ID), data=phyto, family="binomial", nAGQ=0)
car::Anova(m3.occ, test="Chisq")
## Salvia Biomass
m3.bio <- lmer(log(Salvia.biomass) ~ poly(aridity,2) * Microsite * Nutrient +  as.factor(Year) + (1|ID), data=subset(phyto, sal.occ==1))
anova(m3.bio, test="Chisq")

library(emmeans)
emmeans(m3.occ, pairwise~Microsite*Nutrient)

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity","Nutrient"), m3.bio, xlevels=list(aridity=0:11))

plot2 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Nutrient)) + theme_Publication() +ylab(expression(italic("S. columbariae")*"biomass")) +
  geom_jitter(data=phyto, aes(x=aridity, y=log(Salvia.biomass), color=Nutrient), size=2, width = 0.2, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Nutrient), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
   annotate("text", x=0,y=1, label="b", size=8)

grid.arrange(plot1, plot2, nrow=2)


### abiotic characteristics
nutrient <- read.csv("Data//ERG.soilnutrients.csv")

nutArid <- merge(nutrient, subset(arid, Year==2017), by=c("Site"))

nutMean <- nutArid %>% group_by(aridity, microsite) %>%  summarize(n=mean(N),p=mean(P),k=mean(K), n.se=se(N), p.se=se(P), k.se=se(K))

m1 <- lm(log(N) ~ microsite * aridity, data=nutArid)
anova(m1)
shapiro.test(m1$residuals)
r.squared(m1)

ggplot(data=nutArid, aes(x=aridity, y= log(N), color=microsite)) + theme_Publication() + ylab("nitrogen")+
  geom_jitter(size=2, width = 0.2, alpha=0.4) + 
  geom_smooth(method="lm", lwd=2) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) 


m2 <- lm(log(P) ~ microsite * poly(aridity,2), data=nutArid)
anova(m2)
shapiro.test(m2$residuals)
r.squared(m2)

plot3 <- ggplot(data=nutArid, aes(x=aridity, y= log(P), color=microsite)) + 
  theme_Publication() + ylab("phosphorus (ppm)")+
  geom_jitter(size=2, width = 0.2, alpha=0.4) + 
  geom_smooth(method="lm", lwd=2, formula= y ~ poly(x,2)) +
  scale_color_manual(values=c(scol,ocol)) +
  annotate("text", x=10,y=4, label="c", size=8)


m3 <- lm(log(K) ~ microsite * poly(aridity,1), data=nutArid)
anova(m3)
shapiro.test(m3$residuals)
r.squared(m3)

ggplot(data=nutArid, aes(x=aridity, y= log(K))) + 
  theme_Publication() + ylab("potassium (ppm)")+
  geom_jitter(size=2, width = 0.2, alpha=0.4, aes(color=microsite)) + 
  geom_smooth(method="lm", lwd=2, formula= y ~ poly(x,1), color="black") +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) 
+annotate("text", x=.8,y=260, label="d", size=8)


## calculate VPD and temperature variation
library(plantecophys)

## calculate VPD
HOBOdata[,"VPD"] <- RHtoVPD(HOBOdata$RH, HOBOdata$Temp, Pa = 101)

##generate seasons
season1 <- subset(HOBOdata, Year==2015 & Month > 10 | Year==2016 & Month < 5)
season2 <- subset(HOBOdata, Year==2016 & Month > 10 | Year==2017 & Month < 5)
data <- rbind(season1,season2)
data[,"season"] <- c(rep("season.1",nrow(season1)),rep("season.2",nrow(season2)))
data <- na.omit(data)

hoboArid <- merge(data, arid[,c("season","aridity","Site")], by=c("season","Site"))

dailyAvg <- function(x) { (max(x) + min(x))/2}

daily <- hoboArid %>% group_by(season, Site, Microsite, aridity, Rep, Year, Month, Day) %>% summarize(vpd=mean(VPD),tempDaily=dailyAvg(Temp),tempVar=var(Temp, na.rm=T))
siteAvg <- daily %>% group_by(season, aridity,Microsite, Site, Rep) %>% 
  summarize(VPD = mean(vpd), temp=mean(tempDaily),vartemp = mean(tempVar, na.rm=T))

siteAvg$season <- as.factor(siteAvg$season)

## generate unique plot ID
siteAvg[,"plotID"] <- paste0(siteAvg$Microsite,siteAvg$Rep, siteAvg$Site)

m1 <- lmer(VPD ~ Microsite * aridity + season +  (1|plotID), data=siteAvg)
anova(m1)
shapiro.test(residuals(m1))
r.squared(m1)

## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity"), m1, xlevels=list(aridity=1:11))


## VPD along gradient
plot1 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit)) + theme_Publication() + ylab("Vapour pressure deficit")+
  geom_jitter(data=siteAvg, aes(x=aridity, y=VPD, color=Microsite), size=2, width = 0.2, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol,ocol)) + scale_fill_manual(values=c(scol, ocol)) + 
  annotate("text", x=10,y=1.5, label="a", size=8)

## Temperature variation
m2 <- lmer(vartemp ~ aridity * Microsite + season + (1|plotID), data=siteAvg)
anova(m2)
shapiro.test(residuals(m2))
r.squared(m2)


## calcualte partial coefficents and confidence interval
ee <- Effect(c("aridity", "Microsite"), m2, xlevels=list(aridity=1:11))

## Temp Var
plot2 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("temperature variation")+
  geom_jitter(data=siteAvg, aes(x=aridity, y=vartemp, color=Microsite), size=2, width = 0.2, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol)) +ylim(0,80) +
  annotate("text", x=10,y=80, label="b", size=8)

### swc

end <- subset(census, Census=="emergence")
endArid <- merge(end, arid, by=c("Site","Year"))

m3 <- lmer(log(swc+0.1) ~ aridity * Microsite + as.factor(Year) + (1|ID), data=endArid)
anova(m3)
shapiro.test(residuals(m3))
r.squared(m3)

ee <- Effect(c("aridity","Microsite"), m3, xlevels=list(aridity=1:11))

## Soil moisture
ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab("vapour pressure deficit")+
  geom_jitter(data=endArid, aes(x=aridity, y=log(swc+1), color=Microsite), size=2, width = 0.2, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))


## soil compaction
siteArid <- merge(sitevars, subset(arid, Year==2017), by="Site")

m4 <- lm(log(Compaction+1) ~ aridity * Microsite,  data=siteArid)
anova(m4)
shapiro.test(residuals(m4))
r.squared(m4)


ee <- Effect(c("aridity","Microsite"), m4, xlevels=list(aridity=1:11))

## soil compaction
plot4 <- ggplot(data=as.data.frame(ee), aes(x=aridity, y=fit, color=Microsite)) + theme_Publication() + ylab(expression("soil compaction (kg / cm"^2*")"))+
  geom_jitter(data=siteArid, aes(x=aridity, y=log(Compaction+1), color=Microsite), size=2, width = 0.3, alpha=0.4) +
  geom_line(lwd=2) +   geom_ribbon(aes(ymin = lower, ymax=upper, fill=Microsite), alpha=0.3, color=NA) + 
  scale_color_manual(values=c(scol, ocol)) + scale_fill_manual(values=c(scol, ocol))+
  annotate("text", x=10,y=1.8, label="d", size=8)


grid.arrange(plot1, plot2, plot3, plot4, nrow=2)





### Ordination Analysis


comm.sum <- community %>%  group_by(Year, Microsite, Site) %>% summarize_if(is.numeric, funs(sum))
comm.sum <- data.frame(comm.sum)
commSum <- merge(comm.sum, arid, by=c("Year","Site"))

## transform data
comm.trans <- decostand(commSum[,13:53], "hell")

## clean data for ordination
## see distribution of spp
boxplot(comm.trans, xaxt="n")
labs <- colnames(comm.trans)
text(cex=0.8, x=1:41-1, y=-0.12, labs, xpd=TRUE, srt=45)

## remove spp with only one instance
comm.trans <- comm.trans[,!colSums(comm.trans)==apply(comm.trans, 2, max)]

commTrans <- comm.trans[rowSums(comm.trans)!=0,]
commSum <- commSum[rowSums(comm.trans)!=0,]

## CA or PCA
dca1 <- decorana(commTrans) ## length of gradient >2 & determine relative differences in community composition 



## conduct ordination
ord <- cca(commTrans ~ Microsite + aridity, Z = Year, data=commSum)
summary(ord)

RsquareAdj(ord)

## calculate priority
spp.priority <- colSums(commTrans)

## native status
colnative <- merge( data.frame(Species.shorthand = names(commTrans)), status, by="Species.shorthand", sort=FALSE)
colnative[,"color"] <- ifelse(colnative$status=="native", "darkorange3", "red3")



## plot ordination
par(mar=c(4.5,4.5,0.5,0.5))
plot(ord, type="n", xlab="CCA1", ylab="CCA2", xlim=c(-2,2))
orditorp(ord, display = "species", cex = 0.7, col = colnative[,"color"], priority=spp.priority, air=0.5)
orditorp(ord, display = "sites", cex = 1, col = "darkslateblue", air=30, pch=2)


##collect environmental variables for site
nutrients <- read.csv("Data/ERG.soilnutrients.csv")
nutrients.mean <- nutrients %>% group_by(microsite, Site) %>% summarise_all(funs(mean))
nutrients.vars <- data.frame(nutrients.mean)


## summarize daily means across sites
HOBO <- read.csv("Data/ERG.logger.data.csv")
HOBO.data <- HOBO %>%  group_by(Site, Microsite, Year, Month, Day) %>% summarize(temp.var=var(Temp),Temp=mean(Temp),RH=mean(RH))

means <- HOBO.data %>% group_by(Microsite, Site) %>% summarize(Temp.=mean(Temp, na.rm=T),rh=mean(RH, na.rm=T),temp.se=se(Temp),rh.se=se(RH))
means <- data.frame(means)

## soil moisture
smc.early <- subset(smc, Census=="emergence")
smc.avg <- smc.early %>% group_by(Microsite, Site) %>% summarize(SMC=mean(swc, na.rm=T))

envs <- data.frame(swc=smc.avg[,"SMC"],nutrients.vars[,c("N","P","K")], means[,c("Temp.","rh")])

## Check for collinearity
cor(envs)
envs[,"K"] <- NULL ## drop potassium for being correlated
envs[,"rh"] <- NULL ## drop humidity for being correlated

ord.env <- envfit(ord, envs)
plot(ord.env)

```
