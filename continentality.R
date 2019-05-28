## Continentality
## Conrad

library(tidyverse)

clim <- read.csv("Data//ERG.climatedata.csv")

## generate seasons
season1 <- subset(clim, year==2015 & month > 9 | year==2016 & month < 5)
season2 <- subset(clim, year==2016 & month > 9 | year==2017 & month < 5)
clim <- rbind(season1,season2)
clim[,"season"] <- c(rep("season.1",nrow(season1)),rep("season.2",nrow(season2)))



# 
# ## Calculate mean daily range
# A.range <- clim  %>% mutate(diurnal=(min.temp+max.temp)/2)  %>% group_by(Site, season) %>% ## calculate diurnal range
#     summarize(DiA=mean(diurnal, na.rm=T)) ## average min/max per site
# 
# ## get lat lon
# gps <- read.csv("Data//ERGsites.csv")
# names(gps)[1] <- c("Site")
# 
# A.gps <- merge(A.range, gps, by="Site")
# 
# cont <- A.gps %>% mutate(k= (1.7*DiA)/sin(lat+10)-14)


## avg variables
climate <- clim %>%  group_by(Site, season) %>%  summarize(pre=sum(Precip, na.rm=T), min=mean(min.temp, na.rm=T), max= mean(max.temp, na.rm=T), RH=mean(RH.avg, na.rm=T))

arid <- climate %>%  mutate(aridity= pre/((min+max)/2+10)) %>%  data.frame(.)
arid[,"Year"] <- rep(c("2016","2017"),7)