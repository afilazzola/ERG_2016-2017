## Calculate PET

install.packages("Evapotranspiration")


## load library and contants
library(Evapotranspiration)


data("processeddata")
data("constants")

## sample
results <- ET.McGuinnessBordne(data, constants, ts="daily", message="yes")


## Data real
data.station <- read.csv("Data//ERG.climatedata.csv")

date <- as.Date(paste(data.station$year,data.station$month, data.station$days, sep="-"), "%Y-%m-%d")
Tmax <- data.frame(Data=data.station$max.temp, index=date)
Tmin <- data.frame(Data=data.station$min.temp, index=date)


data2 <- list(data.station$max.temp, data.station$min.temp)
names(data2) <- c("Tmax","Tmin")

results <- ET.McGuinnessBordne(data2, constants, ts="daily", message="yes")