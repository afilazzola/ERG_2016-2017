## Calculate aridity from CRU data

library(raster)

sites <- read.csv("Data//ERGsites.csv")

pet <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc")

pre <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc")

r <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = 58)
r@z

pet2016 <- stack()
for(i in 58:64){
  r <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = i)
  names(r) <- paste0("pet",as.character(r@z[[1]]))
  pet2016 <- stack(pet2016, r)
}


arid2016 <- stack()
for(i in 58:64){
  r1 <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = i)
  names(r1) <- paste0("pet",as.character(r1@z[[1]]))
  r2 <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc", band = i)
  names(r2) <- paste0("pre",as.character(r2@z[[1]]))
  arid <- r2/r1
  arid2016 <- stack(arid2016, arid)
  print(i)
}




pet2017 <- stack()
for(i in 65:76){
  r <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = i)
  names(r) <- paste0("pet",as.character(r@z[[1]]))
  pet2017 <- stack(pet2017, r)
}


pre2016 <- stack()
for(i in 53:64){
  r <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc", band = i)
  names(r) <- paste0("pre",as.character(r@z[[1]]))
  pre2016 <- stack(pre2016, r)
}

pre2017 <- stack()
for(i in 65:76){
  r <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc", band = i)
  names(r) <- paste0("pre",as.character(r@z[[1]]))
  pre2017 <- stack(pre2017, r)
}


arid2016 <- stack()
for(i in 58:64){
  r1 <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = i)
  names(r1) <- paste0("pet",as.character(r1@z[[1]]))
  r2 <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc", band = i)
  names(r2) <- paste0("pre",as.character(r2@z[[1]]))
  arid <- r2/r1
  arid2016 <- stack(arid2016, arid)
  print(i)
}

arid2017 <- stack()
for(i in 70:76){
  r1 <- raster("C:\\Users\\Fitz\\Downloads\\pet\\CRUpet.nc", band = i)
  names(r1) <- paste0("pet",as.character(r1@z[[1]]))
  r2 <- raster("C:\\Users\\Fitz\\Downloads\\pre\\CRUpre.nc", band = i)
  names(r2) <- paste0("pre",as.character(r2@z[[1]]))
  arid <- r2/r1
  arid2017 <- stack(arid2017, arid)
  print(i)
}


## spatial points
gps <- sites
crs.world <- arid2016@crs
coordinates(gps) <- ~long+lat
proj4string(gps) <- crs.world

## arid values
gps.2016 <- apply(extract(arid2016, gps), 1, mean)
gps.2017 <- apply(extract(arid2017, gps), 1, mean)


sites[,"arid2016"] <- gps.2016
sites[,"arid2017"] <- gps.2017


write.csv(sites, "Data//aridityCRU.csv", row.names = FALSE)



## arid values
gps.2016 <- apply(raster::extract(pet2016, gps), 1, mean)
gps.2017 <- apply(raster::extract(pre2016, gps), 1, mean)



sites[,"pre2016"] <- gps.2016
sites[,"pre2017"] <- gps.2017

