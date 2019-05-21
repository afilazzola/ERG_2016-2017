
library(picante)
library(ape)
library(brranching)
library(taxize)

## load species list to create tree
spp.list <- read.csv("Data/ERG.specieslist.csv")

## create a combined species name column and drop erinoginum because not to species level
taxa <- paste(spp.list$Genus, spp.list$Species.name)
taxa <- taxa[taxa!="Erinoginum spp"]

## create tree of phylogeny
tree <- phylomatic(taxa=taxa, get = 'POST')

## calculate branch lengths using  Grafen's method
## Grafen, A. (1989) The phylogenetic regression. Philosophical Transactions of the Royal society of London. Series B. Biological Sciences, 326, 119–157.
tree2 <- compute.brlen(tree, method="Grafen")

## view tree and outputs to verify branches
plot(tree2)
tree2

## format community to only have site and microsite
community <- read.csv("Data//ERG.communitydata.csv")
community[is.na(community)] <- 0
comm <- community %>% group_by(Year, Site, Microsite) %>% summarise_if(is.numeric, funs(sum(., na.rm=T)))
comm <- comm[,c(1:3,13:53)] ## extract species abundances only
comm <- data.frame(comm) ## convert to dataframe
comm[,"site.micro"] <- paste(comm$Site,comm$Microsite,as.character(comm$Year)) ## create site by microsite column
rownames(comm) <- comm[,length(comm)] ## replace row names with site by microsite
comm[,c("Year","Site","Microsite","Erinoginum.spp","Boraginaceae.sp.","site.micro")] <- NULL ## drop all columns but ones for analysis

## match species name format
spp.list[,"spp.name"] <- paste(spp.list$Genus, spp.list$Species.name)
colnames(comm) <- spp.list$spp.name[match(colnames(comm), spp.list$Species.shorthand)]

## match species format
names(comm) <- gsub(" ", "_", names(comm), fixed=T) ## underscores instead of spaces
names(comm) <- gsub("__", "_", names(comm), fixed=T) ## replace double underscores with one
names(comm) <- tolower(names(comm)) ## lower case names


## mean phylogenetic distance
mpd.data <- mpd(comm, cophenetic(tree2), abundance.weighted=T)

## format dataframe
mpd.data <- data.frame(Year=c(rep(2016,14),rep(2017,14)),Microsite=c("open","shrub"),site=rownames(comm), mpd=mpd.data)
mpd.data[,"Site"] <- gsub(" open", "", mpd.data[,"site"]) 
mpd.data[,"Site"] <- gsub(" shrub", "", mpd.data[,"Site"]) 
mpd.data[,"Site"] <- gsub(" 2016", "", mpd.data[,"Site"]) 
mpd.data[,"Site"] <- gsub(" 2017", "", mpd.data[,"Site"]) 
mpd.data[,"site"] <- NULL
mpd.data[is.na(mpd.data)] <- 0 ## barstow no plants


## getspecies richness for sites
spp.rich <- community %>% group_by(Year, Site, Microsite) %>% summarize(richness=mean(Richness))
spp.rich <- data.frame(spp.rich)

## join data
mpd.rich <- merge(mpd.data, spp.rich, by=c("Site","Microsite","Year"))
