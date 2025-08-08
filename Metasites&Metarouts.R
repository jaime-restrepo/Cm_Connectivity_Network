#NEEDS:
#Add code to verify that the end number of metaroutes makes sense relative to the original number of routes and the number of routes lost when sites are aggregated together in a metasite


# load required libraries 
library(tidyverse)
library(assertthat)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(geosphere)
library(readxl)


##ALL SITES HAVE BEEN AGGREGATED OUTSIDE OF R BEFORE RUNNING THIS SCRIPT.  nEW FIELDS ADDED TO THE "georeferenced_sites" FILE INCLUDE:
##MetaID, MetaSiteID, Duplicate, Rationale, and MonthsPresent.  Refer to scaling up notes for the methods used to aggregate sites.



# set wd
# ORIGINAL setwd("C:\\Dropbox\\MiCO\\NetworkModels\\")
setwd("C:/LEGOLAS/R/Lit. Rev/Aggregation")

#Specify species
# ORIGINAL spp <- "PHNI"
spp <- "CHMY"

#sitesfile <- paste("AggregatedSites\\georeferenced_sites_", spp, "_unsupervised_v1.csv", sep="")  #uses the old v1 aggregation
# ORIGINAL  sitesfile <- paste("AggregatedSites\\Testing\\georeferenced_sites_", spp, "_unsupervised_v13.4.csv", sep="") #note: I have changed the name to go with the standardised version control
sitesfile <- paste("C:/LEGOLAS/R/Lit. Rev/Aggregation/", spp, "_unsupervised_v13.4.JRpo.csv", sep="")
## ORIGINAL routesfile <- "Routes\\route_details_multiple-species_20220829.csv"
routesfile <- "route_details_Chmy_V1&2.csv"


# load georeferenced sites & routes
sites <-  read.csv(sitesfile, header = T, sep = ",", dec = ".", stringsAsFactors = F, colClasses=c(Method="character"))
#sites <- read_excel(sitesfile, sheet=NULL, range = NULL, col_names=TRUE, col_types=NULL, trim_ws=TRUE)
routes <- read.csv(routesfile, header=T, sep = ",", dec = ".", stringsAsFactors = F) 
routes <- routes %>% subset(Species == spp) #select only routes associated with this species

#Cleanup sites file
if("Commonname" %in% colnames(sites)){sites <- sites %>% rename(CommonName = "Commonname")}
species <- sites$CommonName[1]

##Clean up NumIDs and Activites fields 
sites$NumIDs <- sites$NumIDs %>% na_if("Unk") %>% as.numeric() 
sites$Activities <- str_replace_all(sites$Activities, " ", "")

##get rid of any sites that have no lats or longs
sites <- sites %>% subset(!(is.na(Lat) | is.na(Long)))

##update sites file UniqueID to include species info
sites$UniqueID <- paste(sites$Species, sites$ZoteroID, sites$SiteID, sep = "_")

#Cleanup routes file
##standardise column names
colnames(routes) <- str_replace_all(colnames(routes),c("[.]" = "", "Xanimals" = "NumIndividuals"))
routes$NumIndividuals <- routes$NumIndividuals %>% as.numeric() 

##generate ZoteroID and Unique_IDs for routes file
routes$ZoteroID <- str_extract(routes$ReviewID, "^.{8}")
routes$UniqueID_from <- paste(routes$ZoteroID, routes$SiteFrom, sep = "_")
routes$UniqueID_to <- paste(routes$ZoteroID, routes$SiteTo, sep = "_")

#Generate new unique ID info for metasites info... and redundant "to" and "from" IDs to decrease errors during joins later on
sites$MetaUniqueID <- paste(sites$MetaID, sites$MetaSiteID, sep = "_")
sites$MetaUniqueID_from <- sites$MetaUniqueID
sites$MetaUniqueID_to <- sites$MetaUniqueID


##RECALL THAT ALL SITES HAVE HAD A FIELD ADDED ("MetaSiteID") OUTSIDE OF R THAT ALLOWS THIS SCRIPT TO AGGREGATE THEM.
#Aggregate sites
metasites <- sites %>%
  group_by(MetaUniqueID) %>%
  mutate(val =  ifelse(n() > 1, toString(geomean(cbind(Long,Lat))),toString(cbind(Long,Lat))),
         MaxRadius = max(Radius, na.rm=TRUE), 
         IndividualSites = paste(unique(UniqueID), collapse=";"), 
         Methods = paste(unique(Method), collapse=";"), 
         SamplingTech = paste(unique(Sampling), collapse=";"), 
         # ORIGINAL SiteLocations = paste(unique(SiteLocation), collapse=";"),
         SiteLocations = paste(unique(Location), collapse=";"),
         SiteYears = paste(unique(na.omit(Years)), collapse=";")) %>%
  separate(val, c("MeanLon", "MeanLat"), sep = ",") %>%
  mutate_at(c("MeanLon", "MeanLat"), as.numeric)

metasites <- ungroup(metasites)


##WORK ON ROUTES
#Remove routes associated with duplicate records
metaroutes <- inner_join(routes,subset(metasites, select = c(UniqueID_from, MetaUniqueID_from)), by = "UniqueID_from")
dplyr::setdiff(routes$UniqueID_from,metaroutes$UniqueID_from)
metaroutes <- inner_join(metaroutes,subset(metasites, select = c(UniqueID_to, MetaUniqueID_to)), by = "UniqueID_to")

#Remove routes that now start and end at the same site
metaroutes <- metaroutes %>% subset(!MetaUniqueID_from == MetaUniqueID_to)

#aggregate routes between the same two nodes
metaroutes <- metaroutes %>% group_by(MetaUniqueID_from, MetaUniqueID_to) %>%
  mutate(SumNumIndividuals = sum(NumIndividuals, na.rm = TRUE)) %>%
  distinct(MetaUniqueID_from, MetaUniqueID_to, .keep_all = TRUE)


#Remove duplicate records (i.e., records that used the same data to describe the same sites)
#NOTE: we can't do this before the previous code working on routes or it eliminates routes
#that are associated with duplicate sites, even if they are not duplicate routes
metasites$Duplicate <- as.character(metasites$Duplicate)
metasites <- metasites %>% replace_na(list(Duplicate = "FALSE"))
metasites <- metasites %>% filter(!Duplicate == TRUE)


#Now that the routes have been aggregated, we can finish aggregating the sites
metasites <- metasites %>% group_by(MetaUniqueID) %>%
  mutate(MinNumIndividuals = sum(NumIDs, na.rm = TRUE), #THIS SEEMS TO BE NOT COUNTING WHEN THERE ARE NAs
         SiteType = paste0(Activities, collapse=";")) %>%
  distinct(MetaUniqueID, .keep_all = TRUE)



#Code to aggregate and carry forward ALL Activities
activitytypes <- c("SPA", "NES", "BRE", "FEE", "STA", "STO",  "MIG","WIN", "NBR", "OBS", "")
fullactivity <- c("Spawning", "Nesting", "Breeding", "Feeding/Foraging", "Staging", "Stopover",  "Migrating","Wintering", "Nonbreeding", "Observations", "")


##Clean up SiteType
metasites$SiteType <- metasites$SiteType %>%
  str_replace_all(" ","") %>%
  str_split(";")
metasites <- metasites %>% unnest(SiteType)
metasites$SiteType <- factor(metasites$SiteType, levels = activitytypes, ordered = TRUE) #%>%
metasites <- metasites %>%
  group_by(MetaUniqueID) %>%
  arrange(SiteType, .by_group = TRUE) %>%
  mutate(SiteType = paste(unique(as.character(SiteType)), collapse=";")) %>%
  distinct(MetaUniqueID, .keep_all = TRUE) %>% #order the groups by activity type so that we can pull out the first (most important activity type for use in the symbology)
  mutate(BehaviorForSymbology = str_sub(SiteType, 1, 3))


#PREPARE AND EXPORT METASITE AND METAROUTE FILES
#get a subset of fields to export
metasites_exported <- metasites %>% subset(select = c(CommonName, Species, MetaUniqueID, Methods, SamplingTech, IndividualSites, MinNumIndividuals, SiteType, BehaviorForSymbology, MeanLat, MeanLon, MaxRadius, SiteLocations, Basin)) %>% rename(SiteLocation = SiteLocations)
metaroutes_exported <- metaroutes %>% subset(select = c(MetaUniqueID_from, MetaUniqueID_to, SumNumIndividuals, Method))

## ORIGINAL spp_metasite_file <- paste(c("MetasitesAndMetaroutes\\",spp,"_metasites_testing_dcd.csv"), collapse = "")
## ORIGINAL spp_metaroute_file <- paste(c("MetasitesAndMetaroutes\\",spp,"_metaroutes_testing_dcd.csv"), collapse = "")

spp_metasite_file <- paste(c("/LEGOLAS/R/Lit. Rev/Aggregation/",spp,"_metasites_testing_dcd.csv"), collapse = "")
spp_metaroute_file <- paste(c("/LEGOLAS/R/Lit. Rev/Aggregation/",spp,"_metaroutes_testing_dcd.csv"), collapse = "")

write.csv(metasites_exported, spp_metasite_file, row.names = FALSE)
write.csv(metaroutes_exported, spp_metaroute_file, row.names = FALSE)




