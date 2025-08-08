###Major changes:
####Uses a degree buffer rather than a distance around a site to find neighbours.
####This solves the issue of weird aggregations at high latitudes.
####12.1 adds:
####   a check on radii to confirm they are standardised to 1, 5 or 10 degrees
####   a fix for how longroute sites are kept or dropped (uses a 1 degree buffer instead of a 110km distance)
####13 adds:
####   a fix that stops larger sites from being aggregated together with higher resolution reproduction sites
####   a change to how the non-reproduction sites are ordered before aggregation, so that the script aggregates sites with many neighbors first.
####13.1 adds:
####   updates to all distance measurements to ensure they are using the correct unit
####   small fixes to make sure all sf objects have a coordinate system defined
####   standardisation of code used for calculating distances
####   a fix to sites in longroutes and those sites close to a longroute site are aggregated with the lookup_site
####13.2 adds:
####   a fix to the way the longroute sites are aggregated or not by standardising around use of st_intersection with a degree buffer rather than distances
####   improves reporting on how longroute sites (and ones close by) are handle during aggregation
####   moves all longroute site reporting to Comments 
####13.3 adds:
####   random fix to the duplicate script where it was reading empty dup_test_sites dataframe as having a row when determining whether to tun the loop around line 347
####13.4 adds:
####   calculates sites that should be removed from an aggregation because they are close to a longroute_site based on buffer of longroute_site and whether distance to lookup_site is longer than distance to longroute_site

# load required libraries 
library(tidyverse)
library(assertthat)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(geosphere)
library(readxl)
library(data.table)
library(sp)
library(rgdal)
library(rgeos)
library(sf)
library(maps)


##ALL SITES HAVE BEEN AGGREGATED OUTSIDE OF R BEFORE RUNNING THIS SCRIPT.  NEW FIELDS ADDED TO THE "georeferenced_sites" FILE INCLUDE:
##MetaID, MetaSiteID, Duplicate, Rationale, and MonthsPresent.  Refer to scaling up notes for the methods used to aggregate sites.

#Function to run a conditional filter in dplyr
priority_filter <- function(df, group){
  if(group == 1)
    dplyr::filter(df, Radius == lookup_site$Radius)
  else df
}




#Specify species

spp <- "CHMY"
last_clust_num <- 1
metasites <- data.frame() #makes an empty dataframe to hold the metasites

sitesfile <- paste("AggregatedSites\\Testing\\georeferenced_sites_", spp, "_draft.csv", sep="") #normal line of code for pulling in the site file
routesfile <- "Routes\\route_details_multiple-species_20220829.csv"


# load georeferenced sites & routes

sites <-  read.csv("georeferenced_sites_CHMY_draft.csv", header = T, sep = ",", dec = ".", stringsAsFactors = F, colClasses=c(Method="character"))

routes <- read.csv("route_details_Chmy_V1&2.csv", header=T, sep = ",", dec = ".", stringsAsFactors = F)


#Cleanup sites file
if("Commonname" %in% colnames(sites)){sites <- sites %>% rename(CommonName = "Commonname")}



##Clean up NumIDs and Activities fields 
sites$NumIDs <- sites$NumIDs %>% na_if("Unk") %>% as.numeric() 
sites$Activities <- str_replace_all(sites$Activities, " ", "")
sites$Longroute <- NA
sites <- sites %>% mutate(UniqueID_method = paste(UniqueID, Method, sep = "_"))

##Get rid of any sites that have no lats or longs
sites <- sites %>% subset(!(is.na(Lat) | is.na(Long)))
sites$Activities <- sites$Activities %>% replace_na("OBS")
site_count <- nrow(sites)

#Confirm radii are standardised to 1, 5 or 10 degrees
std_radii <- c(1,5,10)
sites <- sites %>% subset(Radius %in% std_radii)
assert_that(site_count == nrow(sites), 
            msg = paste("Error: Sites were lost because their radius was not 1, 5 or 10 degrees")) #stops the script if sites were due to non-standard radius in datasheet

#Cleanup routes file
##Select only routes for this species and remove those from unused methods
routes <- routes %>% subset(Species == spp & !(Method %in% c("Isotope","Genetic","Acoustic")))

##Standardise column names
colnames(routes) <- str_replace_all(colnames(routes),c("[.]" = "", "Xanimals" = "NumIndividuals"))
routes$NumIndividuals <- routes$NumIndividuals %>% na_if("Unk") %>% as.numeric()  

##Generate ZoteroID and Unique_IDs for routes file
routes$ZoteroID <- str_extract(routes$ReviewID, "^.{8}")
routes$UniqueID_from <- paste(routes$ZoteroID, routes$SiteFrom, sep = "_")
routes$UniqueID_to <- paste(routes$ZoteroID, routes$SiteTo, sep = "_")

routes <- routes %>% left_join(unique(subset(sites, select = c(UniqueID_from, Long, Lat))), by = "UniqueID_from") %>%
  rename(FromLong = Long, FromLat = Lat)
routes <- routes %>% left_join(unique(subset(sites, select = c(UniqueID_to, Long, Lat))), by = "UniqueID_to") %>% 
  rename(ToLong = Long, ToLat = Lat)


#GENERATE LIST OF POTENTIALSITE IDs.  NOTE: MAX WITH THE CURRENT 2-LETTER ID IS 676.
site_ids = list()
for (i in 1:length(letters)) {
  for (j in 1:length(letters)) {
    k <- str_to_upper(paste(letters[i],letters[j], sep=""))
    site_ids <- append(site_ids,k)
  }
}

#Start with the sites with small nesting/breeding/calving (etc) sites and then those with large radii
#Nesting and other important sites are ordered starting with smaller radii to try to keep as many other high resolution 
#sites aggregated together as possible instead of having them get aggregated with sites with larger radii

sites <- sites %>% mutate(Important = ifelse(str_detect(Activities, c("NES|BRE|NUR|SPA|CAL|FLE|HAT")),1,2), NumCloseSites = NA)
sort_sites1 <- sites %>% filter(Important == 1) %>% arrange(Radius)
sort_sites2 <- sites %>% filter(Important == 2)

sort_sites2_sf <- sf::st_as_sf(sort_sites2, coords = c("Long","Lat")) #create an sf object of all sites

for (m in 1:nrow(sort_sites2)){ #iterate through sites to determine which have the most sites close by to aggregate those first
  next_site <- sort_sites2[m,]
  
  #Find nearby sites  
  next_site_sf <- sf::st_as_sf(next_site, coords = c("Long","Lat")) %>% st_set_crs("4326") #create sf object of the lookup_site
  next_buffer_sf <- sf::st_buffer(next_site_sf, next_site$Radius) %>% st_set_crs("4326") #Create an sf object of the buffer of the lookup_site by its radius
  
  num_close_sites <- sf::st_intersects(sort_sites2_sf, next_buffer_sf) %>% unlist() %>% as.numeric() %>% sum()
  sort_sites2_sf[m, "NumCloseSites"] <- num_close_sites
}

sort_sites2 <- sort_sites2_sf %>%
  mutate(Long = sf::st_coordinates(.)[,1], #add coordinates as fields
         Lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  arrange(-NumCloseSites) #orders sort_sites2 by number of close sites (within the radius)

sites <- rbind(sort_sites1,sort_sites2)
assert_that(site_count == nrow(sites), 
            msg = paste("Error: Sites were lost when arranging non-reproductive order for aggregation")) #stops the script if sites were due to non-standard radius in datasheet


while (nrow(sites) > 0){
  i=1
  
  #Select the next UniqueID to find sites to aggregate
  lookup_id <- sites$UniqueID_method[i]
  lookup_site <- sites %>% subset(UniqueID_method == lookup_id)
  #assert_that(lookup_site$UniqueID_method != "M7WEN83D_D_T") #Code for stopping the script at a particular site
  if (nrow(lookup_site) > 1) 
    lookup_site <- lookup_site[1,]
  
  #Find nearby sites  
  lookup_site_sf <- sf::st_as_sf(lookup_site, coords = c("Long","Lat"), crs = st_crs(4326)) #create sf object of the lookup_site
  lookup_buffer_sf <- sf::st_buffer(lookup_site_sf, units::set_units(lookup_site$Radius, degree)) #Create an sf object of the buffer of the lookup_site by its radius
  min_buffer_sf <- sf::st_buffer(lookup_site_sf, units::set_units(1, degree)) #Create an sf object with a 1 degree buffer of the lookup_site (used to in the long-route identification to keep close sites together)
  sites_sf <- sf::st_as_sf(sites, coords = c("Long","Lat"), crs = st_crs(4326)) #create an sf object of all sites... needs to be inside the loop so the sites available for intersection is updated each time
  
  nearby_sites <- sf::st_intersection(sites_sf, lookup_buffer_sf) %>% #identify sites within the buffer distance (radius) of the lookup-site
    select(!ends_with(".1")) %>% #remove duplicate fields
    mutate(Long = sf::st_coordinates(.)[,1], #add coordinates as fields
           Lat = sf::st_coordinates(.)[,2])
  
  ##If the site has a reproductive activity, only retain nearby_sites qith equal or smaller radius
  ##I.e., don't lose the resolution of the reproductive sites
  priority_activities <- c("NES","BRE","NUR","SPA","CAL","FLE","HAT") #List important behaviors which should not be aggregated with larger sites
  
  sitetype <- lookup_site$Activities %>%
    str_replace_all(" ","") %>%
    str_split(";") %>%
    unlist()
  
  if(any(sitetype %in% priority_activities)){
    nearby_sites <- nearby_sites %>% filter(Radius <= lookup_site$Radius)
  } 
  
  #RUN IF YOU WANT TO LOOK AT HOW A PARTICULAR LOOKUP SITE IS BEING AGGREGATED
  theme_set(theme_bw())
  mapcoords <- coord_fixed(xlim = c(150, 180), ylim = c(-50, -10))
  country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                                 data = map_data('world'),
                                 fill = "#FFFFFF", color = "#DDDDDD", #CECECE
                                 size = 0.15)
  ggplot() + country_shapes +
    geom_sf() + xlim(lookup_site$Long - 20,lookup_site$Long + 20) + ylim(lookup_site$Lat - 20, lookup_site$Lat + 20) +
    geom_sf(data = lookup_buffer_sf, alpha = .3) +
    geom_sf(data = min_buffer_sf, colour = "yellow") +
    geom_sf(data = sites_sf, colour = "blue") +
    geom_sf(data = lookup_site_sf, colour = "black") +
    geom_sf(data = nearby_sites)
  
  
  ##Filter sites that aren't in the same basin (e.g., opposite sides of Panama) 
  if(lookup_site$Basin[1] == "ATL"){
    nearby_sites <- nearby_sites %>% filter(Basin != "PAC")
  } else if(lookup_site$Basin[1] == "PAC") {
    nearby_sites <- nearby_sites %>% filter(Basin != "ATL")      
  }
  
  
  #If there are routes who have both endpoints in the nearby_sites dataframe, these routes will be lost when they are aggregated.
  #This loop checks the length of those routes and removes the site from nearby_sites if we do not want to lose routes over a certain length.
  #It also looks for nearby_sites which are close to any sites removed because they are part of long routes, and
  #removes them from the aggregation as well if they are close to a removed longroute site
  id_list <- as.character(nearby_sites$UniqueID) #get UniqueIDs of the nearby sites
  test_routes <- routes %>% subset(UniqueID_from %in% id_list & UniqueID_to %in% id_list) #Find routes for which both their FROM and TO sites are in the nearby_sites
  
  if (nrow(test_routes) > 0) {
    
    #Generate a new column with distances between To and From points in routes associate with nearby_sites
    test_routes_from_sf <- sf::st_as_sf(test_routes, coords = c("FromLong","FromLat"), crs = st_crs(4326)) #create sf object of the from points in test_routes 
    test_routes_to_sf <- sf::st_as_sf(test_routes, coords = c("ToLong","ToLat"), crs = st_crs(4326)) #create sf object of the to points in test_routes
    test_routes <- test_routes %>% 
      mutate(Dist = st_distance(test_routes_from_sf, test_routes_to_sf, by_element = TRUE, which = "Great Circle"))  #Calculates the distance between to and from points in test_routes
    
    
    #SUPERIMPORTANT!!!! This next line of code determines the cutoff for what length route we are willing to lose in the aggregation process.
    #This is also the only strict km distance measure/threshold left in the script, rather than degree buffers
    route_threshold <- units::set_units(500000, meter) 
    long_routes <- test_routes %>% subset(Dist > route_threshold) #set to aggregate pairs of sites that have a route <500km long
    
    if (nrow(long_routes) > 0) {
      long_routes_all_ids <- unique(c(long_routes$UniqueID_from, long_routes$UniqueID_to)) %>% .[. != lookup_site$UniqueID] #makes a vector of all ids in the remaining routes, and then removes the original lookup_site from that list
      long_route_sites_sf <- nearby_sites %>% subset(UniqueID_method %in% long_routes_all_ids)#pull out sites for UniqueIDs that are part of long routes
      long_route_sites_to_keep <- sf::st_intersection(long_route_sites_sf, min_buffer_sf) %>% #identify sites within the buffer distance (radius) of the lookup_site
        select(!ends_with(".1")) %>% #remove duplicate fields
        mutate(Long = sf::st_coordinates(.)[,1], #add coordinates as fields
               Lat = sf::st_coordinates(.)[,2])
      long_route_ids_to_keep <- st_drop_geometry(long_route_sites_to_keep) %>%
        select(UniqueID_method) %>%
        unlist() #spit out the UniqueIDs by themselves as a vector      
      
      long_route_sites_removed_ids <- long_routes_all_ids[!long_routes_all_ids %in% long_route_ids_to_keep]
      
      
      ##Identify sites that are very close to sites being removed to preserve long routes, and remove them as well so they can be aggregated with the closer site.
      for (q in long_route_sites_removed_ids) {
        removed_site_sf <- nearby_sites %>% subset(UniqueID_method == q) #iteratively grab a row from nearby_sites that matches the UniqueID of each value in long_route_sites_removed_ids
        #removed_site_sf <- removed_site_sf[1,] #gets rid of sites with duplicate UniqueIDs but different methods
        removed_site_min_buffer <- sf::st_buffer(removed_site_sf, units::set_units(removed_site_sf$Radius, degree)) #Create an sf object with degree radius buffer of the removed_site
        
        sites_near_removed_sites <- sf::st_intersection(nearby_sites, removed_site_min_buffer) %>% #identify sites within the buffer distance (radius) of the removed_site
          subset(!(UniqueID_method %in% c(lookup_site$UniqueID_method, removed_site_sf$UniqueID_method))) %>%#which are not either the original lookup_site or the removed site itself
          select(!ends_with(".1")) %>% #remove duplicate fields
          mutate(Long = sf::st_coordinates(.)[,1], #add coordinates as fields
                 Lat = sf::st_coordinates(.)[,2])
        
        #RUN IF YOU WANT TO LOOK AT HOW A PARTICULAR REMOVED SITE IS BEING AGGREGATED
        # theme_set(theme_bw())
        # country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
        #                                data = map_data('world'),
        #                                fill = "#FFFFFF", color = "#DDDDDD", #CECECE
        #                                size = 0.15)
        # ggplot() + country_shapes +
        #   geom_sf() + xlim(removed_site_sf$Long - 20,removed_site_sf$Long + 20) + ylim(removed_site_sf$Lat - 20, removed_site_sf$Lat + 20) +
        #   geom_sf(data = lookup_buffer_sf, alpha = .5, colour = "red") +
        #   geom_sf(data = removed_site_min_buffer, colour = "green", alpha = .5) +
        #   geom_sf(data = nearby_sites, colour = "blue") +
        #   geom_sf(data = removed_site_sf, colour = "green") +
        #   geom_sf(data = lookup_site_sf, colour = "red") +
        #   geom_sf(data = sites_near_removed_sites, colour = "yellow")
        
        
        sites_near_removed_sites <- sites_near_removed_sites %>% 
          mutate(DistToRemovedSite = as.numeric(st_distance(removed_site_sf, sites_near_removed_sites, by_element=TRUE, which = "Great Circle")), #Calculates the distance between all sites and the removed_site
                 DistToLookupSite = as.numeric(st_distance(lookup_site_sf, sites_near_removed_sites, by_element=TRUE, which = "Great Circle")))  #Calculates the distance between all sites and the lookup_site
        
        removed_neighbour_sites <- sites_near_removed_sites %>% subset(DistToRemovedSite < DistToLookupSite)
        
        sites_near_removed_sites_ids <- st_drop_geometry(removed_neighbour_sites) %>% 
          select(UniqueID) %>% unlist() #spit out the UniqueIDs of sites_near_removed_sites
        
        sites[sites$UniqueID %in% sites_near_removed_sites_ids,"Longroute"] <- ifelse(is.na(sites[sites$UniqueID %in% sites_near_removed_sites_ids,"Longroute"]),
                                                                                      paste("This site was not aggregated with", lookup_id, "because it was close to ", q," which was removed from that aggregation to preserve a long route"), 
                                                                                      paste(sites[sites$UniqueID %in% sites_near_removed_sites_ids,"Longroute"],
                                                                                            paste("This site was not aggregated with", lookup_id, "because it was close to ", q," which was removed from that aggregation to preserve a long route"), sep="; "))
      }
      
      nearby_sites <- nearby_sites %>% subset(!(UniqueID %in% c(long_route_sites_removed_ids, sites_near_removed_sites_ids))) #remove selected sites from nearby_sites, and allow them to be aggregated separately  
      
      sites[sites$UniqueID %in% long_route_sites_removed_ids,"Longroute"] <- ifelse(is.na(sites[sites$UniqueID %in% long_route_sites_removed_ids,"Longroute"]),
                                                                                    paste("This site was not aggregated with", lookup_id, "because it would result in the loss of a route longer than ", route_threshold, "m"), 
                                                                                    paste(sites[sites$UniqueID %in% long_route_sites_removed_ids,"Longroute"],
                                                                                          paste("This site was not aggregated with", lookup_id, "because it would result in the loss of a route longer than ", route_threshold, "m"), sep="; "))
      
      sites[sites$UniqueID %in% long_route_ids_to_keep,"Longroute"] <- ifelse(is.na(sites[sites$UniqueID %in% long_route_ids_to_keep,"Longroute"]),
                                                                              paste("This site was part of a long route, but was aggregated with", lookup_id, "because it is within 1 degree of the lookup_site"), 
                                                                              paste(sites[sites$UniqueID %in% long_route_ids_to_keep,"Longroute"],
                                                                                    paste("This site was part of a long route, but was aggregated with", lookup_id, "because it is within 1 degree of the lookup_site"), sep="; "))
    }    
  }
  
  
  #Code to find the centroid and aggregate sites if there are more than one nearby_sites
  if (nrow(nearby_sites) > 1) {
    nearby_centroid <- data.frame(geomean(cbind(nearby_sites$Long,nearby_sites$Lat))) %>% #Calculate the centroid of the new metasite based on remaining nearby_sites
      sf::st_as_sf(coords = c("x","y"), crs = st_crs(4326)) %>%
      mutate(Long = sf::st_coordinates(.)[,1],
             Lat = sf::st_coordinates(.)[,2])
    
    new_num_nearby_sites = nrow(nearby_sites)
    old_num_nearby_sites = 0
    
    #Run a loop that stops when there is no change in the number of nearby_sites due to a change in the location of the metasite centroid
    while (new_num_nearby_sites != old_num_nearby_sites) {
      old_num_nearby_sites = new_num_nearby_sites
      centroid_buffer_sf <- sf::st_buffer(nearby_centroid, units::set_units(lookup_site$Radius, degree))
      
      nearby_sites <- sf::st_intersection(nearby_sites, centroid_buffer_sf) %>%
        select(!ends_with(".1")) %>%  #gets rid of duplicate columns
        mutate(Long = sf::st_coordinates(.)[,1], #adds coordinates as columns
               Lat = sf::st_coordinates(.)[,2]) 
      
      if (nrow(nearby_sites) > 1){
        nearby_centroid <- data.frame(geomean(cbind(nearby_sites$Long,nearby_sites$Lat))) %>% #Calculate the centroid of the new metasite based on remaining nearby_sites
          sf::st_as_sf(coords = c("x","y"), crs = st_crs(4326)) %>%
          mutate(Long = sf::st_coordinates(.)[,1],
                 Lat = sf::st_coordinates(.)[,2])
      } else {
        nearby_centroid <- cbind(nearby_sites$Long[1],nearby_sites$Lat[1])
      }
      new_num_nearby_sites = nrow(nearby_sites)
    }
  }
  
  nearby_sites$clust <- last_clust_num
  nearby_sites <- nearby_sites %>% mutate(MetasiteNum = clust, MetaSiteID = site_ids[MetasiteNum])
  nearby_sites$MetaSiteID <- as.character(nearby_sites$MetaSiteID)
  
  #Add sites to output data.frame
  metasites <- rbind(metasites,nearby_sites)
  
  #Remove sites from "sites" dataframe that were added to the "metasites" dataframe
  sites <- sites %>% subset(!(UniqueID_method %in% nearby_sites$UniqueID_method))
  assert_that((nrow(metasites) + nrow(sites)) == site_count, msg = print("The script lost or gained sites in the aggregation process. Script stopped."))
  
  last_clust_num <- last_clust_num +1
  print(paste(nrow(sites),"left to aggregate.")) #added so the person running the script can follow along and know if it gets hung
}

assert_that(nrow(metasites) == site_count, msg = print("The number of final sites in the metasites file is not the same as the original number of sites in the sitesfile"))
metasites <- st_drop_geometry(metasites) %>% mutate(Comments = Longroute)

##CODE TO IDENTIFY DUPLICATE DATA
#Find sites that literally have the same info in key columns. 
#Generally happens when we have Telemetry and Mark-Recapture at the same site, but they have been broken out in to two rows
samedata_all_sites <- metasites %>% data.table() %>% arrange(desc(NumIDs)) 
samedata_samesites <- samedata_all_sites %>% mutate(Duplicate = duplicated(samedata_all_sites, by = c("ZoteroID","Long","Lat"))) #,"Method"
samedata_samesites_ids <- samedata_samesites %>% filter(Duplicate == TRUE) %>% select(UniqueID_method) %>%
  as.matrix() %>% as.vector() #just spitting out a clean list of ids

metasites[metasites$UniqueID_method %in% samedata_samesites_ids,"Duplicate"] <- TRUE
metasites[metasites$UniqueID_method %in% samedata_samesites_ids,"Rationale"] <- "Same site"
metasites[metasites$UniqueID_method %in% samedata_samesites_ids,"Comments"] <- ifelse(is.na(metasites[metasites$UniqueID_method %in% samedata_samesites_ids,"Comments"]),
                                                                                      paste("This site is considered a duplicate because it had the same coordinates and ZoteroID as another site."), 
                                                                                      paste(metasites[metasites$UniqueID_method %in% samedata_samesites_ids,"Comments"],
                                                                                            paste("This site is considered a duplicate because it had the same coordinates and ZoteroID as another site."), sep=";"))

#Identify sites that have a route between them, but have been aggregated in the same metasite
routes <- routes %>% left_join(subset(metasites, select = c(UniqueID_from, MetaSiteID)), by = "UniqueID_from", keep=FALSE) %>%
  rename(MetaSiteID_from = MetaSiteID)
routes <- routes %>% left_join(subset(metasites, select = c(UniqueID_to, MetaSiteID)), by = "UniqueID_to", keep=FALSE) %>%
  rename(MetaSiteID_to = MetaSiteID)
routes <- routes %>% mutate(Self = (MetaSiteID_from == MetaSiteID_to))
self_route <- subset(routes, Self == TRUE)
duplicate_list = c()
if(nrow(self_route) > 0) {
  for (x in 1:nrow(self_route)) {
    dup_test_ids <- c(self_route$UniqueID_from[x], self_route$UniqueID_to[x])
    dup_test_sites <- metasites %>% subset(UniqueID %in% dup_test_ids)
    if (is.na(max(dup_test_sites$NumIDs))) { #If all NUMIDs are NAs, then just choose the first one
      kept_duplicate_site_id <- dup_test_sites[1,"UniqueID_method"] %>% as.character()
    } else {
      kept_duplicate_site_id <- dup_test_sites[(which.max(dup_test_sites$NumIDs)),"UniqueID_method"][1]  %>% as.character()
    }
    duplicate_site_ids <- dup_test_sites %>% subset(!UniqueID_method == kept_duplicate_site_id, select = "UniqueID_method") %>%
      as.matrix() %>% as.vector() #just spitting out a clean list of ids
    for (g in duplicate_site_ids) {
      if (!g %in% duplicate_list){
        metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Duplicate"] <- TRUE
        metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Rationale"] <- ifelse(is.na(metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Rationale"]),
                                                                                           "Self-route", 
                                                                                           paste(metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Rationale"],"Self-route", sep=";"))
        metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Comments"] <- ifelse(is.na(metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Comments"]),
                                                                                          paste("This site was aggregated with a site it had a route to:", kept_duplicate_site_id), 
                                                                                          paste(metasites[metasites$UniqueID_method %in% duplicate_site_ids,"Comments"],
                                                                                                paste("This site was aggregated with a site it had a route to:", kept_duplicate_site_id), sep=";"))
      }
      duplicate_list <- c(duplicate_list, duplicate_site_ids)
    }
  }
}

#Identify sites that might be duplicates, but didn't share a route and are in different papers
samedata_all_sites <- data.frame(samedata_all_sites)
#for (z in 1:nrow(samedata_all_sites)) { #used for testing
duplicate_list = c()
metasites$Duplicate <- as.character(metasites$Duplicate) 
while (nrow(samedata_all_sites) > 0) {
  samedata_test_site <- metasites %>% filter(UniqueID_method == samedata_all_sites[1,"UniqueID_method"]) #Make sure the final bit of this line reads: "samedata_all_sites[1,"UniqueID_method"]". Sometimes I change the number during testing
  samedata_test_meta_id <- samedata_test_site$MetaSiteID
  samedata_test_id_method <- samedata_test_site$UniqueID_method
  samedata_test_sampling_method <- samedata_test_site$Method
  
  year_list <- unlist(
    str_split(
      gsub(pattern = "\\s",
           replacement = "",
           as.character(samedata_test_site$Years)),
      ";")
  )
  
  lastnames <- function(author) {
    lastname <- str_split(author,",")[[1]][1]
    return(lastname)
  }    
  
  author_list <- unlist(
    str_split(
      gsub(pattern = "\\s",
           replacement = "",
           as.character(samedata_test_site$Author)),
      ";")
  )
  author_names <- as.character(lapply(author_list, lastnames))
  
  
  
  #Breakout the non-duplicated sites to compare against the original test site for other signs of duplicated data  
  samedata_agg_sites <- metasites %>% subset(MetaSiteID == samedata_test_meta_id & #get all sites that have been aggregated in the same metasite
                                               !(UniqueID_method == samedata_test_id_method) & #remove one with the exact same UniqueID including method
                                               (Method == samedata_test_sampling_method) & #keep only sites with the same sampling method 
                                               !(ZoteroID == samedata_test_site$ZoteroID)) #remove those sites from the same paper, as they would only be duplicates if they shared a route (see above)
  
  #ifelse(nrow(samedata_agg_sites)>0, print(paste("Found one:",samedata_test_site$UniqueID)), print(paste(samedata_test_site$UniqueID, "has no duplicates"))) #used for testing
  #} #used for testing
  for (z in 1:nrow(samedata_agg_sites)){
    #z=11 #Used for testing
    agg_year_list <- unlist(
      str_split(
        gsub(pattern = "\\s",
             replacement = "",
             as.character(samedata_agg_sites$Years[z])),
        ";")
    )
    
    agg_author_list <- unlist(
      str_split(
        gsub(pattern = "\\s",
             replacement = "",
             as.character(samedata_agg_sites$Author[z])),
        ";")
    )
    agg_author_names <- as.character(lapply(agg_author_list, lastnames))
    
    #same_years <- sum(year_list %in% agg_year_list) > 0
    same_years <- ifelse(all(year_list %in% agg_year_list) | all(agg_year_list %in% year_list), TRUE, FALSE)
    same_authors <- sum(author_names %in% agg_author_names, na.rm = TRUE) > 0
    which_authors <- author_names[which(author_names %in% agg_author_names)]
    
    if (same_years == TRUE & same_authors == TRUE) {
      samedata_duplicate_sites <- rbind(samedata_test_site,samedata_agg_sites[z,]) %>% replace_na(list(NumIDs = 0))
      if (samedata_duplicate_sites[1,"NumIDs"] == samedata_duplicate_sites[2,"NumIDs"]) {
        samedata_kept_site_id <- samedata_duplicate_sites[1,"UniqueID_method"] %>% as.character()
        samedata_dropped_site_id <- samedata_duplicate_sites[2,"UniqueID_method"] %>% as.character()    
      } else {
        samedata_kept_site_id <- samedata_duplicate_sites[(which.max(samedata_duplicate_sites$NumIDs)),"UniqueID_method"] %>% as.character()
        samedata_dropped_site_id <- samedata_duplicate_sites[(which.min(samedata_duplicate_sites$NumIDs)),"UniqueID_method"] %>% as.character()
      }
      #print(paste("Duplicate site:", samedata_dropped_site_id, sep="")) #used for testing
      #print(c("Duplicate_list:", duplicate_list)) #used for testing
      if (!samedata_dropped_site_id %in% duplicate_list){
        #print(paste("Test site:", samedata_test_id_method, "; Duplicate site:", samedata_dropped_site_id, sep="")) #used for testing
        metasites[metasites$UniqueID_method == samedata_dropped_site_id,"Duplicate"] <- TRUE
        metasites[metasites$UniqueID_method == samedata_dropped_site_id,"Rationale"] <- ifelse(is.na(metasites[metasites$UniqueID_method %in% samedata_dropped_site_id,"Rationale"]),
                                                                                               "Same data", 
                                                                                               paste(metasites[metasites$UniqueID_method %in% samedata_dropped_site_id,"Rationale"],"Same data", sep=";"))
        samedata_comment <- paste("This site was aggregated with", samedata_kept_site_id, "that had the same year(s) and one or more of the same authors")
        metasites[metasites$UniqueID_method == samedata_dropped_site_id,"Comments"] <- ifelse(is.na(metasites[metasites$UniqueID_method %in% samedata_dropped_site_id,"Comments"]),
                                                                                              samedata_comment, 
                                                                                              paste(metasites[metasites$UniqueID_method %in% samedata_dropped_site_id,"Comments"],
                                                                                                    samedata_comment, sep=";"))
        samedata_all_sites <- samedata_all_sites %>% filter(!UniqueID_method %in% samedata_duplicate_sites$UniqueID_method)
      }
      duplicate_list = unique(c(duplicate_list, samedata_dropped_site_id))
    }
  }
  samedata_all_sites <- samedata_all_sites %>% filter(!UniqueID_method %in% c(samedata_test_id_method,duplicate_list))
}


metasites <- metasites %>% mutate(across(c("Duplicate"), ~replace_na(., "FALSE")))
#metasites <- metasites %>% mutate(Comments = paste(Comments,Longroute))


metasites$MetaSiteID <- as.character(metasites$MetaSiteID)
metasites$MetaID <- paste("META", spp, sep = "")
metasites <- metasites %>% relocate(c(MetaID, MetaSiteID, Duplicate, Rationale, Comments), .before = Location) %>% 
  select(!c("clust","MetasiteNum","Longroute","NumCloseSites"))

## exportfile <- paste("AggregatedSites\\Testing\\georeferenced_sites_", spp, "_unsupervised_v13.4.csv", sep="")
exportfile <- paste("C:/LEGOLAS/R/Lit. Rev/Aggregation/", spp, "_unsupervised_v13.4.JRpo.csv", sep="")
write.csv(metasites, exportfile, quote=TRUE, row.names = FALSE)

