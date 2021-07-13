#'The script reads the downloaded and provided data and combines them into
#' one file 'ocss.gpkg'. Than recent records from the areal margin, and 
#' reccurrence records are marked in the 'occ_type' column. The documented
#' human introduced records are removed.

library(raster)
library(rgdal)

# define CRS
crs_5514 <- "+proj=krovak +lat_0=49.5 +lon_0=24.83333333333333 +alpha=30.28813972222222 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +towgs84=589,76,480,0,0,0,0 +units=m +no_defs" 

crs_4326 <- "+proj=longlat +datum=WGS84 +no_defs"

# odra_survey
odra_survey <-  read.csv(file.path("..","source_data","pk","odra_survey.csv"))

odra_survey <- SpatialPointsDataFrame(data = odra_survey["year"],
                                coords = cbind(odra_survey$lon, odra_survey$lat),
                                proj4string = CRS(crs_4326)
                                )

odra_survey$source <- "odra_mapping"
odra_survey$unc <- NA
odra_survey$occ_type <- "occs"
odra_survey$occ_type[40] <- "recent_B"


# national mapping
nat_mapping <- readOGR(file.path("..","source_data","pk","orthotpera_map.gpkg"))

nat_mapping$source <- "nat_mapping"
nat_mapping$unc <- NA
nat_mapping$occ_type <- "occs"

#NDOP
ndop_occs <- read.csv(file.path("..","source_data","ndop","ruspolia_ndop_tab.csv"),dec=",")
ndop_occs$occ_type <- "occs"
ndop_occs[ndop_occs$ID_NALEZ==49888668,"occ_type"] <- "recent_A"

ndop_occs <- SpatialPointsDataFrame(
                        coords=ndop_occs[c("X","Y")],
                        data=data.frame(unc=ndop_occs$CXPRESNOST,year=ndop_occs$DATUM_OD,occ_type=ndop_occs$occ_type),
                        proj4string=CRS(crs_5514)
                    )
                    
ndop_occs$source = "ndop"
ndop_occs$year = as.numeric(substring(ndop_occs$year,1,4))

#remove human introduction occurences on the west
ndop_occs <- ndop_occs[ndop_occs@coords[,1]>-700000,]

ndop_occs <- spTransform(ndop_occs,CRS(crs_4326))
# Holusa2007

hol2007 <- read.csv(file.path("..","source_data","pk","Holusa_et_al_2007.csv"))
hol2007 <- SpatialPointsDataFrame(hol2007[,2:1],
                                    hol2007["year"],
                                proj4string=CRS(crs_4326))

hol2007$source <- "hol2007"
hol2007$unc <- NA
hol2007$occ_type <- "occs"
hol2007$occ_type[1] <- "initial_A"
hol2007$occ_type[3] <- "initial_B"

occs <- rbind(odra_survey,ndop_occs,nat_mapping,hol2007)

writeOGR(occs,file.path("..","processed_data","ocss.gpkg"),layer="occs",driver="GPKG",overwrite=T)
