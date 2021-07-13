#' The script create 'processed_data/env' dircetory, reads all environmental
#' data used in HSM, and reproject and resample them to 15 arc sec 
#' (ca.500 at the equator ) resolution. Resulted data are stored in 
#' 'processed_data/env' as geotiffs ('.tif')

library(raster)
library(rgdal)

proc_path <- file.path("..","processed_data","env")
dir.create(proc_path,showWarnings=F)

cz_bounds <- readRDS(file.path("..","source_data","admin","gadm36_CZE_0_sp.rds"))
# create small buffer for clipping raster layers to ensure include occurence points near boarders
cz_bounds_b <- raster::buffer(cz_bounds,0.1)

#' Geom
geom_src <- raster(file.path("..","source_data","env","geom.vrt"))

geom <- crop(geom_src,cz_bounds_b)
geom_1 <- aggregate(geom, 5, fun=function(x, na.rm=T) {mean(x==1, na.rm=na.rm)})

writeRaster(geom_1,file.path("..","processed_data","env","geom_1.tif"),overwrite=T)

geom_1 <- raster(file.path("..","processed_data","env","geom_1.tif"))

#' Geom - slope, tci
geom_st_src <- stack(list.files(file.path("..","source_data","env"),pattern="tri.vrt|cti.vrt|slope.vrt|elev-stdev.vrt",full.names=T))
geom_st <- crop(geom_st_src,cz_bounds_b)
geom_st <- aggregate(geom_st, 5)

writeRaster(geom_st,file.path("..","processed_data","env",paste0(names(geom_st),".tif")),bylayer=T,overwrite=T)

#' Tree Cover Density
tcd_src <- raster(file.path("..","source_data","env","tcd.vrt"))
tcd <- crop(tcd_src,spTransform(cz_bounds_b,crs(tcd_src)))
tcd_rp <- aggregate(tcd, 50,mean)
tcd_rp <- projectRaster(tcd_rp, geom_1, method = "bilinear")
tcd_rp <- resample(tcd_rp, geom_1, method = "bilinear")

writeRaster(tcd_rp,file.path("..","processed_data","env","tcd.tif"),overwrite=T)

#' Grassland
gra <- raster(file.path("..","source_data","env","gra.vrt"))
gra_1 <- aggregate(gra, 50, fun=function(x, na.rm=T) {mean(x==1, na.rm=na.rm)})
gra_1 <- projectRaster(gra_1, geom_1, method = "bilinear")
gra_1 <- resample(gra_1, geom_1, method = "bilinear")

writeRaster(gra_1,file.path("..","processed_data","env","gra.tif"),overwrite=T)
