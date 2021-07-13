#' Serve only for downloading raw data. Script create 'source_data/env' 
#' directory. Script downloads Geomorpho90, GADM, and SRTM data. Data 
#' from Copernicus Land Monitoring Service have to be obtained manually as
#' is described. Script makes virtual rasters ('.vrt') from downloaded
#' data if applicable.
library(gdalUtils)
library(raster)
env_path <- file.path("..","..","source_data","env")
dir.create(env_path,showWarnings=F)

if_download.file <- function (url,destfile,method){
    if (!file.exists(destfile)) {
    download.file(url ,destfile,method=method) }
}

#' Tree cover density & Grasslands
dir.create(file.path(env_path,"tree_cover_dens"),showWarnings=F)
dir.create(file.path(env_path,"gra"),showWarnings=F)

#' Get the data manually from Copernicus Land Monitoring Service (10m resolution): 
#https://land.copernicus.eu/pan-european/high-resolution-layers/forests/tree-cover-density/status-maps/tree-cover-density-2018
#https://land.copernicus.eu/pan-european/high-resolution-layers/grassland/status-maps/grassland-2018

# Extract the `DATA` directory into 'source_data/env/tree_cover_dens' in project directory, to get:

# source_data/env/tree_cover_dens/DATA/*.tif
# source_data/env/gra/DATA/*.tif

tcd_l <- list.files(file.path(env_path,"tree_cover_dens","DATA"),pattern=".tif$",full.names=T)
gdalbuildvrt(tcd_l,file.path(env_path,"tcd.vrt"),overwrite=T)

gra_src <- list.files(file.path(env_path,"gra","DATA"),full.names=T,pattern=".tif")
gdalbuildvrt(gra_src,file.path(env_path,"gra.vrt"),overwrite=T)

#' # Geomorpho90

geom_vars <- c("slope","cti","geom","tri","elev-stdev")

tiles <- c("n50e015","n45e015","n50e010","n45e010")

for (gv in geom_vars){
    dir.create(file.path(env_path,gv),showWarnings=F)
    if_download.file(paste0("https://cloud.sdsc.edu/v1/AUTH_opentopography/hosted_data/OTDS.012020.4326.1/raster/",gv,"/",gv,"_90M_n30e000.tar.gz"),file.path(env_path,gv,paste0(gv,"_90M_n30e000.tar.gz")),"wget")
    l <- normalizePath(list.files(file.path(env_path,gv),full.names=T))
    l <- l[1]
    tif_l <- file.path(l,untar(l,list=T))
    tif_tiles <- paste0(gv,"_90M_",tiles,".tif")
    untar(l,files=tif_tiles,exdir=file.path(env_path,gv))
    gdalbuildvrt(file.path(env_path,gv,tif_tiles),file.path(env_path,paste0(gv,".vrt")))
}

#' # Admin boundaries
admin_path <- file.path("..","..","source_data","admin")
dir.create(admin_path,showWarnings=F)
getData("GADM",country="CZE",level=0,path=admin_path)

#' # SRTM (used only for map background in QGIS)
srtm_path <- file.path(env_path,"srtm")
dir.create(srtm_path,showWarnings=F)

getData('SRTM', lon=13, lat=49,path=srtm_path)
getData('SRTM', lon=16, lat=49,path=srtm_path)
getData('SRTM', lon=13, lat=51,path=srtm_path)
getData('SRTM', lon=16, lat=51,path=srtm_path)

gdalbuildvrt(list.files(srtm_path,pattern="*.tif",full.names=T),file.path(env_path,paste0("srtm",".vrt")))
