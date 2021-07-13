# This script calculates Geographic, LCP and passage LCP distances.
library(gdistance)
library(raster)
library(rgdal)
library(rgeos)

crs_5514 <- "+proj=krovak +lat_0=49.5 +lon_0=24.83333333333333 +alpha=30.28813972222222 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +towgs84=589,76,480,0,0,0,0 +units=m +no_defs" 

#' read and transform occurence data
occs <- readOGR(file.path("..","processed_data","ocss.gpkg"),layer="occs")
occs$id <- as.numeric(row.names(occs))
occs <- occs[occs$year>2005,]
occs <- spTransform(occs,crs("+init=epsg:4326"))

#' read HSM output and make transition object for gDsitance functions
base_r <- raster(file.path("..","processed_data","hsm.tif"))
values(base_r)[is.na(values(base_r))] <- 0.0001

trCost <- transition(base_r, mean, directions=16)
trCost <- geoCorrection(trCost, type="c")

#' function for calculating LCP, returns lines with distances
lcp_dist <- function(origin,dest,trans = trCost){
    lcp <- shortestPath(trCost, origin, dest, output="SpatialLines")
    dists <- gLength(spTransform(lcp,crs(crs_5514)),byid=T)/1000
    lcp <- SpatialLinesDataFrame(sl=lcp,data=data.frame(dist=dists), match.ID = F)
    lcp@data <- cbind(lcp@data,dest@data)
    return(lcp)
}

#' function for shifting occurrence points to centroids of the HSM pixels
#' , this ensure to claculate geographic distances between points to be 
#' comparable with LCP based distances, which are calcualted from centroids 
rasterize_points <- function (p){
    r <- rasterize(coordinates(p),base_r)
    r_p <- rasterToPoints(r,spatial=T)
    new_coords <- matrix(ncol=2)
    for(i in 1:nrow(p)) {#inspired by https://gis.stackexchange.com/a/121732/36001
        d <- spDistsN1(r_p, p[i,])
          new_coords <- rbind(new_coords,coordinates(r_p)[which.min(d),])
        }  
    shift_r_p <- SpatialPointsDataFrame(list(new_coords[-1,1],new_coords[-1,2]),data=p@data)
    crs(shift_r_p)<-   crs(p)
    return(shift_r_p)
}

#' function for calculating Geograhical distances, returns lines with distances
geo_dist <- function(origin,dest){
    r_origin <- rasterize_points(origin)
    r_dest <- rasterize_points(dest)
    geo <- lapply(split(r_dest, r_dest$id), function(x) Lines(list(Line(rbind(coordinates(x),coordinates(r_origin)))), x$id[1]))
    geo <- SpatialLinesDataFrame(sl=SpatialLines(geo),data=data.frame(dist=rep(NA,length(geo))), match.ID = F)
    crs(geo) <- crs(r_origin)
    geo$dist <- gLength(spTransform(geo,crs("+init=epsg:5514")),byid=T)/1000
    geo@data <- cbind(geo@data,r_dest@data)
    return(geo)
}
#' function for calculating passage LCP, returns lines with distances, 
#' and passage rasters
pass_dist <- function (origin,dest,b_r = base_r){
    for (i in 1:length(dest)){
        st <- Sys.time()
        print(length(dest))
        print(i)
        cr_base_r <- crop(b_r,raster::buffer(rbind(origin, dest[i,]),10000))
        plot(cr_base_r)
        cr_trCost <- transition(cr_base_r, mean, directions=16)
        cr_trCost <- geoCorrection(cr_trCost, type="c")
        pass_r <- passage(cr_trCost, origin, dest[i,],0.0001)
        trCost_pass <- transition(pass_r, mean, directions=16)
        trCost_pass <- geoCorrection(trCost_pass, type="c")
        pass_lcp <- shortestPath(trCost_pass, origin, dest[i,], output="SpatialLines")
        pass_lcp_dist <- gLength(spTransform(pass_lcp,crs("+init=epsg:5514")),byid=T)/1000
        pass_lcp <- SpatialLinesDataFrame(sl=pass_lcp,data=data.frame(dist=pass_lcp_dist), match.ID = F)
        pass_r <- resample(pass_r,b_r)
        if (i==1){
            pass_stack <- stack(pass_r)
            pass_lcps <- pass_lcp
        } else {
            pass_stack <- stack(pass_stack,pass_r)
            pass_lcps <- rbind(pass_lcps,pass_lcp)
        }
        names(pass_stack)[i] <- paste0("pass_",dest$id[i])
        gc()
        print(Sys.time()-st)
    }
    pass_lcps@data <- cbind(pass_lcps@data,dest@data)
    return(c(pass_lcps=pass_lcps,pass_stack=pass_stack))
}

#' function for removing points near by origin point when calculating 
#' distances. This is necessary due to is impossible calculate distances
#' between points if there are within one raster cell
rm_close <- function(origin,occs){
    x <- extract(raster::buffer(origin,300),occs)
    rem_pts <- x[!is.na(x$poly.ID == 1),"point.ID"]
    occs <- occs[-rem_pts,]
    return(occs)
}


#' # 1 Odra -regional scale

#' function for stepwise annual calculation of distances in Odra basin 
#' with selected method ('fd' is argument for function described above).
#' Returns lines with distances, and origin points for each year. The
#' procedure is described in the paper: "We defined the origin site
#' for each year as the point that had the lowest average distance to all sites discovered
#' in the next year. In addition, in each step we removed the occurrences that were more
#' distant than 3 km in the opposite direction (west) than the spread direction. We repeated
#' this for every year."
recurs_dist <- function(fd,years=2017:2020){
    for (i in years){
        if (i==2017){
            origins <- occs_odra[occs_odra$year==i-1,]
        } else {
            origins <- rm_close(occs_odra[occs_odra$year==i,],occs_odra[occs_odra$year==i-1,])
        }
        dests <- occs_odra[occs_odra$year==i,]
        dests <- dests[dests@coords[,1] > max(origins@coords[,1])-0.027 ,]
        len_or <- length(origins)
        df_or_len <- data.frame(feat=rep(NA,len_or))
        for (f in 1:len_or){
            l <- fd(origins[f,],dests)
            if (typeof(l)=="list"){
                l <- l[[1]]
            }
            df_or_len$feat[f] <- f
            df_or_len$mean[f] <- mean(l$dist)
            df_or_len$max[f] <- max(l$dist)
        }
        o <- origins[which(df_or_len$mean==min(df_or_len$mean)),]
        l <- fd(o,dests)
        if (typeof(l)=="list"){
               l <- l[[1]]
        }
        if (i==min(years)){
            lines_o <- l
            origin <- o 
        } else {
            lines_o <- rbind(lines_o,l)
            origin <- rbind(origin,o)
        }
        gc()
    }
    return(c(lines=lines_o,origins=origin))
}

#' filter only data from Odra survey, perform calculations and writes results
occs_odra <- occs[occs$source=="odra_mapping",]
origin <- occs_odra[occs_odra$year==2016,]

lcps_o <- recurs_dist(lcp_dist)
geo_o <- recurs_dist(geo_dist)
pass_o <- recurs_dist(pass_dist)

names(lcps_o$lines)[1] <- "lcp"
names(geo_o$lines)[1] <- "geo"
names(pass_o$lines)[1] <- "pass"

df_o <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), list(lcps_o$lines, geo_o$lines[,c(1,6)], pass_o$lines[,c(1,6)]))

write.csv(df_o,file.path("..","results","df_o.csv"))
writeOGR(pass_o$lines,file.path("..","processed_data","pass_o.gpkg"),driver="GPKG",layer="lcps",overwrite=T)
writeOGR(geo_o$lines,file.path("..","processed_data","geo_o.gpkg"),driver="GPKG",layer="lcps",overwrite=T)
writeOGR(lcps_o$lines,file.path("..","processed_data","lcps_o.gpkg"),driver="GPKG",layer="lcps",overwrite=T)
writeOGR(pass_o$origins,file.path("..","processed_data","pass_o_origins.gpkg"),driver="GPKG",layer="origins",overwrite=T)
writeOGR(geo_o$origins,file.path("..","processed_data","geo_o_origins.gpkg"),driver="GPKG",layer="origins",overwrite=T)
writeOGR(lcps_o$origins,file.path("..","processed_data","lcps_o_origins.gpkg"),driver="GPKG",layer="origins",overwrite=T)

#' # 2 National - landscape scale

#' function for dividing passage raster stack. This creates raster stack
#' with paths only belongs to defined origin ('A' or 'B'), this is only 
#' written for use in 'nat_dist' function
div_pass_stack <- function (pass_stack,div_lines,ori_name){
    print("start")
    list_ids <- div_lines[div_lines$origin==ori_name,]$id  
    pass_names <- paste0("pass_",list_ids)
    print(pass_names)
    pass_stack <- pass_stack[[pass_names]]
    print(pass_stack)
    values(pass_stack)[is.na(values(pass_stack))] <- 0
    rDiv <- max(max(pass_stack) * (1 - min(pass_stack)) - min(pass_stack), 0)
    plot(rDiv)
    plot(div_lines[div_lines$origin==ori_name,],add=T)
    writeRaster(rDiv,file.path("..","processed_data",paste0("rDiv_",ori_name,".tif")),overwrite=T)
}
#' function for calculating paths and distances on landscape scale. 
#' Paths are further divided between initial points 'A' or 'B', depend 
#' on from which initial points was shorter distance to other certain point.
nat_dist <- function (fd,fname){
    for (p in c("A","B")){
        ori <- fd(occs_init[grepl(p,occs_init$occ_type),],occs_dest)
        if (typeof(ori)=="list"){
            ori_stack <- ori[[2]]
            saveRDS(ori_stack,file.path("..","processed_data",paste0("stack_",p,".RDS")))
            ori <- ori[[1]]
            if (p=="A"){
                A_stack <- ori_stack
            } else {
                B_stack <- ori_stack
            }
        }
        writeOGR(ori,file.path("..","processed_data",paste0(fname,"_",p,".gpkg")),driver="GPKG",layer="lcps",overwrite=T)
        ori$origin <- p
        if (p=="A"){
            A <- ori
        } else {
            B <- ori
        }
    }
    div_A <- A[which(A$dist < B$dist),]
    div_B <- B[which(A$dist > B$dist),]
    div <- rbind(div_A,div_B)
    writeOGR(div,file.path("..","processed_data",paste0(fname,"_div.gpkg")),driver="GPKG",layer="lcps",overwrite=T)
    if (exists("ori_stack")){
        div_pass_stack(A_stack,div,"A")
        div_pass_stack(B_stack,div,"B")
    }
    return(div)
}

occs_init <- occs[occs$occ_type %in% c("initial_A","initial_B"),]

occs_dest <- rm_close(occs_init,occs)#remove points close to initial
#' filter out Odra basin data except one recent areal margin occurence
occs_rec_B <- occs_dest[occs_dest$occ_type=="recent_B",]
occs_dest <- occs_dest[occs_dest$source!="odra_mapping",]
occs_dest <- rbind(occs_rec_B,occs_dest)


lcps_nat <- nat_dist(lcp_dist,"lcp")
geo_nat <- nat_dist(geo_dist,"geo")
pass_nat <- nat_dist(pass_dist,"pass")

names(lcps_nat)[1] <- "lcp"
names(geo_nat)[1] <- "geo"
names(pass_nat)[1] <- "pass"

df_nat <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), list(lcps_nat, geo_nat[,c(1,6)], pass_nat[,c(1,6)]))

write.csv(df_nat,file.path("..","results","df_nat.csv"))
