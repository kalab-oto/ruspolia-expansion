#' The script executes the entire HSM, including filtering the finding 
#' data (i.e. spatial thinning), selecting environmental variables with 
#' VIFm, background data generation, choosing features and beta multiplier 
#' for maxent (ENMeval), performing and evaluating maxent model. Final 
#' model output is exported as geotiff (.tif), and ENMeval and model 
#' evaluation results are exported to .csv.

library(raster)
library(rgdal)
library(gdalUtils)
library(spThin)
library(sdm)
library(usdm)
library(dismo)
library(ecospat)
library(ENMeval)
require(devtools)
    source_url("https://raw.githubusercontent.com/RWunderlich/SEDI/master/R/sedi.R")

#' read boundaries data
cz_bounds <- getData("GADM",country='CZE',level=0,path=file.path("..","source_data","admin"))

#' # OCC data

occs <- readOGR(file.path("..","processed_data","ocss.gpkg"),layer="occs")

plot(occs,pch=20,col="red")
plot(cz_bounds, add = TRUE)

occs$species <- "ruspolia"

rownames(occs@data) <- 1:nrow(occs)
occs$id <- rownames(occs@data)

#' perform spatial thinning

thinned_sp <- spThin::thin(loc.data =  as.data.frame(occs),
                lat.col = "coords.x1",
                long.col = "coords.x2",
                spec.col = "species",
                thin.par = 5,#km
                reps = 1,
                out.dir = "../processed_data",
                locs.thinned.list.return = TRUE,
                write.files = FALSE
            )

thinned_sp <- as.data.frame(thinned_sp)
thinned_sp$id <- rownames(thinned_sp)
thinned_sp <- merge(thinned_sp,occs,by="id")

thinned_sp <- SpatialPointsDataFrame(
                    coords=thinned_sp[c("coords.x1","coords.x2")],
                    data=data.frame(thinned_sp[c("source","species")]),
                    proj4string=CRS("EPSG:4326")
                )
plot(occs,pch=20,col="red")
plot(thinned_sp,pch=20,col="green",add=T)
plot(cz_bounds, add = TRUE)

#' # ENV data
env_source <- list.files(file.path("..","processed_data","env"),full.names=T,pattern=".tif$")[-c(6)]
env_data <- mask(stack(env_source),raster::buffer(cz_bounds,0.1))

# excluding variables based on VIF
na_mask <- sum(is.na(env_data))>0
na_mask[na_mask==1] <- NA
env_data <- mask(env_data,na_mask)

vif <- vifcor(env_data,.7)
env_data <- exclude(env_data,vif)

#' # Model

colnames(thinned_sp@coords) <- c("x","y")
thinned_sp$sp_name <- 1
sp_train <- thinned_sp[,"sp_name"]

# generate background data
bg_sp <- randomPoints(env_data, 10000)
bg <- extract(env_data,bg_sp)
bg <- cbind(bg_sp,bg)
bg <- as.data.frame(bg)


# Model calibration
rms <- seq(0.5,10,0.5)
fcs <- c("LQ", "LQP", "LQT", "LQH", "LQHT", "LQTP", "LQHP", "LQHPT")

enm_eval_results <- data.frame()
for (i in fcs){
    enm_eval <- ENMevaluate(occ=sp_train@coords, env = env_data, bg.coords=bg_sp,method='randomkfold', kfolds=5, fc=i, RMvalues=rms, algorithm='maxent.jar',parallel=F)
    enm_eval_results <- rbind(enm_eval@results,enm_eval_results)
    gc()
    print(unique(enm_eval_results$features))
    }
enm_eval_results <- enm_eval_results[which(!is.na(enm_eval_results$AICc)),]
enm_eval_results$delta.AICc <- enm_eval_results$AICc - min(enm_eval_results$AICc)
comb <- as.character(enm_eval_results$features[enm_eval_results$delta.AICc==0])
b_beta <- enm_eval_results$rm[enm_eval_results$delta.AICc==0] 

b_args <- c()
if (!grepl("L", comb)){
    b_args <- c(b_args,"nolinear")
}
if (!grepl("H", comb)){
    b_args <- c(b_args,"nohinge")
}
if (!grepl("Q", comb)){
    b_args <- c(b_args,"noquadratic")
}
if (!grepl("P", comb)){
    b_args <- c(b_args,"noproduct")
}
if (grepl("T", comb)){
    b_args <- c(b_args,"threshold")
}

#' Model

d <- sdmData(train=sp_train, predictors=env_data,bg=bg)

m <- sdm(sp_name~.,data=d,
        methods=c('maxent'),
        replication='cv',
        cv.folds=5,
        modelSettings=list(maxent=list(beta=b_beta,args=c(b_args,"noremoveDuplicates", "noautofeature"))),
        n=10
        )

roc(m,smooth=T)
rcurve(m)

#' Evaluation
mean_ev <- function (mdl = m,testing_set="test.dep",measures=measur){
    df <- as.data.frame(apply(getEvaluation(mdl,w=(getModelId(mdl)),stat=measures,opt=2,wtest= testing_set)[-1], 2, mean))
    df <- cbind(df,as.data.frame(apply(getEvaluation(mdl,w=(getModelId(mdl)),stat=measures,opt=2,wtest= testing_set)[-1], 2, sd)))
    df <- cbind(df,as.data.frame(apply(getEvaluation(mdl,w=(getModelId(mdl)),stat=measures,opt=2,wtest= testing_set)[-1], 2, max)))
    df <- cbind(df,as.data.frame(apply(getEvaluation(mdl,w=(getModelId(mdl)),stat=measures,opt=2,wtest= testing_set)[-1], 2, min)))
    
    colnames(df) <- c("mean","sd","max","min")
    if ("boyce" %in% measures){
        print("calculating boyce")
        boyce_li <- c()
        for (repl in getModelId(mdl)){
            cv_bg <- c(mdl@models$sp_name$maxent[[repl]]@evaluation$train@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation$train@observed==0],mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed==0])
            cv_pres <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed==1]
            b <- ecospat.boyce(cv_bg, cv_pres,PEplot=F)
            boyce_li <- c(boyce_li,b$Spearman.cor)
        }
        df <- rbind(df,c(mean(boyce_li),sd(boyce_li),max(boyce_li),min(boyce_li)))
        rownames(df)[nrow(df)] <- "boyce" 
    }
    if ("SEDI" %in% measures){
        print("calculating SEDI")
        sedi_li <- c()

        for (repl in getModelId(mdl)){
            th_val <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@threshold_based[2,2]
            o <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed
            p <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted
            cmx <- sdm:::.cmx(o>th_val,p>th_val) 
            sedi_li <- c(sedi_li,sedi(cmx[1,1],cmx[1,2],cmx[2,1],cmx[2,2])[[2]])
        }
        df <- rbind(df,c(mean(sedi_li),sd(sedi_li),max(sedi_li),min(sedi_li)))
        rownames(df)[nrow(df)] <- "sedi_li" 
    }

    return(df)
}

measur <- c('boyce', 'SEDI')

#' Variable imprtance
var_imp <- getVarImp(m,id=getModelId(m),method="maxent",wtest="test.dep")@varImportanceMean$corTest[2]*100
var_imp <- round(var_imp,0)

var_imp

eval_dep <- mean_ev(m,"test.dep")
eval_dep

#' Predict
env_pr_data <- exclude(stack(env_source),vif)
env_pr_data <- crop(env_pr_data,raster::buffer(occs,10000))

pr_eu <- ensemble(m,env_pr_data, setting=list(method='weighted',stat='AUC',opt=2))

plot(pr_eu)
plot(cz_bounds, add = TRUE)
plot(thinned_sp,pch=20,col="red",add=T)

#' Export results
writeRaster(pr_eu,file.path("..","processed_data","hsm.tif"),overwrite=T)
write.csv(enm_eval_results, file.path("..","processed_data","ENMeval.csv")) 
write.csv(eval_dep,file.path("..","processed_data","eval_dep.csv"))


