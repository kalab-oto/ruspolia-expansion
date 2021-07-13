library(drc)

#' function for calculating conditional maximum - cumulative maximum 
#' distance - only the maximum distance from every year was used, if 
#' maximum is lwoer than# the maximum from previous years, it is replaced 
#' with the previous maximum
y_cond_max <- function (df,dist_var,st_y=2007){
    df <- aggregate(x=df[dist_var], by = list(df$year),FUN= max)
    names(df) <- c("year",dist_var)
    y <- st_y:2020
    if (length(y)!=length(unique(df$year))){
        ndf <- data.frame(year = y[!y %in% df$year],dist=0)
        names(ndf)[2] <- dist_var
        df <- rbind(df,ndf)
    }
    df <- df[order(df$year),]
    df[dist_var] <- cummax(df[dist_var])
    return(df)
}
#' this function only wraps 'lm' function printing a 'plot'
lm_plot <- function (df){
    lm_df <- lm(df[,2]~df$year)
    print(summary(lm_df))
    print(confint(lm_df))
    newx <- df$year
    plot(df$year,df[,2], pch = 16, xlab = "year", ylab = "dist")
    abline(lm_df)
    conf_interval <- predict(lm_df, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
    matlines(newx,conf_interval[,2:3])
    return(lm_df)
}

#' function for explorative analysis of distances in particular yers 
#' in Odra basin. Returns annual annual mean and sd
signle_dists <- function (lm_o){
    odr_dists <- lm_o$model[[1]]
    for (i in 1:length(odr_dists)){
        odr_dists[i] <- sum(odr_dists[i])- sum(odr_dists[1:i-1])
    }
    df <- c(dist=odr_dists,mean=mean(odr_dists),sd=sd(odr_dists))
    return(round(df,2))
    }
#' function for explorative analysis of two paths from origin points to 
#' areal margin (landscape)
most_dist_paths <- function (lcps,dist_var){
    max_dist_A <- max(lcps[lcps$origin=="A",dist_var])
    max_dist_B <- max(lcps[lcps$origin=="B",dist_var])
    df <- data.frame(path=c("A","B"),distance=round(c(max_dist_A,max_dist_B),1),km_y=round(c(max_dist_A,max_dist_B)/14,1))
    return(df)
}

#' # Odra - regional scale
lcps_odr <- read.csv(file.path("..","results","df_o.csv"))
df_odr <- Reduce(function(x,y) merge(x = x, y = y, by = "year"), list(aggregate(lcp ~ year,lcps_odr, max), aggregate(geo ~ year,lcps_odr, max), aggregate(pass ~ year,lcps_odr, max)))
df_odr <- cumsum(df_odr)
df_odr$year <- 2017:2020

lm_o_geo <- lm_plot(df_odr[c(1,3)])
signle_dists(lm_o_geo)

lm_o_lcp <- lm_plot(df_odr[c(1,2)])
signle_dists(lm_o_lcp)

lm_o_pass <- lm_plot(df_odr[c(1,4)])
signle_dists(lm_o_pass)

#' # National - landscape scale
lcps_nat <- read.csv(file.path("..","results","df_nat.csv"))
lcps_nat <- lcps_nat[lcps_nat$year>2006,]

lcps_nat_rec_B <- lcps_nat[lcps_nat$occ_type=="recent_B",]
lcps_nat <- lcps_nat[lcps_nat$source!= "odra_mapping",]
lcps_nat_rec_B <- rbind(lcps_nat,lcps_nat_rec_B)

df_nat <- Reduce(function(x,y) merge(x = x, y = y, by = "year"), list(y_cond_max(lcps_nat,"lcp"), y_cond_max(lcps_nat,"geo"), y_cond_max(lcps_nat,"pass")))

lm_nat_geo <- lm_plot(df_nat[c(1,3)])
most_dist_paths(lcps_nat_rec_B,"geo")

lm_nat_lcp <- lm_plot(df_nat[c(1,2)])
most_dist_paths(lcps_nat_rec_B,"lcp")

lm_nat_pass <- lm_plot(df_nat[c(1,4)])
most_dist_paths(lcps_nat_rec_B,"pass")

#' ## Weibull's model
nrep=10000
pr=seq(10,90,10)
rep_lms <- data.frame(id_seed=rep(NA,nrep))
lms_sum <- data.frame(repl=nrep,sample_pr=pr,mean_lcp=NA,sd_lcp=NA,mean_geo=NA,sd_geo=NA,mean_pass=NA,sd_pass=NA)

#' replication with subsequently removing part of the data.
for (p in pr){
    for (i in 1:nrep){
        set.seed(i)
        rep_lms$id_seed[i] <- i
        df_sub <- df_nat[sample(nrow(df_nat),round(nrow(df_nat)*p/100)),]
        df_sub_cond_max <- Reduce(function(x,y) merge(x = x, y = y, by = "year"), list(y_cond_max(df_sub,"lcp"), y_cond_max(df_sub,"geo"), y_cond_max(df_sub,"pass")))
        rep_lms$geo[i] <- lm(df_sub_cond_max$geo~df_sub_cond_max$year)$coefficients[2] 
        rep_lms$lcp[i] <- lm(df_sub_cond_max$lcp~df_sub_cond_max$year)$coefficients[2] 
        rep_lms$pass[i] <- lm(df_sub_cond_max$pass~df_sub_cond_max$year)$coefficients[2] 
    }
    lms_sum[lms_sum$sample_pr==p,]$mean_lcp <- mean(rep_lms$lcp)
    lms_sum[lms_sum$sample_pr==p,]$sd_lcp <- sd(rep_lms$lcp)
    lms_sum[lms_sum$sample_pr==p,]$mean_geo <- mean(rep_lms$geo)
    lms_sum[lms_sum$sample_pr==p,]$sd_geo <- sd(rep_lms$geo)
    lms_sum[lms_sum$sample_pr==p,]$mean_pass <- mean(rep_lms$pass)
    lms_sum[lms_sum$sample_pr==p,]$sd_pass <- sd(rep_lms$pass)
}

write.csv(lms_sum,file.path("..","results","lms_subs.csv"))

#' perform Weibull's model
new_pr=1:1000
pr_df <- data.frame(new_pr=new_pr,mean_lcp = NA, mean_geo = NA, mean_pass = NA)

for (i in names(lms_sum[c(3,5,7)])){
    model <- drm(lms_sum[,i]~lms_sum$sample_pr,fct=LL.3())
    model$call$formula <-as.formula(paste0("lms_sum$",i,"~lms_sum$sample_pr"))
    model_select <- mselect(model, list(LL.2(), LL.3(), LL.3u(), LL.4(),LL.5(),W1.2(),W1.3(),W1.4(),W2.2(),W2.3(),W2.4(),BC.4(),BC.5(),LL2.2(),LL2.3(),LL2.3u(),LL2.4(),LL2.5(),AR.2(),AR.3(),MM.2(),MM.3()), linreg=TRUE, icfct=AIC)
    fname <- names(model_select[,"Res var"][1])
    model <- drm(lms_sum[,i]~lms_sum$sample_pr,fct=do.call(fname,list()))
    plot(model)
    summary(model)
    print(model)
    pr_df[,i] <-predict(model,data.frame(sample_pr=new_pr))
}

pr_df$mean_geo[10]/lm_nat_geo$coefficients[2][[1]]
pr_df$mean_lcp[10]/lm_nat_lcp$coefficients[2][[1]]
pr_df$mean_pass[10]/lm_nat_pass$coefficients[2][[1]]

lm_nat_geo$coefficients[2][[1]]/max(pr_df$mean_geo)
lm_nat_lcp$coefficients[2][[1]]/max(pr_df$mean_lcp)
lm_nat_pass$coefficients[2][[1]]/max(pr_df$mean_pass)

write.csv(pr_df,file.path("..","results","w_model.csv"))

####################
dist_ests <- mget(c("lm_o_lcp","lm_o_geo","lm_o_pass","lm_nat_lcp","lm_nat_geo","lm_nat_pass"))

saveRDS(dist_ests,file.path("..","results","estimates.RDS"))
