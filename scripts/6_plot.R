#' This script used ggplot2 to generate all the necessary graphic output 
#' directly for the manuscript or to generate graphs for maps in QGIS 
#' (use '7_qgis_maps.py').


library(ggplot2)
library(ggforce)
library(reshape)

plot_path <- file.path("..","plots")
dir.create(plot_path,showWarnings=F)

y_cond_max <- function (df,dist_var){
    df <- aggregate(x=df[dist_var], by = list(df$year),FUN= max)
    names(df) <- c("year",dist_var)
    y <- min(df$year):max(df$year)
    if (length(y)!=length(unique(df$year))){
        ndf <- data.frame(year = y[!y %in% df$year],dist=0)
        names(ndf)[2] <- dist_var
        df <- rbind(df,ndf)
    }
    df <- df[order(df$year),]
    df[dist_var] <- cummax(df[dist_var])
    return(df)
}

dist_ests <- readRDS(file.path("..","results","estimates.RDS"))

lm_eqn <- function(m){#inspired by https://stackoverflow.com/a/7549819/3984070
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)[italic(Adj)]^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 4),
             r2 = format(summary(m)$adj.r.squared , digits = 3)))
    as.character(as.expression(eq))
}

geo_col <- "#FA817D"
lcp_col <- "#FCFE83"
pass_col <- "#80FC7C"
stroke_col <- "#585858"
exp_dpi <- 600
l <- c("Geographic distance (km)","LCP distance (km)","Passage LCP distance (km)")

# Odra
lcps_odr <- read.csv(file.path("..","results","df_o.csv"))
df_odr <- Reduce(function(x,y) merge(x = x, y = y, by = "year"), list(aggregate(lcp ~ year,lcps_odr, max), aggregate(geo ~ year,lcps_odr, max), aggregate(pass ~ year,lcps_odr, max)))
df_odr <- cumsum(df_odr)
df_odr$year <- 2017:2020

geo_lm_plot <- ggplot(lcps_odr,aes(as.factor(year),geo,ymin=10,ymax=20))+
    geom_violin(fill=geo_col,col=stroke_col,alpha=0.25,lwd=0.2,adjust=0.7,draw_quantiles = 0.5)+
    geom_sina(fill=geo_col,color=stroke_col,shape = 21, size =2 , stroke = 0.2)+
    ylab(l[1]) + 
    xlab("Year") +
    theme_light()

ggsave(
  file.path("..","plots","geo_lm_plot.svg"),
  plot = geo_lm_plot,
  width = 93,
  height = 70,
  units="mm"
)


lcp_lm_plot <- ggplot(lcps_odr,aes(as.factor(year),lcp,ymin=10,ymax=20))+
    geom_violin(fill=lcp_col,col=stroke_col,alpha=0.25, lwd = 0.2,adjust=0.7,draw_quantiles = 0.5)+
    geom_sina(fill=lcp_col,color=stroke_col,shape = 21, size =2 , stroke = 0.2)+
    ylab(l[2]) + 
    xlab("Year") +
    theme_light()

ggsave(
  file.path("..","plots","lcp_lm_plot.svg"),
  plot = lcp_lm_plot,
  width = 93,
  height = 70,
  units="mm"
)

pass_lm_plot <- ggplot(lcps_odr,aes(as.factor(year),pass,ymin=10,ymax=20))+
    geom_violin(fill=pass_col,col=stroke_col,alpha=0.25, lwd = 0.2,adjust=0.7,draw_quantiles = 0.5)+
    geom_sina(fill=pass_col,color=stroke_col,shape = 21, size =2 , stroke = 0.2)+
    ylab(l[3]) + 
    xlab("Year") +
    theme_light()

ggsave(
  file.path("..","plots","pass_lm_plot.svg"),
  plot = pass_lm_plot,
  width = 93,
  height = 70,
  units="mm"
)

names(df_odr)[c(3,2,4)] <- gsub(" distance \\(km)","",l)
df_o_melt <- melt(df_odr,"year")
names(df_o_melt)[3] <- "dist"
names(df_o_melt)[2] <- "method"
df_o_melt$method <- factor(df_o_melt$method,levels(df_o_melt$method)[c(2,1,3)])
max_dist_plot_o <- ggplot(df_o_melt,aes(year,dist,colour = method))+
    scale_color_manual(values=c(geo_col,lcp_col,pass_col))+
    stat_smooth(method = "lm",fill=NA,size=1)+
    geom_point(size =2)+
    geom_point(color=stroke_col,shape = 21, size =2 , stroke = 0.2) +
    geom_label(x = 2017.75, y = 55 , label = lm_eqn(dist_ests$lm_o_geo),fill=geo_col,col=stroke_col, parse = TRUE)+
    geom_label(x = 2017.75, y = 50 , label = lm_eqn(dist_ests$lm_o_lcp),fill=lcp_col,col=stroke_col, parse = TRUE)+
    geom_label(x = 2017.75, y = 45, label = lm_eqn(dist_ests$lm_o_pass),fill=pass_col,col=stroke_col, parse = TRUE)+
    ylab("Distance (km)")+ 
    xlab("Year") +
    labs(color='Method') +
    theme_light()

ggsave(
  file.path("..","plots","max_dist_plot_o.svg"),
  plot = max_dist_plot_o,
  width = 24,
  height = 12,
  units="cm"
)

ggsave(
  file.path("..","plots","max_dist_plot_o.png"),
  plot = max_dist_plot_o,
  width = 24,
  height = 12,
  units = "cm",
  dpi = exp_dpi
)


# Landscape

lcps_nat <- read.csv(file.path("..","results","df_nat.csv"))
df_nat <- Reduce(function(x,y) merge(x = x, y = y, by = "year"), list(y_cond_max(lcps_nat,"lcp"), y_cond_max(lcps_nat,"geo"), y_cond_max(lcps_nat,"pass")))
lcps_nat <- lcps_nat[lcps_nat$year>2006,]
lcps_nat[nrow(lcps_nat)+1,"year"] <- 2008

geo_lm_plot_nat <- ggplot(lcps_nat,aes(as.factor(year),geo,ymin=10,ymax=60)) +
    geom_violin(fill=geo_col,col=stroke_col,alpha=0.25,lwd=0.2,scale="width",draw_quantiles = 0.5)+
    geom_sina(fill=geo_col,color=stroke_col,shape = 21, size =1 , stroke = 0.2,scale="width",maxwidth=0.5)+
    ylab(l[1]) + 
    xlab("Year") +
    scale_x_discrete(labels=c("","","","2010","","","","","2015","","","","","2020")) +
    theme_light()

ggsave(
  file.path("..","plots","geo_lm_plot_nat.svg"),
  plot = geo_lm_plot_nat,
  width = 93,
  height = 70,
  units="mm"
)

lcp_lm_plot_nat <- ggplot(lcps_nat,aes(as.factor(year),lcp,ymin=10,ymax=60)) +
    geom_violin(fill=lcp_col,col=stroke_col,alpha=0.25,lwd=0.2,scale="width",draw_quantiles = 0.5)+
    geom_sina(fill=lcp_col,color=stroke_col,shape = 21, size =1 , stroke = 0.2,scale="width",maxwidth=0.5)+
    ylab(l[2]) + 
    xlab("Year") +
    scale_x_discrete(labels=c("","","","2010","","","","","2015","","","","","2020")) +
    theme_light()
    
ggsave(
  file.path("..","plots","lcp_lm_plot_nat.svg"),
  plot = lcp_lm_plot_nat,
  width = 93,
  height = 70,
  units="mm"
)

pass_lm_plot_nat <- ggplot(lcps_nat,aes(as.factor(year),pass,ymin=10,ymax=60)) +
    geom_violin(fill=pass_col,col=stroke_col,alpha=0.25,lwd=0.2,scale="width",draw_quantiles = 0.5)+
    geom_sina(fill=pass_col,color=stroke_col,shape = 21, size =1 , stroke = 0.2,scale="width",maxwidth=0.5)+
    ylab(l[3]) + 
    xlab("Year") +
    scale_x_discrete(labels=c("","","","2010","","","","","2015","","","","","2020")) +
    theme_light()

ggsave(
  file.path("..","plots","pass_lm_plot_nat.svg"),
  plot = pass_lm_plot_nat,
  width = 93,
  height = 70,
  units="mm"
)

names(df_nat)[c(3,2,4)] <- gsub(" distance \\(km)","",l)
df_nat_melt <- melt(df_nat,"year")
names(df_nat_melt)[3] <- "dist"
names(df_nat_melt)[2] <- "method"
df_nat_melt$method <- factor(df_nat_melt$method,levels(df_nat_melt$method)[c(2,1,3)])
max_dist_plot_lsc <- ggplot(df_nat_melt,aes(year,dist,colour = method))+
    scale_color_manual(values=c(geo_col,lcp_col,pass_col))+
    stat_smooth(method = "lm",fill=NA,size=1)+
    geom_point(size =2)+
    geom_point(color=stroke_col,shape = 21, size =2 , stroke = 0.2)+
    ylim(0,200) +
    geom_label(x = 2008.75, y = 190 , label = lm_eqn(dist_ests$lm_nat_geo),fill=geo_col,col=stroke_col, parse = TRUE)+
    geom_label(x = 2008.75, y = 175 , label = lm_eqn(dist_ests$lm_nat_lcp),fill=lcp_col,col=stroke_col, parse = TRUE)+
    geom_label(x = 2008.75, y = 160, label = lm_eqn(dist_ests$lm_nat_pass),fill=pass_col,col=stroke_col, parse = TRUE)+ 
    ylab("Distance (km)") + 
    xlab("Year") +
    labs(color='Method') +
    theme_light()
    
ggsave(
  file.path("..","plots","max_dist_plot_lsc.svg"),
  plot = max_dist_plot_lsc,
  width = 24,
  height = 16,
  units="cm"
)

ggsave(
  file.path("..","plots","max_dist_plot_lsc.png"),
  plot = max_dist_plot_lsc,
  width = 24,
  height = 16,
  units = "cm",
  dpi = exp_dpi
)

# Landscape - Weibull's model
w_df <- read.csv(file.path("..","results","w_model.csv"),row.names=1)
lms_sum <- read.csv(file.path("..","results","lms_subs.csv"),row.names=1)
c("Geographic","LCP","Passage LCP")

w_df[1] <- log10(2:1001)
names(w_df)[c(3,2,4)] <- gsub(" distance \\(km)","",l)
names(w_df)[1] <- "l_scale"
w_df <- melt(w_df,"l_scale")
w_df$variable <- factor(w_df$variable ,levels(w_df$variable)[c(2,1,3)])

names(lms_sum)[c(5,3,7)] <- gsub(" distance \\(km)","",l)
lms_sum <- melt(lms_sum[c(2:3,5,7)],"sample_pr")
lms_sum <- rbind(lms_sum,list(100,"Passage LCP",dist_ests$lm_nat_pass$coefficients[2] ),list(100,"Geographic",dist_ests$lm_nat_geo$coefficients[2] ),list(100,"LCP",dist_ests$lm_nat_lcp$coefficients[2] ))
lms_sum$variable <- factor(lms_sum$variable ,levels(lms_sum$variable)[c(2,1,3)])

w_sampling <- ggplot(w_df,aes(l_scale,value,colour=variable))+
    geom_line(lwd=1)+
    geom_point(data=lms_sum[lms_sum$sample_pr!=100,],mapping=aes(log10(sample_pr+1),value,fill=variable),color=stroke_col,shape = 21, size = 2.5 , stroke = 0.5)+
    geom_point(data=lms_sum[lms_sum$sample_pr==100,],mapping=aes(log10(sample_pr+1),value,fill=variable),color=stroke_col,shape = 22, size = 6 , stroke = 0.5)+
    ylim(5,17)+
    scale_fill_manual(values=c(geo_col,lcp_col,pass_col)) +
    scale_color_manual(values=c(geo_col,lcp_col,pass_col)) +
    scale_x_continuous(breaks=log10(c(11,31,101,1001)), labels=c(10,30,100,1000),limit=c(1,3))+
    ylab("Predicted expansion rate (km/y)") + 
    xlab("Sampling effort (%)")+
    labs(fill='Method',color='Method') +
    theme_light()

ggsave(
  file.path("..","plots","w_sampling_model.png"),
  plot = w_sampling,
  width = 24,
  height = 16,
  units = "cm",
  dpi = exp_dpi
)
