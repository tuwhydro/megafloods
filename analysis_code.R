##############################################################################
#This script is meant to reproduce the main results of: 
#Megafloods in Europe can be anticipated from observations in hydrologically similar catchments
#Bertola et al. (2023)
##############################################################################
# Reading in the data

stations <- read.csv("gauges.csv",header=T)
flood_data <- read.csv("floods.csv",header=T)


##############################################################################
# Fig. 1 & EDT 1 & EDT 2


# Megafloods selection
find.mf.func <- function(stations.tab, amax.tab, start.year){
  
  library(robustbase)
  library(extRemes)

  megafloods.df <- data.frame(id=character(), LAT=double(), LON=double(), area=double(),
                              Q_MF=double(), year_MF=integer(), date_MF=lubridate::ymd())
  
  nn<-0 # count number of outliers
  for(i in 1: dim(stations.tab)[1]){
    
    dummy <- amax.tab[amax.tab$id==stations.tab$id[i],]
    dummy <- dummy[!is.na(dummy$q),]
    if(length(dummy$q)<10) next #exclude short series
    
    # megafloods criteria:

    # 1. outlier 
    dummy.quant <- quantile(dummy$q, probs = c(0.25,0.75))
    dummy.upper <- dummy.quant[2]+3*(dummy.quant[2]-dummy.quant[1])
    index.outlier <- which(dummy$q>=dummy.upper)
    
    if(length(index.outlier)==0) next
    index.MF <- 0 
    for(tt in index.outlier){
      ty <- dummy$year[tt]
      nyears.before <- length(dummy$q[dummy$year<ty])
      
      # 3. after start.year and with at least 20 years of data before event
      if(nyears.before<20) next
      if(ty<start.year) next
      
      # 2. surprising (Tl/Tsl>3) record-breaking event 
      qmax.old <- dummy$q[dummy$year<=ty]
      local.gev <- fevd(qmax.old, type = "GEV", method = "Lmoments")
      neP <- pextRemes(local.gev,q=qmax.old)
      RP <- 1/(1-neP)
      RP[RP==Inf] <- 10^5
      rRP <- RP[length(RP)]/max(RP[-length(RP)])
      if(rRP>3) index.MF <- c(index.MF,tt)
    }
    index.outlier <- index.MF[-1]
    
    if(length(index.outlier)>0){
      
      
      for(oo in index.outlier){
        nn <- nn+1
        megafloods.df[nn,1] <- stations.tab$id[i]
        megafloods.df[nn,2] <- stations.tab$lat[i]
        megafloods.df[nn,3] <- stations.tab$lon[i]
        megafloods.df[nn,4] <- as.numeric(stations.tab$area[i])
        megafloods.df[nn,5] <- dummy$q[oo] #outlier discharge
        megafloods.df[nn,6] <- dummy$year[oo] # outlier year
        megafloods.df[nn,7] <- as.Date(paste(dummy$year[oo],dummy$month[oo],dummy$day[oo], sep="-"),format="%Y-%m-%d") # outlier date

      }
    }
  }
  return(megafloods.df)
}
megafloods.df <- find.mf.func(stations, flood_data, 2000)

# Mean and Max discharge
stations$mean_q <- rep(NA,dim(stations)[1])
stations$max_q <- rep(NA,dim(stations)[1])
for(i in 1:dim(stations)[1]){
  qq <- flood_data[flood_data$id == stations$id[i],2]
  stations$mean_q[i] <- round(mean(qq),5)
  stations$max_q[i] <- max(qq)
  rm(qq)
}

# Envelope curve for 5 regions
env.func <- function(stations.tab, amax.tab){
  library(scales)
  library(quantreg)
  ENV_tab <- matrix(data=NA, ncol=6,nrow=5)
  colnames(ENV_tab) <- c("slope_qr_median","intercept_qr_median","slope_qr_env","intercept_qr_env","intercept_env", "n_stat_years")
  rownames(ENV_tab) <- 1:5
  ENV_tab <- as.data.frame(ENV_tab)
  for(rr in 1:5){
    print(rr)

    dummy.stations <- stations.tab[stations.tab$regions1 == rr,]
    dummy.q <- amax.tab[amax.tab$id %in% dummy.stations$id,]
    dummy.q <- dummy.q[dummy.q$q>0,]
    dummy.q$qlog10 <- log10(dummy.q$q)
    dummy.q$Alog10 <- rep(NA, length(dummy.q$id))
    for(i in 1:length(dummy.stations$id)){
      ind.stat <- which(dummy.q$id == dummy.stations$id[i])
      dummy.q$Alog10[ind.stat] <- log10(dummy.stations$area[i])
    }    
    # quantile regression z=0.5
    p<-1
    dummy.rq <- rq(dummy.q$qlog10 ~ dummy.q$Alog10, tau=0.5, data=dummy.q, na.action=na.omit,
                   method="fnc", R=cbind(0,rbind(diag(p),-diag(p))), r=c(rep(-999,p),-rep(0,p)))
    ENV_tab[rr,1] <- dummy.rq$coefficients[2] #slope of qr
    ENV_tab[rr,2] <- dummy.rq$coefficients[1] #intercept of qr
    # quantile regression z=0.999
    dummy.rq <- rq(dummy.q$qlog10 ~ dummy.q$Alog10, tau=0.999, data=dummy.q, na.action=na.omit,
                   method="fnc", R=cbind(0,rbind(diag(p),-diag(p))), r=c(rep(-999,p),-rep(0,p)))
    ENV_tab[rr,3] <- dummy.rq$coefficients[2] #slope of qr
    ENV_tab[rr,4] <- dummy.rq$coefficients[1] #intercept of qr
    intercept <- -unlist(ENV_tab[rr,3])*log10(dummy.stations$area)+log10(dummy.stations$max_q)
    ENV_tab[rr,5] <- max(intercept)   #intercept of envelope
    ENV_tab[rr,6] <- dim(dummy.q)[1] #nstatyears
  }
  return(ENV_tab)
  
}
ENV_reg <- env.func(stations, flood_data)

# EDT 1
ext.data.tab1.func <- function(envelope.tab){
  ExtDataTable1 <- as.data.frame(matrix(data=NA, ncol=5, nrow=5))
  colnames(ExtDataTable1) <- c("Slope_Envelope","Intercept_Envelope","Intercept_Envelope_1000km2","Slope_median_regression","Intercept_median_regression")
  rownames(ExtDataTable1) <- rownames(envelope.tab)
  ExtDataTable1$Slope_Envelope <- envelope.tab$slope_qr_env
  ExtDataTable1$Intercept_Envelope <- envelope.tab$intercept_env
  ExtDataTable1$Intercept_Envelope_1000km2 <- 10^(envelope.tab$intercept_env+envelope.tab$slope_qr_env*log10(1000))
  ExtDataTable1$Slope_median_regression <- envelope.tab$slope_qr_median
  ExtDataTable1$Intercept_median_regression <- envelope.tab$intercept_qr_median
  ExtDataTable1 <- round(ExtDataTable1,2)
  return(ExtDataTable1)
}
ExtDataTable1 <- ext.data.tab1.func(ENV_reg)

# EDT 2
ExtDataTable2 <- matrix(data=NA, ncol=5, nrow=5)
colnames(ExtDataTable2) <- c("Number_of_gauges","Avg_Ratio_max_mean_Q","Number_of_targets", "Percentage_of_targets","Avg_Ratio_max_mean_Q_targets")
rownames(ExtDataTable2) <- 1:5
ExtDataTable2 <- as.data.frame(ExtDataTable2)
for(rr in 1:5){
  dummy.stations <- stations[stations$regions1==rr,]
  ExtDataTable2$Number_of_gauges[rr] <- dim(dummy.stations)[1] 
  ExtDataTable2$Avg_Ratio_max_mean_Q[rr] <- round(mean(dummy.stations$max_q/dummy.stations$mean_q),2)
  dummy.mf <- dummy.stations[dummy.stations$id %in% megafloods.df$id,]
  ExtDataTable2$Number_of_targets[rr] <- dim(dummy.mf)[1]
  ExtDataTable2$Percentage_of_targets[rr] <-  round(ExtDataTable2$Number_of_targets[rr]*100/ExtDataTable2$Number_of_gauges[rr],1) 
  ExtDataTable2$Avg_Ratio_max_mean_Q_targets[rr] <- round(mean(dummy.mf$max_q/dummy.mf$mean_q),2) 
  rm(dummy.mf,dummy.stations)
}



##############################################################################
# Fig. 2 & Fig. 4

# mean and cv up to target years
mhq.cv.func <- function(stations.tab, amax.tab, target.years){
  library(robustbase)
  EU_MHQ_CV_list <- vector(mode = "list", length = length(target.years))
  names(EU_MHQ_CV_list) <- target.years
  for(YY in 1:length(target.years)){ #loop on target years
    print(target.years[YY])
    TARGET_YEAR <- target.years[YY]
    EU_MHQ_CV <- as.data.frame(matrix(data=NA, ncol=3, nrow = dim(stations.tab)[1]))
    colnames(EU_MHQ_CV) <- c("MHQ","CV","MHQ_100")
    EU_MHQ_CV <- cbind(stations.tab[1:4],EU_MHQ_CV)
    for(rrr in 1:dim(EU_MHQ_CV)[1]){ #loop on stations
      dummy.tab <- amax.tab[(amax.tab$id==EU_MHQ_CV[rrr,1])&(amax.tab$year<TARGET_YEAR),]
      dummy.tab <- dummy.tab[!is.na(dummy.tab$q),]
      if(dim(dummy.tab)[1]<10) next #at least 10 years before target
      dummy.Q <- dummy.tab$q
      dummy.Y <- dummy.tab$year
      #exclude outliers
      dummy.quant <- quantile(dummy.Q, probs = c(0.25,0.75))
      dummy.upper <- dummy.quant[2]+3*(dummy.quant[2]-dummy.quant[1])
      dummy.Y<-dummy.Y[dummy.Q<=dummy.upper]
      dummy.Q<-dummy.Q[dummy.Q<=dummy.upper]
      EU_MHQ_CV[rrr,5] <- mean(dummy.Q)
      EU_MHQ_CV[rrr,6] <- sd(dummy.Q)/EU_MHQ_CV[rrr,5]
      rm(dummy.tab,dummy.Q,dummy.Y,dummy.quant,dummy.upper)
    }
    #normalize MHQ to 100km2
    area.regr <- lm(log10(EU_MHQ_CV$MHQ)~log10(EU_MHQ_CV$area))
    EU_MHQ_CV[,7] <- 10^(log10(EU_MHQ_CV$MHQ)-area.regr$coefficients[2]*log10(EU_MHQ_CV$area)+area.regr$coefficients[2]*log10(100))
    EU_MHQ_CV_list[[YY]] <- EU_MHQ_CV
  }
  return(EU_MHQ_CV_list)
}
EU_MHQ_CV_list  <- mhq.cv.func(stations, flood_data, 2000:2022)


# predict mf
predict.mf.function <- function(mf.list, stations.tab, target.years, amax.tab, mhq.cv.tab, env.tab){
  
  #specify weights and distance limit for finding donor group
  w1<-1
  w2<-1
  w3<-1
  DIST <-1
  
  library(quantreg)
  
  #prepare lists for saving results
  REG_MEGAFLOODS_list <- vector(mode = "list", length = 5) #List megafloods per each region
  REG_ENVELOPE_param <- vector(mode = "list", length = 5) #Envelope parameters
  REG_DONORS_list <- vector(mode = "list", length = 5) #Donor list per region per megaflood

  for(REG in 1:5){ #loop on Regions
    REG_list.stations <- stations.tab$id[stations.tab$regions1==REG]
    MEGAFLOODS_list <- mf.list[mf.list$id %in% REG_list.stations,]
    DONORS_list <- vector(mode = "list", length = dim(MEGAFLOODS_list)[1]) 
    ENVELOPE_param <- as.data.frame(matrix(data=NA, ncol=6,nrow=dim(MEGAFLOODS_list)[1]))
    colnames(ENVELOPE_param) <- c("slope","intercept.region","intercept.mf","q.predicted","n_statyears_donors","n_stat_donors")
    
    for(i in 1:dim(MEGAFLOODS_list)[1]){ # loop on megafloods
      print(i)
      
      TARGET_YEAR <- MEGAFLOODS_list$year_MF[i]
      TARGET_Q <- MEGAFLOODS_list$Q_MF[i]
      TARGET_DATE <- MEGAFLOODS_list$date_MF[i]
      code.mf <- MEGAFLOODS_list$id[i]
      EU_MHQ_CV <- mhq.cv.tab[[which(target.years==TARGET_YEAR)]]
      EU_MHQ_CV <- EU_MHQ_CV[!is.na(EU_MHQ_CV$MHQ),]
      EU_MHQ_CV <- EU_MHQ_CV[EU_MHQ_CV$id %in% REG_list.stations,]
      
      #distances 
      ind.stat <- which(EU_MHQ_CV$id==code.mf)
      sd1 <- (log10(EU_MHQ_CV$area)-mean(log10(EU_MHQ_CV$area)))/sd(log10(EU_MHQ_CV$area))
      sd2 <- (log10(EU_MHQ_CV$MHQ_100)-mean(log10(EU_MHQ_CV$MHQ_100)))/sd(log10(EU_MHQ_CV$MHQ_100))
      sd3 <- (log10(EU_MHQ_CV$CV)-mean(log10(EU_MHQ_CV$CV)))/sd(log10(EU_MHQ_CV$CV))
      d1 <- sd1-sd1[ind.stat]
      d2 <- sd2-sd2[ind.stat]
      d3 <- sd3-sd3[ind.stat]
      dd <- ((d1^2)*w1+(d2^2)*w2+(d3^2)*w3)^0.5
      #donor group
      group <- EU_MHQ_CV$id[dd<=DIST]
      stations.Reg.group <-stations.tab[stations.tab$id%in%group,]
      EU_MHQ_CV_group <- EU_MHQ_CV[EU_MHQ_CV$id %in% group,]
      #characteristics of donor catchments
      DONORS <- data.frame(id=character(), lat=double(), lon=double(), area=double(),
                           #donor record flood before MF
                           SFOR.B=double(), SFOR.year.B=integer(), SFOR.date.B=lubridate::ymd()) 
      for(s in 1: dim(EU_MHQ_CV_group)[1]){
        DONORS[s,1] <- as.numeric(EU_MHQ_CV_group$id[s])
        DONORS[s,2] <- as.numeric(EU_MHQ_CV_group$lat[s])
        DONORS[s,3] <- as.numeric(EU_MHQ_CV_group$lon[s])
        DONORS[s,4] <- as.numeric(EU_MHQ_CV_group$area[s])
        dummy <- amax.tab[amax.tab$id==EU_MHQ_CV_group$id[s],]
        QQ <- dummy$q[dummy$year<TARGET_YEAR]
        DD <- as.Date(paste(dummy$year,dummy$month,dummy$day, sep="-"),format="%Y-%m-%d") [dummy$year<TARGET_YEAR]
        YY <- dummy$year[dummy$year<TARGET_YEAR]
        ind.SFOR <- which.max(QQ)
        if(length(ind.SFOR)>1) ind.SFOR <- max(ind.SFOR) 
        DONORS[s,5] <- QQ[ind.SFOR] 
        DONORS[s,6] <- YY[ind.SFOR] 
        DONORS[s,7] <- DD[ind.SFOR] 
        rm(dummy, ind.SFOR, QQ, YY, DD)
      }
      
      #envelope with region slope
      ENVELOPE_param$slope[i] <- env.tab$slope_qr_env[REG] 
      ENVELOPE_param$intercept.region[i] <- env.tab$intercept_qr_env[REG] 
      intercept <- -unlist(ENVELOPE_param$slope[i])*log10(DONORS$area)+log10(DONORS$SFOR.B)
      DONORS$intercept <- intercept
      DONORS <- DONORS[order(DONORS$intercept, decreasing = T),]
      DONORS_list[[i]] <- DONORS
      ENVELOPE_param$intercept.mf[i] <- max(intercept)
      ENVELOPE_param$q.predicted[i] <- ENVELOPE_param$intercept.mf[i] + ENVELOPE_param$slope[i]*log10(MEGAFLOODS_list$area[i]) 
      ENVELOPE_param$n_statyears_donors[i] <- dim(amax.tab[amax.tab$id %in% DONORS$id,])[1] # n station-years in donor group
      ENVELOPE_param$n_stat_donors[i] <- length(group) 
    }
    
    REG_MEGAFLOODS_list[[REG]] <- MEGAFLOODS_list
    REG_ENVELOPE_param[[REG]] <- ENVELOPE_param
    REG_DONORS_list[[REG]] <-  DONORS_list
  } #end loop Regions.map
  
  names(REG_MEGAFLOODS_list) <- 1:5
  names(REG_ENVELOPE_param) <- 1:5
  names(REG_DONORS_list) <- 1:5
  
  RES <- list(REG_MEGAFLOODS_list, REG_ENVELOPE_param, REG_DONORS_list)
  names(RES) <- c("mf.list", "EC.param", "donor.list")
  return(RES)
  
}
MF.predicted <- predict.mf.function(megafloods.df, stations, 2000:2022, flood_data, EU_MHQ_CV_list, ENV_reg)

PREDvsOBS.plot <- function(mf.pred.list){
  
  library(lubridate)
  library(plotrix)
  library(circular)
  
  REG_MEGAFLOODS_list <- mf.pred.list[[1]]
  REG_ENVELOPE_param <- mf.pred.list[[2]]
  REG_DONORS_list <- mf.pred.list[[3]]

  for(rr in 1:5){
    ENVELOPE_param <- REG_ENVELOPE_param[[rr]]
    MEGAFLOODS_list <- REG_MEGAFLOODS_list[[rr]]
    DONORS_list <- REG_DONORS_list[[rr]]
    

    id.mf.r <- MEGAFLOODS_list$id
    lat.r <- MEGAFLOODS_list$LAT
    lon.r <- MEGAFLOODS_list$LON
    region.r <- rep(rr, length(id.mf.r))
      
    #q
    Q.obs.r <- MEGAFLOODS_list$Q_MF
    Q.pred.r <- round(10^ENVELOPE_param$q.predicted,5)
    
    #timing
    m<-365.25
    doy.target.r <- yday(as.Date(MEGAFLOODS_list$date_MF))
    teta.target.r <- doy.target.r*2*pi/m
    teta.donors.r <- rep(NA,length(teta.target.r))
    doy.donors.r <- rep(NA,length(teta.target.r))
    R.donors.r <- rep(NA,length(teta.target.r))
    for(mm in 1:length(DONORS_list)){
      #10 donors
      DONORS <- DONORS_list[[mm]] 
      uu <-min(10,dim(DONORS)[1])
      DONORS <- DONORS[1:uu,]
      doy<-yday(as.Date(DONORS$SFOR.date.B))
      teta <- doy*2*pi/m
      x <- cos(teta)
      y <- sin(teta)
      x_donor <- mean(x, na.rm=T)
      y_donor <- mean(y, na.rm=T)
      teta.donors.r[mm] <- atan2(y_donor, x_donor)
      doy.donors.r[mm] <- round(teta.donors.r[mm]*m/(2*pi),0)  
      if(doy.donors.r[mm]<0) doy.donors.r[mm] <- round(doy.donors.r[mm]+m,0)  
      R.donors.r[mm] <- round((x_donor^2+y_donor^2)^0.5,2)
    }
    if(rr==1){
      id.mf<-id.mf.r
      lat<-lat.r
      lon<-lon.r
      region<-region.r
      Q.obs<-Q.obs.r
      Q.pred<-Q.pred.r
      doy.target<-doy.target.r
      doy.donors<-doy.donors.r
      R.donors<-R.donors.r
    }else{
      id.mf<-c(id.mf,id.mf.r)
      lat<-c(lat,lat.r)
      lon<-c(lon,lon.r)
      region<-c(region,region.r)
      Q.obs<-c(Q.obs,Q.obs.r)
      Q.pred<-c(Q.pred,Q.pred.r)
      doy.target<-c(doy.target,doy.target.r)
      doy.donors<-c(doy.donors,doy.donors.r)
      R.donors<-c(R.donors,R.donors.r)
    }

  }
  Ratio <- round(Q.obs/Q.pred,2)
  PREDvsOBS.mf <- data.frame(id.mf, lat, lon, region, Q.obs, Q.pred, Ratio, doy.target, doy.donors,R.donors)
  return(PREDvsOBS.mf)
}
PREDvsOBS.mf <- PREDvsOBS.plot(MF.predicted)




