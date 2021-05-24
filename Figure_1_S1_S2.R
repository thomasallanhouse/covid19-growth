## Copyright (c) 2020 Ian Hall
## See LICENCE for licensing information

## Growth rate estimates for confirmed cases in Europe and for different metrics in Italy using GAM
## Figure 1 (main text) and figures S1 and S2 (electronic supplementary material) of:
## 
## Pellis L, Scarabel F, Stage HB, Overton CE, Chappell LHK, Fearon E, Bennett E, 
## University of Manchester COVID-19 Modelling Group, Lythgoe KA, House TA and Hall I, 
## "Challenges in control of COVID-19: short doubling time and long delay to effect of interventions", 
## Philosophical Transactions of the Royal Society B (2021)


DoubleTime <- function(dat, timev, npts=200, meth="GCV.Cp", FE='None', plt=FALSE, thck=1.5, figtitle="", origin="2020-01-09"){
  
  res <- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts),doub=rep(0,npts),doubup=rep(0,npts),doublow=rep(0,npts))
  Tv <- timev
  
  if(FE=='None'){
    MGAM <- gam(dat~s(Tv), family=quasipoisson, method=meth)
  }else{
  DW <- weekdays(as.Date(Tv, origin = origin))
  if(FE=='WE'){
    DW <- ifelse(DW=='Sunday','WE',ifelse(DW=='Saturday','WE','WD'))
  }
  MGAM <- gam(dat~s(Tv)+DW, family=quasipoisson, method=meth)
  }
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  if(FE=='None'){
    newd <- data.frame(Tv=xv)
  }else{
    dow <- weekdays(as.Date(xv, origin = origin))
    if(FE=='WE'){
      dow <- ifelse(dow=='Sunday','WE',ifelse(dow=='Saturday','WE','WD'))
    }
    newd <- data.frame(Tv=xv, DW=dow)
  }
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  if(FE=='None'){
    newd <- data.frame(Tv=xv)
  }else{
    dow <- weekdays(as.Date(xv, origin = origin))
    if(FE=='WE'){
      dow <- ifelse(dow=='Sunday','WE',ifelse(dow=='Saturday','WE','WD'))
    }
    newd <- data.frame(Tv=xv, DW=dow)
  }
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

  off <- ifelse(FE=='None',1,ifelse(FE=='WE',2,7))  
  Xi <- Xp*0 
  Xi[,1:9+off] <- Xp[,1:9+off] ## weekend Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  ## derivative calculation, pers comm S. N. Wood, found in mgcv:  Mixed  GAM  Computation  Vehicle  with  automatic  smoothness  estimation.  R  packageversion 1.8-31 (2019) https://CRAN.R-project.org/package=mgcv.

  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$doub <- ifelse(res$sdt < 0, 100, log(2)/res$sdt)
  res$doubup <- ifelse(res$sdtup < 0, 100, log(2)/res$sdtup)
  res$doublow <- ifelse(res$sdtlow < 0, 100, log(2)/res$sdtlow)

  MGLM <- glm(dat~(Tv), family=quasipoisson)
  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])

  if(plt==TRUE){
    par(mar = c(5,4,4,4) + 0.1)
    par(mfrow=c(1,2))
    plot(as.Date(xv, origin = origin),df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main=figtitle, lwd=2*thck)
    lines(as.Date(xv, origin = origin),df+2*df.sd,lty=2, lwd=thck);
    lines(as.Date(xv, origin = origin),df-2*df.sd,lty=2, lwd=thck)
    abline(h=0, col=4)
    text(as.Date(43930, origin = origin),max(df),ifelse(df[length(df)]-2*df.sd[length(df)]>0,'Growth',ifelse(df[length(df)]+2*df.sd[length(df)]<0,'Decay','Plateau')))
    axis(4,at=c(-log(2)/7, 0,log(2)/7,log(2)/4,log(2)/3,log(2)/2), labels=c(-7,'Infinite',7,4,3,2))
    mtext("Doubling time", side = 4, line = 3)
    
    plot(as.Date(Tv, origin = origin), dat, main='Fit compared with model', ylab='Number', xlab='Time', pch=16, col=4)
    lines(as.Date(Tv, origin = origin), fitted(MGAM), lwd=2*thck)

    p <- predict(MGAM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    lines(as.Date(Tv, origin = origin), upr, col=1, lty=2, lwd=thck)
    lines(as.Date(Tv, origin = origin), lwr, col=1, lty=2, lwd=thck)
    
  }
  
  res
}

library(mgcv)

## Data:
# The files in the repository contain time series of cumulative numbers, both for confirmed cases in Europe and for different metrics in Italy.
# Under the assumption that data reported on day n are those collected on day n-1, the corresponding dates have been shifted backwards by one day
# Example: to plot the number of daily cases published by WHO on 1 Mar, one needs to read the cumulative data in this dataset from 28 Feb (2 days earlier);
# then use "diff" to compute the new cases collected on 29 Feb, which are plotted against 29 Feb (x-axes in plots are times of real-collection)
# For reference: Italy published latest data at 6pm, so there is one extra day of data compared to WHO data file, 
# so Italy has collected data on 31 March and published them on 31 March, but WHO will publish it on 1 Apr.

## European countries:
start <- 52 # start = 52 is 29 Feb on both cumulative and daily incidence: this is the data published by WHO on 1 Mar. 
# To compute the incidence on day 52 (29 Feb), the file is read from day 51 (the "start-1" appearing in lines 127 and 151)
Lag <- 0 # Use data till end of dataset

## Italy:
startIT <- 2 # startIT = 2 is 24 Feb on Italian data
# As above, startIT cannot be 1, because incidence on day startIT requires computing the difference between startIT and startIT-1 (lines 135 and 172)
LagIT <- 15 # Plot only 3 weeks for Italy (discard second half of March), because hospital/ICU data are prevalence, and incidence cannot be reconstructed. 
# We still calculate daily counts, but the approximation is sensible only on a short enough period that we can assume not many individuals have left hospital/ICU.

## Open data files
datD <-data.frame(read.csv(paste("WHO_data_18.csv", sep=''))) 
datITD<-data.frame(read.csv(paste("Italy.csv", sep='')))
head(datD)
head(datITD)
ncountries = ncol(datD)-2

# Default methods as used in the paper (GCV.Cp with day-of-the-week effect, WD) 
# Further possibilities available for FitMethod, FixedEff
FitMethod <-c('GCV.Cp') # c('GCV.Cp','ML') # Explore fit with different options for the GAM
FixedEff <-c('WD') # c('WE', 'WD', 'None') # Explore alternative fit with weekend effect (WE), day-of-the-week effect (WD), or nothing
colval <-c('black', 'red', 'green', 'blue', 'purple')
clist <- c(18,5, 6, 8, 15) # List of the columns for the 5 countries in Figure 1(a)
clistIT <- c(1:3,5) # Relevant columns from Italian data for Figure 1(b)
thin <- 2 # Width of thin line
thick <- 3 # Width of thick line

print('Generating figures for the Electronic Supplementary Material:')
for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('SuppWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
      for(i in 1:ncountries){
          DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+i]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datD)[2+i])
      }
    dev.off()
    
    pdf(file=paste('SuppItaly_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    for(i in 1:5){
      if(i!=4){
        DoubleTime(diff(datITD[(startIT-1):(dim(datITD)[1]-LagIT),2+i]),datITD$Day[startIT:(dim(datITD)[1]-LagIT)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datITD)[2+i], origin="2020-02-23")
      }
    }
    dev.off()
  }
}

print('Generating figures for the Main text:')
for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('MainDoublingTimeWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    par(mfrow=c(1,1))
    par(mar = c(5,4,4,4) + 0.1)
    plot(as.Date(seq(datD$Day[start], datD$Day[(dim(datD)[1]-Lag)], 1), origin="2020-01-09"), 
         rep(0,(1+dim(datD)[1]-Lag-start)),ylab='Instantaneous doubling time', xlab='Date', 
         ylim=c(1.5,7), col='white', main='Confirmed cases - Multiple countries')
    for(i in 1:length(clist)){
      temp<-DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+clist[i]]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j])
      lines(as.Date(seq(datD$Day[start], datD$Day[(dim(datD)[1]-Lag)], length.out=200), origin="2020-01-09"), 
            temp[,4], col=colval[i], lwd=thick)
      lines(as.Date(seq(datD$Day[start], datD$Day[(dim(datD)[1]-Lag)], length.out=200), origin="2020-01-09"), 
            temp[,5], lty=2, col=colval[i], lwd=thin)
      lines(as.Date(seq(datD$Day[start], datD$Day[(dim(datD)[1]-Lag)], length.out=200), origin="2020-01-09"), 
            temp[,6], lty=2, col=colval[i], lwd=thin)
    }
    mtext("Instantaneous growth rate", side = 4, line = 3)
    axis(4,at=log(2)/seq(0.1,0.4,by=0.05), labels=seq(0.1,0.4,by=0.05),ylab ='Instantaneous growth rate')
    legend( as.Date(datD$Day[(dim(datD)[1]-Lag)-8], origin="2020-01-09"), 2.85,col=colval, lty=1, lwd=thick, legend=colnames(datD)[2+clist])
    dev.off()
    
    pdf(file=paste('MainDoublingTimeItaly_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    par(mfrow=c(1,1))
    par(mar = c(5,4,4,4) + 0.1)
    plot(as.Date(seq(datITD$Day[startIT], datITD$Day[(dim(datITD)[1]-LagIT)], 1), origin="2020-02-23"), 
         rep(0,(1+dim(datITD)[1]-LagIT-startIT)),ylab='Instantaneous doubling time', xlab='Date', 
         ylim=c(1.5,7), col='white', main='Multiple data streams - Italy')
    for(i in 1:length(clistIT)){
      temp<-DoubleTime(diff(datITD[(startIT-1):(dim(datITD)[1]-LagIT),2+clistIT[i]]),datITD$Day[startIT:(dim(datITD)[1]-LagIT)], meth=FitMethod[k], FE=FixedEff[j], origin="2020-02-23")
      lines(as.Date(seq(datITD$Day[startIT], datITD$Day[(dim(datITD)[1]-LagIT)], length.out=200), origin="2020-02-23"), 
            temp[,4], col=colval[i], lwd=thick)
      lines(as.Date(seq(datITD$Day[startIT], datITD$Day[(dim(datITD)[1]-LagIT)], length.out=200), origin="2020-02-23"), 
            temp[,5], lty=2, col=colval[i], lwd=thin)
      lines(as.Date(seq(datITD$Day[startIT], datITD$Day[(dim(datITD)[1]-LagIT)], length.out=200), origin="2020-02-23"), 
            temp[,6], lty=2, col=colval[i], lwd=thin)
    }
    mtext("Instantaneous growth rate", side = 4, line = 3)
    axis(4,at=log(2)/seq(0.1,0.4,by=0.05), labels=seq(0.1,0.4,by=0.05),ylab ='Instantaneous growth rate')
    legend( as.Date(datITD$Day[(dim(datITD)[1]-LagIT)-5], origin="2020-02-23"), 2.56,col=colval, lty=1, lwd=thick, legend=colnames(datITD)[2+clistIT])
    dev.off()
  }
}


