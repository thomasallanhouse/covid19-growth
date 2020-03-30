DoubleTime <- function(dat, npts=200, plt=FALSE, figtitle=""){
  
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))#,expdt=rep(0,3))
  Tv <- seq(1,length(dat),1)
  MGAM <- gam(dat~s(Tv), family=quasipoisson)

  xv<-seq(1,length(dat), length=npts)
  newd <- data.frame(Tv=xv)
  X0 <- predict(MGAM, newd,type="lpmatrix")

  eps <- 1e-7 ## finite difference interval
  xv <- seq(1,length(dat), length=npts) +eps
  newd <- data.frame(Tv=xv)
  X1 <- predict(MGAM, newd,type="lpmatrix")

  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

  Xi <- Xp*0 
  Xi[,1:9+1] <- Xp[,1:9+1] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$sdt <- ifelse(res$sdt < 0, 0.0001, res$sdt)
  res$sdt <- log(2)/res$sdt
  res$sdtup <- ifelse(res$sdtup < 0, 0.0001, res$sdtup)
  res$sdtup <- log(2)/res$sdtup
  res$sdtlow <- ifelse(res$sdtlow < 0, 0.0001, res$sdtlow)
  res$sdtlow <- log(2)/res$sdtlow
  
  MGLM <- glm(dat~(Tv), family=quasipoisson)
  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
#  res$expdt <- Doubling
  
  if(plt==TRUE){
  par(mfrow=c(1,3))
    plot(Tv, dat, main='Fit compared with model', ylab='Cases', xlab='Time')
    lines(Tv, fitted(MGAM))
    lines(Tv, fitted(MGLM), col=2)

    p <- predict(MGAM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    lines(Tv, upr, col=1, lty=2)
    lines(Tv, lwr, col=1, lty=2)

    p <- predict(MGLM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGLM$family$linkinv(upr)
    lwr <- MGLM$family$linkinv(lwr)
    lines(Tv, upr, col=2, lty=2)
    lines(Tv, lwr, col=2, lty=2)

    plot(xv,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main=figtitle)
    lines(xv,df+2*df.sd,lty=2);
    lines(xv,df-2*df.sd,lty=2)
    abline(h=Doubling, lty=c(1,2,2), col=2)

    plot(xv,log(2)/df,type="l", ylab='Days', xlab='Time', main='Doubling time', ylim=c(0,10)) #,ylim=range(log(2)/c(df+2*df.sd,df-2*df.sd)))
    lines(xv,log(2)/(df+2*df.sd),lty=2);
    lines(xv,log(2)/ifelse(df<2*df.sd,0.00001, df-2*df.sd),lty=2)
    Doubling<-log(2)/Doubling
    abline(h=Doubling, lty=c(1,2,2), col=2)
    
#    mtext(figtitle, outer=TRUE,  cex=1, line=0)
  }
  
  res
}

library(mgcv)

# generated full analysis
# now look at main countries
##

datWHO<-data.frame(read.csv("who_data.csv"))
# convert data to incidence from raw cumulative
dat1 <- cbind(datWHO[2:79,1:2], apply(datWHO[,3:17],2,diff))
TrustVal<-c(78-35, 78) 
##
## have been setting the upper limit to length of data frame (most recent day), but earlier date is fairly arbitrary.
##

majorCountryID<-seq(3,17,1) 

pdf(file="Majordoubling.pdf")
for(i in majorCountryID){
  DoubleTime(datWHO[TrustVal[1]:TrustVal[2],i], plt=TRUE, figtitle = colnames(dat1)[i])
}
dev.off()

pdf(file="MajordoublingIncidence.pdf")
for(i in majorCountryID){
  DoubleTime(dat1[TrustVal[1]:TrustVal[2],i], plt=TRUE, figtitle = colnames(dat1)[i])
}
dev.off()


dtUK <- DoubleTime(dat1$UK[TrustVal[1]:TrustVal[2]])

datItaly<-data.frame(read.csv("italy_data.csv"))
dtffC <- DoubleTime(datItaly$Confirmed)
dtffH <- DoubleTime(c(datItaly$Hospital[1],diff(datItaly$Hospital)))
dtffI <- DoubleTime(c(datItaly$ICU[1]     ,diff(datItaly$ICU)))
dtffS <- DoubleTime(datItaly$Self.Isolated)
dtffD <- DoubleTime(datItaly$Deaths)

pdf(file="GAMItaly.pdf")
par(mfrow=c(1,1))
cols <- rainbow(4)
plot(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffC$sdt, lwd=2, type='l',xaxs='i',main="Multiple Metrics, Italy", yaxs='i', yaxp=c(1,9,8), col=cols[1], ylim=c(1,9), ylab='Doubling Time', xlab='Date')
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffC$sdtlow, lwd=2, col=cols[1], lty=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffC$sdtup, lwd=2, col=cols[1], lty=2)

lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffH$sdt, col=cols[2], lwd=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffH$sdtlow, lwd=2, col=cols[2], lty=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffH$sdtup, lwd=2, col=cols[2], lty=2)

lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffI$sdt, col=cols[3], lwd=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffI$sdtlow, lwd=2, col=cols[3], lty=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffI$sdtup, lwd=2, col=cols[3], lty=2)

lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffD$sdt, col=cols[4], lwd=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffD$sdtlow, lwd=2, col=cols[4], lty=2)
lines(as.Date(seq(datItaly$day[1],datItaly$day[33], length=200), origin = "1899-12-30"), dtffD$sdtup, lwd=2, col=cols[4], lty=2)
 axis(4, yaxp=c(1,9,8), labels=FALSE)
legend(as.Date(66, origin = "2020-01-11"), 3.5,legend=c('Confirmed','Hospitalised','ICU','Deaths'), bg='white', cex=1, lty=1, col=cols, lwd=2)
dev.off()

pdf(file="GAMMultiCountry.pdf")

dtUK <- DoubleTime(dat1$UK[TrustVal[1]:TrustVal[2]])
plot(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdt, lwd=2, type='l',xaxs='i', yaxs='i', yaxp=c(1,9,8), main="Multiple Countries, Confirmed cases", col=1, ylim=c(1,9), ylab='Doubling Time', xlab='Date')
lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdtlow, lwd=2, col=1, lty=2)
lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdtup, lwd=2, col=1,lty=2)
colval<-1
majorCountryID<-c(10, 6,7,15)
cols <- rainbow(length(majorCountryID))

# countries with major outbreaks (defined as 500 or more on last day of report) 
for(i in majorCountryID){
  dt<-DoubleTime(dat1[TrustVal[1]:TrustVal[2],i])
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdt, col=cols[colval], lwd=2)
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdtlow, lwd=2, lty=2, col=cols[colval])
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdtup, lty=2, lwd=2, col=cols[colval])
  colval=colval+1
}
legend(as.Date(77-31, origin = "2020-01-8"), 8.9,legend=c(colnames(dat1)[17],colnames(dat1)[majorCountryID]), bg='white', cex=1, lty=1, col=c(1,cols), lwd=2)
 axis(4, yaxp=c(1,9,8), labels=FALSE)
dev.off()

pdf(file="Italydoubling.pdf")
dtffC <- DoubleTime(datItaly$Confirmed, plt=TRUE, figtitle = 'Confirmed Cases')
dtffH <- DoubleTime(c(datItaly$Hospital[1],diff(datItaly$Hospital)), plt=TRUE, figtitle = 'Hospitalisation') #colnames(datItaly)[i])
dtffI <- DoubleTime(c(datItaly$ICU[1]     ,diff(datItaly$ICU)), plt=TRUE, figtitle = 'ICU')# colnames(datItaly)[i])
dtffD <- DoubleTime(datItaly$Deaths, plt=TRUE, figtitle = 'Deaths')
#for(i in 2:6){
#  DoubleTime(datItaly[,i], plt=TRUE, figtitle = colnames(datItaly)[i])
#}
dev.off()

pdf(file="GAMMultiCountryCumulative.pdf")

dtUK <- DoubleTime(datWHO$UK[TrustVal[1]:TrustVal[2]])
plot(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdt, lwd=2, type='l',xaxs='i', yaxs='i', yaxp=c(1,9,8), main="Multiple Countries, Confirmed cases", col=1, ylim=c(1,9), ylab='Doubling Time', xlab='Date')
lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdtlow, lwd=2, col=1, lty=2)
lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dtUK$sdtup, lwd=2, col=1,lty=2)
colval<-1
majorCountryID<-c(10, 6,7,15)
cols <- rainbow(length(majorCountryID))

for(i in majorCountryID){
  dt<-DoubleTime(datWHO[TrustVal[1]:TrustVal[2],i])
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdt, col=cols[colval], lwd=2)
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdtlow, lwd=2, lty=2, col=cols[colval])
  lines(as.Date(seq(TrustVal[1]+1,TrustVal[2], length=200), origin = "2020-01-9"), dt$sdtup, lty=2, lwd=2, col=cols[colval])
  colval=colval+1
}
legend(as.Date(77-31, origin = "2020-01-8"), 8.9,legend=c(colnames(dat1)[17],colnames(dat1)[majorCountryID]), bg='white', cex=1, lty=1, col=c(1,cols), lwd=2)
 axis(4, yaxp=c(1,9,8), labels=FALSE)
dev.off()

