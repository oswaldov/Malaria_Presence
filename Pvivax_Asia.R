## Load R packages

library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(quantmod)
library(stargazer)
library(rpart.plot)

## cleaning up the malaria-climate data file
pviv<-read.csv("Prevalence_Pvivax.csv", header=TRUE)

## subset by continent
dat<-pviv[which(pviv$continent=="Asia"),]

## select columns of interest
dd<-dat[c(1,4,6,8,10,12:116)]


## Only prevalence cases between 2 months from the start time.
ddd <- dd %>%
  mutate(difference = ifelse(year_end == year_start, month_end - month_start,
                             ifelse(year_end != year_start, ((0 + month_end)+(12 - month_start)) , NA)))


dat1 <- ddd [which(ddd$difference <= 2),]

dat1<-dat1[dat1$year_start>=1990,]

##Remove NA values

ww<- which(!is.na(dat1$pv_pos))
dat1<- dat1[ww,]

ww<- which(!is.na(dat1$gdpvl))
dat1<- dat1[ww,]

ww<- which(!is.na(dat1$pop_den))
dat1<- dat1[ww,]

ww<- which(!is.na(dat1$NDVI_index))
dat1<- dat1[ww,]

dat1<- dat1[ww,]

## summary of presence/absence data
table(dat1$pv_pa)

## log transformation pop_den and GDP

dat1$lGDP<- round(log(dat1$gdpvl), digits = 2)


dat1$lpop_den <- round(log(dat1$pop_den), digits = 2) + 5

## convert NDVI data between 0 and 1
dat1$NDVI <- (dat1$NDVI_byte * 0.004)

## select variables of interests
preds<-c("pv_pa")
age<- c("lower_age","upper_age")
tempmin<- c("tmin1","tmin2","tmin3","tmin4","tmin5","tmin6",
            "tmin7","tmin8","tmin9","tmin10","tmin11","tmin12")
tempmean<- c("tmean1","tmean2","tmean3","tmean4","tmean5","tmean6",
             "tmean7","tmean8","tmean9","tmean10","tmean11","tmean12")
tempmax<- c("tmax1","tmax2","tmax3","tmax4","tmax5","tmax6",
            "tmax7","tmax8","tmax9","tmax10","tmax11","tmax12")
precip<- c("precp1","precp2","precp3","precp4","precp5","precp6",
           "precp7","precp8","precp9","precp10","precp11","precp12")
qtempmin<- c("qtmin1","qtmin2","qtmin3","qtmin4")
qtempmean<- c("qtmean1","qtmean2")##,"qtmean3","qtmean4")
qtempmax<- c("qtmax1","qtmax2","qtmax3","qtmax4")
qprecip<- c("qprecp1","qprecp2")##,"qprecp3","qprecp4")
bio<- c("bio3","bio16","bio17") ##"bio15"
ro<- c("Roq1","Roq2")##,"Roq3","Roq4")
other<- c("year_start", "elev","NDVI")
soc<- c("hindex","lGDP","lpop_den")

all<- c(preds, qtempmean, qprecip, bio, ro, other, soc)

datapvasi<- dat1[,all]

##Build CART trees

## First fit forcing the tree to have lots of branches so we can
## examine it and figure out where to trim
CP=0.0005
MS=50
cnt<-rpart.control(minsplit=MS, cp=CP, xval=100)

f.null<-rpart(pv_pa ~  . , data=datapvasi, method="class", control=cnt)
plotcp(f.null) ## use this to decide on a cp/size for trimming the tree
printcp(f.null)
graphics.off()

x11()
par(mfrow=c(1,1))
#plot(f.null, uniform=TRUE, margin=0.0075)
#text(f.null, digits=1, use.n=TRUE, cex=0.5)

rpart.plot(f.null, type=4, extra=104, box.palette = "GnBu",branch.lty=3,shadow.col="gray",nn=TRUE)

f.null$cptable


mm<-min(signif(f.null$cptable[,4],3))+mean(f.null$cptable[,5])
mm

w<-which(f.null$cptable[,4]<mm)
w

CP.new<-f.null$cptable[min(w)-1,1]
CP.new


## looking for 1 standard error from 

## now trimming based on the above
cnt<-rpart.control(minsplit= MS, cp=CP.new, xval=100)
f<-rpart(pv_pa ~  . , data=datapvasi, method = "class", control=cnt)
plotcp(f)
graphics.off()

x11()
par(mfrow=c(1,1))
#plot(f, uniform=TRUE, margin=0.0075); text(f, digits=1, use.n=TRUE, cex=0.5)
rpart.plot(f, type=4, extra=104, box.palette = "GnBu",branch.lty=3,shadow.col="gray",nn=TRUE, 
           space = 0, tweak = 0.8)

## GLM models 
## squared bio variables

datapvasi$bio03sq<-datapvasi$bio3^2
datapvasi$qtmean1sq<- datapvasi$qtmean1^2
datapvasi$qprecp1sq<- datapvasi$qprecp1^2


head(datapvasi)
dim(datapvasi)
temps2<-c("bio03sq","qtmean1sq","qprecp1sq")
par(mfrow=c(2,2), bty="n")
for(i in temps2) plot(datapvasi[,i], datapvasi$pv_pa, main=i)

##models

m.null<-glm(pv_pa ~ 1, family="binomial", data=datapvasi)
m.null
m.full<-glm(pv_pa ~ . , family="binomial", data=datapvasi)
m.full


modelpvAs<- step(m.null, scope=formula(m.full), direction="forward", criterion = "BIC")
summary(modelpvAs)


##Plot residuals
plot(fitted(modelpvAS), res)
abline(0,0)

par(mfrow=c(1,1))
qqnorm(res, main = expression('Q-Q plot for '*italic(P.~vivax)*' in Asia'), bty = "o")
qqline(res, col = "steelblue", lwd = 3)
legend(-4.5, 3.9, legend = "E", bty = "n", cex = 2.5)

## Plot density 
plot(density(res, adjust = 1.2), main = expression('Density plot for '*italic(P.~vivax)*' in Asia'), lwd = 2, bty = "o")
legend(-4.8, 0.43, legend = "F", bty = "n", cex = 2.5)


## Build prediction curves using outputs from GLM and CART models
## lets look at marginal predictions, based on particular predictors

## for bio3: Isothermality

ey <- expression(italic(P. ~ vivax) ~ predicted ~ presence)

o<-order(datapvasi$bio3)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd = TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$bio3[o], datapvasi$pv_pa, xlab="Temperature oscillation (%)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(10,100),main="GLM fit bio3-Asia")
fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio3=datapvasi$bio3[o])
y.m2<-predict(loess(fit~bio3, data=fit.dat.m, span = 1), data.frame(bio3=seq(10, 100, by=1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio03[o]/10, y.sym, col=2, lwd=3)
lines(seq(10, 100, by=1), y.m2$fit, col=2, lwd=3)
legend(10, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$bio3[o]), datapvasi$pv_pa, ylim=c(0,1),xlim=c(10,100),
     xlab="Temperature oscillation (%)", ylab= ey, cex.lab = 1.1,
     main="CART fit bio3-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio3=datapvasi$bio3[o])
fit.dat.f
y.f2<-predict(loess(fit~bio3, data=fit.dat.f, span = 1), data.frame(bio3=seq(10, 100, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(10, 100, by=1), y.f2$fit, col=4, lwd=3)
legend(10, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=105, y=1.3, legend="E", xpd = NA, bty="n", cex = 1.7)



## for bio16: Min temperature in coldest month (C)

o<-order(datapvasi$bio16)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(datapvasi$bio16[o], datapvasi$pv_pa, xlab="Precipitation wettest qtr (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,3000),main="GLM fit bio16-Asia")
fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio16=datapvasi$bio16[o])
y.m2<-predict(loess(fit~bio16, data=fit.dat.m, span = 1), data.frame(bio16=seq(0, 3000, by=10)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,3000, by=10), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$bio16[o]), datapvasi$pv_pa, ylim=c(0,1),xlim=c(0,3000),
     xlab="Precipitation wettest qtr (mm)", ylab= ey, cex.lab = 1.1,
     main="CART fit bio16-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio16=datapvasi$bio16[o])
fit.dat.f
y.f2<-predict(loess(fit~bio16, data=fit.dat.f, span=1), data.frame(bio16=seq(0, 3000, by=10)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 3000, by=10), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=3050, y=1.3, legend="F", xpd = NA, bty="n", cex = 1.7)


##bio17
## Precipitation of driest quarter

o<-order(datapvasi$bio17)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$bio17[o], datapvasi$pv_pa, xlab="Precipitation driest qtr (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,1000),main="GLM fit bio17-Asia")
fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio17=datapvasi$bio17[o])
y.m2<-predict(loess(fit~bio17, data=fit.dat.m, span=1), data.frame(bio17=seq(0, 1000, by=10)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,1000, by=10), y.m2$fit, col=2, lwd=3)
legend(10, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$bio17[o]), datapvasi$pv_pa,ylim=c(0,1),xlim=c(0,1000),
     xlab="Precipitation driest qtr (mm)", ylab= ey, cex.lab = 1.1, 
     main="CART fit bio17-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio17=datapvasi$bio17[o])
fit.dat.f
y.f2<-predict(loess(fit~bio17, data=fit.dat.f, span=1), data.frame(bio17=seq(0, 1000, by=10)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1000, by=10), y.f2$fit, col=4, lwd=3)
legend(10, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1000, y=1.3, legend="G", xpd = NA, bty="n", cex = 1.7)


##FOR qtmean1 ##########

o<-order(datapvasi$qtmean1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$qtmean1[o], datapvasi$pv_pa, xlab="Temperature (°C)",
     ylab= ey, cex.lab = 1.1,
     xlim=c(0,35), ylim=c(0,1),main="GLM fit qtmean1-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qtmean1=datapvasi$qtmean1[o])
y.m2<-predict(loess(fit~qtmean1, data=fit.dat.m, span=1), data.frame(qtmean1=seq(0, 35, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit <0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 35, by=0.1), y.m2$fit, col=2, lwd=3)
legend(2, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$qtmean1[o]), datapvasi$pv_pa,
     xlab="Temperature (°C)", ylab= ey, cex.lab = 1.1, 
     xlim=c(0,35), ylim=c(0,1),main="CART fit qtmean1-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qtmean1=datapvasi$qtmean1[o])
fit.dat.f
y.f2<-predict(loess(fit~qtmean1, data=fit.dat.f, span=1), data.frame(qtmean1=seq(0, 35, by=0.1)), se=TRUE)
y.f2
y.f2$fit[y.f2$fit <0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 35, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(2, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

legend(x=37, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)


##FOR qtmean2

o<-order(datapvasi$qtmean2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$qtmean2[o],datapvasi$pv_pa, xlab="Temperature (°C)",
     ylab= ey, cex.lab = 1.1, 
     xlim=c(-10,30), ylim=c(0,1),main="GLM fit qtmean2-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qtmean2=datapvasi$qtmean2[o])
y.m2<-predict(loess(fit~qtmean2, data=fit.dat.m, span=1), data.frame(qtmean2=seq(-10, 30, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit < 0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(-10,30, by=0.1), y.m2$fit, col=2, lwd=3)
legend(-5, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$qtmean2[o]), datapvasi$pv_pa,
     xlab="Temperature (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(-10,30), ylim=c(0,1),main="CART fit qtmean2-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qtmean2=datapvasi$qtmean2[o])
fit.dat.f
y.f2<-predict(loess(fit~qtmean2, data=fit.dat.f, span=1), data.frame(qtmean2=seq(-10,30, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit < 0] = 0
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(-10, 30, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(-5, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x= 32, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)


###FOR qprecp1

o<-order(datapvasi$qprecp1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$qprecp1[o], datapvasi$pv_pa[o], xlab="Precipitation (mm)",
     ylab= ey, cex.lab=1.1,  xlim=c(0,1000),ylim=c(0,1),
     main="GLM fit qprecp1-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qprecp1=datapvasi$qprecp1[o])
y.m2<-predict(loess(fit~qprecp1, data=fit.dat.m, span=1), data.frame(qprecp1=seq(0, 1000, by=1)), se=TRUE)
y.m2$fit[y.m2$fit < 0] <- 0

lines(seq(0, 1000, by=1), y.m2$fit, col=2, lwd=3)
legend(50, 0.9, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$qprecp1[o]), datapvasi$pv_pa[o],
     xlab="Precipitation (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1000), ylim=c(0,1),main="CART fit qprecp1-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qprecp1=datapvasi$qprecp1[o])
fit.dat.f
y.f2<-predict(loess(fit~qprecp1, data=fit.dat.f, span=1), data.frame(qprecp1=seq(0, 1000, by=1)), se=TRUE)
y.f2
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(0, 1000, by=1), y.f2$fit, col=4, lwd=3)
legend(50, 0.9, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1050, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)


###FOR qprecp2

o<-order(datapvasi$qprecp2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$qprecp2[o], datapvasi$pv_pa[o], xlab="Precipitation (mm)",
     ylab= ey, cex.lab=1.1, xlim=c(0,700),ylim=c(0,1),
     main="GLM fit qprecp2-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qprecp2=datapvasi$qprecp2[o])
y.m2<-predict(loess(fit~qprecp2, data=fit.dat.m, span=0.96), data.frame(qprecp2=seq(0, 700, by=1)), se=TRUE)
y.m2

lines(seq(0, 700, by=1), y.m2$fit, col=2, lwd=3)
legend(50, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$qprecp2[o]), datapvasi$pv_pa[o],
     xlab="Precipitation (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,700), ylim=c(0,1),main="CART fit qprecp2-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qprecp2=datapvasi$qprecp2[o])
fit.dat.f
y.f2<-predict(loess(fit~qprecp2, data=fit.dat.f,span=0.95), data.frame(qprecp2=seq(0, 700, by=1)), se=TRUE)
y.f2
##lines
lines(seq(0, 700, by=1), y.f2$fit, col=4, lwd=3)
legend(50, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=750, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)


### Elevation (DEM)

o<-order(datapvasi$elev)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$elev[o], datapvasi$pv_pa[o], xlab="Elevation (meters)",
     ylab= ey, cex.lab=1.1, xlim=c(0,4000),ylim=c(0,1),
     main="GLM fit elevation-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, elev=datapvasi$elev[o])
y.m2<-predict(loess(fit~elev, data=fit.dat.m, span=1), data.frame(elev=seq(0,4000, by=10)), se=TRUE)
y.m2

lines(seq(0, 4000, by=10), y.m2$fit, col=2, lwd=3)
legend(0, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$elev[o]), datapvasi$pv_pa[o],
     xlab="Elevation (meters)", ylab= ey, cex.lab=1.1,
     xlim=c(0,4000), ylim=c(0,1),main="CART fit Elevation-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, elev=datapvasi$elev[o])
fit.dat.f
y.f2<-predict(loess(fit~elev, data=fit.dat.f, span=1), data.frame(elev=seq(0,4000, by=10)), se=TRUE)
y.f2
##lines
lines(seq(0, 4000, by=10), y.f2$fit, col=4, lwd=3)
legend(0, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=4150, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)

###NDVI

o<-order(datapvasi$NDVI)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$NDVI[o], datapvasi$pv_pa[o], xlab="NDVI index",
     ylab= ey, cex.lab=1.1, xlim=c(0,1),ylim=c(0,1),
     main="GLM fit NDVI-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, NDVI=datapvasi$NDVI[o])
y.m2<-predict(loess(fit~NDVI, data=fit.dat.m, span=1), data.frame(NDVI=seq(0,1, by=0.1)), se=TRUE)
y.m2

lines(seq(0, 1, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.1, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$NDVI[o]), datapvasi$pv_pa[o],
     xlab="NDVI index", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit NDVI-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, NDVI=datapvasi$NDVI[o])
fit.dat.f
y.f2<-predict(loess(fit~NDVI, data=fit.dat.f, span=1), data.frame(NDVI=seq(0,1, by=0.1)), se=TRUE)
y.f2
##lines
lines(seq(0, 1, by=0.1), y.f2$fit, col=4, lwd=3)
legend(0.1, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1, y=1.3, legend="D", xpd = NA, bty="n", cex = 1.7)



## socio -demographic variables 
###FOR YEAR START 

o<-order(datapvasi$year_start)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$year_start[o], datapvasi$pv_pa[o], xlab="Year of the study",
     ylab= ey, cex.lab=1.1, xlim=c(1990,2015),ylim=c(0,1),
     main="GLM fit year_study-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, year_start=datapvasi$year_start[o])
y.m2<-predict(loess(fit~year_start, data=fit.dat.m, span=1), data.frame(year_start=seq(1990,2015, by=1)), se=TRUE)
y.m2

lines(seq(1990, 2015, by=1), y.m2$fit, col=2, lwd=3)
legend(1990, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$year_start[o]), datapvasi$pv_pa[o],
     xlab="Year of the study", ylab= ey, cex.lab=1.1,
     xlim=c(1990,2015), ylim=c(0,1),main="CART fit year_study-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, year_start=datapvasi$year_start[o])
fit.dat.f
y.f2<-predict(loess(fit~year_start, data=fit.dat.f, span=1), data.frame(year_start=seq(1990, 2015, by=1)), se=TRUE)
y.f2$fit[y.f2$fit < 0] <- 0
##lines
lines(seq(1990, 2015, by=1), y.f2$fit, col=4, lwd=3)
legend(1990, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=2017, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)


###FOR LGDP ADJUSTED 

o<-order(datapvasi$lGDP)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$lGDP[o], datapvasi$pv_pa[o], xlab="log(GDP) per capita",
     ylab= ey, cex.lab=1.1, xlim=c(7,10),ylim=c(0,1),
     main="GLM fit GDPPC-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, lGDP=datapvasi$lGDP[o])
y.m2<-predict(loess(fit~lGDP, data=fit.dat.m, span=1), data.frame(lGDP=seq(7,10, by=0.05)), se=TRUE)
y.m2

##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(7, 10, by=0.05), y.m2$fit, col=2, lwd=3)
legend(7, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$lGDP[o]), datapvasi$pv_pa[o],
     xlab="log(GDP) per capita", ylab= ey, cex.lab=1.1, 
     xlim=c(7,10), ylim=c(0,1),main="CART fit GDPPC-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, lGDP=datapvasi$lGDP[o])
fit.dat.f
y.f2<-predict(loess(fit~lGDP, data=fit.dat.f, span=1), data.frame(lGDP=seq(7,10, by=0.05)), se=TRUE)
y.f2$fit[y.f2$fit < 0] = 0
##lines
lines(seq(7, 10, by=0.05), y.f2$fit, col=4, lwd=3)
legend(7, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=10.2, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)


###FOR population density####

o<-order(datapvasi$lpop_den)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

x <- datapvasi$lpop_den[o] - 5
x[x < 0] <- 0
x

#datapvasi$lpop_den[o]-5

## plot based off predictions from the glm for data pre-2004
plot(x,datapvasi$pv_pa[o], xlab=expression(log~population~density~(Km^2)),
     ylab= ey, cex.lab= 1.1, xlim=c(0,10),ylim=c(0,1),
     main="GLM fit pop density-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, pop_den=x)
y.m2<-predict(loess(fit~pop_den, data=fit.dat.m, span=1), data.frame(pop_den=seq(0, 10, by=0.1)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 10, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.5, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(x), datapvasi$pv_pa[o],
     xlab=expression(log~population~density~(Km^2)), ylab= ey, cex.lab=1.1, 
     xlim=c(0,10), ylim=c(0,1),main="CART fit pop density-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, pop_den=x)
fit.dat.f
y.f2<-predict(loess(fit~pop_den, data=fit.dat.f, span=1), data.frame(pop_den=seq(0, 10, by=0.1)), se=TRUE)
y.f2
y.f2$fit[y.f2$fit > 0.98] = 0.98

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 10, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(0.5, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=10.7, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)

###FOR Human Index####

o<-order(datapvasi$hindex)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$hindex[o],datapvasi$pv_pa[o], xlab="Human development index",
     ylab= ey, cex.lab= 1.1, xlim=c(0.3,0.8),ylim=c(0,1),
     main="GLM fit hindex-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, hindex=datapvasi$hindex[o])
y.m2<-predict(loess(fit~hindex, data=fit.dat.m, span=1), data.frame(hindex=seq(0.3, 0.8, by=0.05)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0.3, 0.8, by=0.05), y.m2$fit, col=2, lwd=3)
legend(0.29, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$hindex[o]), datapvasi$pv_pa[o],
     xlab="Human development index", ylab= ey, cex.lab=1.1, 
     xlim=c(0.3,0.8), ylim=c(0,1),main="CART fit hindex-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, hindex=datapvasi$hindex[o])
fit.dat.f
y.f2<-predict(loess(fit~hindex, data=fit.dat.f, span=1), data.frame(hindex=seq(0.3, 0.8, by=0.05)), se=TRUE)
y.f2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0.3, 0.8, by=0.05), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(0.29, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=0.83, y=1.3, legend="D", xpd = NA, bty="n", cex = 1.7)


###FOR R0q1

o<-order(datapvasi$Roq1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$Roq1[o], datapvasi$pv_pa[o], xlab="R0-Temp 1st quarter",
     ylab= ey, cex.lab=1.1,  xlim=c(0,1),ylim=c(0,1),
     main="GLM fit R0q1-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, Roq1=datapvasi$Roq1[o])
y.m2<-predict(loess(fit~Roq1, data=fit.dat.m, span=1), data.frame(Roq1=seq(0,1, by=0.1)), se=TRUE)
y.m2

lines(seq(0, 1, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.05, 0.26, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 

## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$Roq1[o]), datapvasi$pv_pa[o],
     xlab="R0-Temp 1st quarter", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit R0q1-Asia")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, Roq1=datapvasi$Roq1[o])
fit.dat.f
y.f2<-predict(loess(fit~Roq1, data=fit.dat.f, span=1), data.frame(Roq1=seq(0,1, by=0.1)), se=TRUE)
y.f2
##lines
lines(seq(0, 1, by=0.1), y.f2$fit, col=4, lwd=3)
legend(0.05, 0.26, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1.05, y=1.3, legend="E", xpd = NA, bty="n", cex = 1.7)


######FOR ROq2 #####

o<-order(datapvasi$Roq2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd = TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapvasi$Roq2[o], datapvasi$pv_pa[o], xlab="R0-Temp 2nd quarter",
     ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="GLM fit R0q2-Asia")

fitted<-as.numeric(modelpvAs$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, Roq2=datapvasi$Roq2[o])
y.m2<-predict(loess(fit~Roq2, data=fit.dat.m, span=1), data.frame(Roq2=seq(0,1, by=0.1)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.05, 0.26, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapvasi$Roq2[o]), datapvasi$pv_pa[o],
     xlab="R0-Temp 2nd quarter", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit R0q2-Asia")

fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, Roq2=datapvasi$Roq2[o])
fit.dat.f
y.f2<-predict(loess(fit~Roq2, data=fit.dat.f, span=1.3), data.frame(Roq2=seq(0, 1, by=0.1)), se=TRUE)
y.f2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1, by=0.1), y.f2$fit, col=4, lwd=3)
legend(0.05, 0.26, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1.05, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)



