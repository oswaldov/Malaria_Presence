

## Load the R packages
library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(quantmod)
library(rpart.plot)

setwd()
## Load the malaria-presence data file
pfalcAF<-read.csv("Prevalence_Pfalciparum.csv", header=TRUE)

## subset by continent
datpf<-pfalcAF[which(pfalcAF$continent=="Africa"),]

## select columns of interest
dd<-datpf[c(1,4,6,8,10:115)]


## Only prevalence cases between 2 months from the start time.
ddd <- dd %>%
  mutate(difference = ifelse(year_end == year_start, month_end - month_start,
                             ifelse(year_end != year_start, ((0 + month_end)+(12 - month_start)) , NA)))


dat1 <- ddd [which(ddd$difference <= 2),]

dat1<-dat1[dat1$year_start>=1990,]

### check for NA values

ww<- which(!is.na(dat1$pf_pos))
dat1<- dat1[ww,]
dim(dat1)

ww<- which(!is.na(dat1$gdpvl))
dat1<- dat1[ww,]
dim(dat1)

ww<- which(!is.na(dat1$pop_den))
dat1<- dat1[ww,]
dim(dat1)

ww<- which(!is.na(dat1$NDVI_byte))
dat1<- dat1[ww,]
dim(dat1)

dat1<- dat1[ww,]
head(dat1)

## summary of presence/absence data
table(dat1$pf_pa)

## log transformation pop_den and GDP

dat1$lpop_den <- log(dat1$pop_den) + 4
dat1$lpop_den <- round(dat1$lpop_den, 2)

dat1$lgdppc <- round(log(dat1$gdpvl), 2)

## convert NDVI data between 0 and 1

dat1$NDVI <- (dat1$NDVI_byte * 0.004) 

## select variables of interests

preds<-c("pf_pa")
#age<- c("lower_age","upper_age")
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
bio<- c("bio3","bio15","bio16","bio17")##,"bio16","bio17")##,"bio11"
ro<- c("Roq1","Roq2")##,"Roq3","Roq4")
other<- c("year_start", "elev","NDVI")##,"rural_urban")
soc<- c("hindex","lgdppc","lpop_den")##"Pop"

allPfalcAF<- c(preds, qtempmean, qprecip, bio, ro, other, soc)


datapfAF<- dat1[,allPfalcAF]
head(datapfAF,3)
dim(datapfAF)

##Build CART trees

## First fit forcing the tree to have lots of branches so we can
## examine it and figure out where to trim

CP=0.0005
MS=20
cnt<-rpart.control(minsplit=MS, cp=CP, xval=100)

f.null<-rpart(pf_pa ~  . , data=datapfAF, method = "class", control=cnt)
plotcp(f.null) ## use this to decide on a cp/size for trimming the tree
printcp(f.null)
graphics.off()

x11()
par(mfrow=c(1,1))
plot(f.null, uniform=TRUE, margin=0.0075)
text(f.null, digits = 1, use.n=TRUE, cex=0.5)


f.null$cptable


mm<-min(signif(f.null$cptable[,4],3))+mean(f.null$cptable[,5])
mm

w<-which(f.null$cptable[,4]<mm)
w

CP.new<-f.null$cptable[min(w)-1,1]
CP.new


## looking for 1 standard error from 

## now trimming based on the above
cnt<-rpart.control(minsplit=200, cp=CP.new, xval=100)  ##CP.new = 0.004568528
f<-rpart(pf_pa ~  . , data=datapfAF, method = "class", control=cnt)
plotcp(f)
graphics.off()

## Plot tree
x11()
par(mfrow=c(1,1))
plot(f, uniform=TRUE, margin=0.0075); text(f, digits=1, use.n=TRUE, cex=0.5)

rpart.plot(f, type = 4, extra = 104, box.palette = "GnBu", branch.lty=3,shadow.col="gray",nn=TRUE)


## squared variables of interests

names(datapfAF)

datapfAF$bio03sq<-datapfAF$bio3^2
datapfAF$qtmean1sq<- datapfAF$qtmean1^2
datapfAF$qprecp1sq<- datapfAF$qprecp1^2



## GLM models

m.null<-glm(pf_pa ~ 1, family="binomial", data=datapfAF)
m.null
m.full<-glm(pf_pa ~ . , family="binomial", data=datapfAF)
m.full

modelpfAF<-step(m.null, scope=formula(m.full), direction="forward", criterion = "BIC")
summary(modelpfAF)


## Plot the residuals

plot(fitted(modelpfAF), res)
abline(0,0)
qqnorm(res, main = expression('Q-Q plot for '*italic(P.~falciparum)*' in Africa'), bty="o")
qqline(res, col = "steelblue", lwd = 3)
legend(-4.8, 4, legend = "A", bty = "n", cex = 2.5)

## plot the density curve

plot(density(res, adjust = 1.2), main = expression('Density plot for '*italic(P.~falciparum)*' in Africa'), lwd = 2)
legend(-5.4, 0.43, legend = "B", bty = "n", cex = 2.5)


## Build prediction curves using outputs from GLM and CART models

## lets look at marginal predictions, based on particular predictors
## for bio3: Isothermality

ey <- expression(italic(P. ~ falciparum) ~ predicted ~ presence)

o<-order(datapfAF$bio3)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd = TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$bio3[o], datapfAF$pf_pa, xlab="Temperature oscillation (%)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(30,100),main="GLM fit bio3-Africa")
fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio3=datapfAF$bio3[o])
y.m2<-predict(loess(fit~bio3, data=fit.dat.m, span = 1), data.frame(bio3=seq(30, 100, by=1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio03[o]/10, y.sym, col=2, lwd=3)
lines(seq(30, 100, by=1), y.m2$fit, col=2, lwd=3)
legend(45, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred

f.pred <- ifelse(f.pred=="1",0,1)
f.pred

plot(jitter(datapfAF$bio3[o]), datapfAF$pf_pa, ylim=c(0,1),xlim=c(30,100),
     xlab="Temperature oscillation (%)", ylab= ey, cex.lab = 1.1,
     main="CART fit bio3-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio3=datapfAF$bio3[o])
fit.dat.f
y.f2<-predict(loess(fit~bio3, data=fit.dat.f, span = 1), data.frame(bio3=seq(30, 100, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(30, 100, by=1), y.f2$fit, col=4, lwd=3)
legend(45, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=105, y=1.3, legend="E", xpd = NA, bty="n", cex = 1.7)



## for bio16: Precipitation of wettest quarter (mm)
## lets look at marginal predictions, based on particular predictors

o<-order(datapfAF$bio16)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(datapfAF$bio16[o], datapfAF$pf_pa, xlab="Precipitation wettest qtr (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,2000),main="GLM fit bio16-Africa")
fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio16=datapfAF$bio16[o])
y.m2<-predict(loess(fit~bio16, data=fit.dat.m, span = 1), data.frame(bio16=seq(0, 2000, by=10)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,2000, by=10), y.m2$fit, col=2, lwd=3)
legend(200, 0.3, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)
plot(jitter(datapfAF$bio16[o]), datapfAF$pf_pa, ylim=c(0,1),xlim=c(0,2000),
     xlab="Precipitation wettest qtr (mm)", ylab= ey, cex.lab = 1.1,
     main="CART fit bio16-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio16=datapfAF$bio16[o])
fit.dat.f
y.f2<-predict(loess(fit~bio16, data=fit.dat.f, span=1), data.frame(bio16=seq(0, 2000, by=10)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 2000, by=10), y.f2$fit, col=4, lwd=3)
legend(200, 0.3, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=2150, y=1.3, legend="F", xpd = NA, bty="n", cex = 1.7)



##bio17
## Precipitation of driest quarter

o<-order(datapfAF$bio17)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$bio17[o], datapfAF$pf_pa, xlab="Precipitation driest qtr (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,400),main="GLM fit bio17-Africa")
fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, bio17=datapfAF$bio17[o])
y.m2<-predict(loess(fit~bio17, data=fit.dat.m, span=1), data.frame(bio17=seq(0, 400, by=1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,400, by=1), y.m2$fit, col=2, lwd=3)
legend(70, 0.3, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$bio17[o]), datapfAF$pf_pa,ylim=c(0,1),xlim=c(0,400),
     xlab="Precipitation driest qtr(mm)", ylab= ey, cex.lab = 1.1, 
     main="CART fit bio17-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, bio17=datapfAF$bio17[o])
fit.dat.f
y.f2<-predict(loess(fit~bio17, data=fit.dat.f, span=1), data.frame(bio17=seq(0, 400, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 400, by=1), y.f2$fit, col=4, lwd=3)
legend(70, 0.3, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=415, y=1.3, legend="G", xpd = NA, bty="n", cex = 1.7)



## SOCIO-DEMOGRAPHIC VARIABLES #####

###FOR YEAR START 

o<-order(datapfAF$year_start)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$year_start[o], datapfAF$pf_pa[o], xlab="Year of the study",
     ylab= ey, cex.lab=1.1, xlim=c(1990,2020),ylim=c(0,1),
     main="GLM fit year_study-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, year_start=datapfAF$year_start[o])
y.m2<-predict(loess(fit~year_start, data=fit.dat.m, span=0.85), data.frame(year_start=seq(1990,2020, by=1)), se=TRUE)
y.m2

lines(seq(1990, 2020, by=1), y.m2$fit, col=2, lwd=3)
legend(1995, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$year_start[o]), datapfAF$pf_pa[o],
     xlab="Year of the study", ylab= ey, cex.lab=1.1,
     xlim=c(1990,2020), ylim=c(0,1),main="CART fit year_study-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, year_start=datapfAF$year_start[o])
fit.dat.f
y.f2<-predict(loess(fit~year_start, data=fit.dat.f, span=1), data.frame(year_start=seq(1990, 2020, by=1)), se=TRUE)
y.f2$fit[y.f2$fit < 0] <- 0
##lines
lines(seq(1990, 2020, by=1), y.f2$fit, col=4, lwd=3)
legend(1995, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=2023, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)



###FOR LGDP ADJUSTED 

o<-order(datapfAF$lgdppc)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$lgdppc[o], datapfAF$pf_pa[o], xlab="log(GDP) per capita",
     ylab= ey, cex.lab=1.1, xlim=c(5,11),ylim=c(0,1),
     main="GLM fit GDPPC-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, sGDP=datapfAF$lgdppc[o])
y.m2<-predict(loess(fit~sGDP, data=fit.dat.m, span=1), data.frame(sGDP=seq(5,11, by=0.1)), se=TRUE)
y.m2

##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(5, 11, by=0.1), y.m2$fit, col=2, lwd=3)
legend(6, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$lgdppc[o]), datapfAF$pf_pa[o],
     xlab="log(GDP) per capita", ylab= ey, cex.lab=1.1, 
     xlim=c(5,11), ylim=c(0,1),main="CART fit GDP-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, GDP=datapfAF$lgdppc[o])
fit.dat.f
y.f2<-predict(loess(fit~GDP, data=fit.dat.f, span=1), data.frame(GDP=seq(5,11, by=0.1)), se=TRUE)
y.f2
##lines
lines(seq(5, 11, by=0.1), y.f2$fit, col=4, lwd=3)
legend(6, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=11.5, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)



###FOR population density####

o<-order(datapfAF$lpop_den)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

x <- datapfAF$lpop_den[o] - 4
x[x < 0] <- 0
x

## plot based off predictions from the glm for data pre-2004
plot(x,datapfAF$pf_pa[o], xlab=" log population density (Km2)",
     ylab= ey, cex.lab= 1.1, xlim=c(0,12),ylim=c(0,1),
     main="GLM fit pop density-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, pop_den= (datapfAF$lpop_den[o] - 4))
y.m2<-predict(loess(fit~pop_den, data=fit.dat.m, span=1), data.frame(pop_den=seq(0, 12, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit < 0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 12, by=0.1), y.m2$fit, col=2, lwd=3)
legend(2, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
x1 <- datapfAF$lpop_den[o] - 4
x1[x1 < 0] <- 0

f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(x1), datapfAF$pf_pa[o],
     xlab="log population density (Km2)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,12), ylim=c(0,1),main="CART fit pop density-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, pop_den=(datapfAF$lpop_den[o] - 4))
fit.dat.f
y.f2<-predict(loess(fit~pop_den, data=fit.dat.f, span=1), data.frame(pop_den=seq(0, 12, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit > 1] <- 1
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 12, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(2, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=13, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)


###FOR human density####

o<-order(datapfAF$hindex)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$hindex[o],datapfAF$pf_pa[o], xlab="Human development index",
     ylab= ey, cex.lab= 1.1, xlim=c(0,0.8),ylim=c(0,1),
     main="GLM fit HDI-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, hindex=datapfAF$hindex[o])
y.m2<-predict(loess(fit~hindex, data=fit.dat.m, span=1), data.frame(hindex=seq(0, 0.8, by=0.05)), se=TRUE)
y.m2$fit[y.m2$fit >1] <- 1

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 0.8, by=0.05), y.m2$fit, col=2, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$hindex[o]), datapfAF$pf_pa[o],
     xlab="Human development index", ylab= ey, cex.lab=1.1, 
     xlim=c(0,0.8), ylim=c(0,1),main="CART fit HDI-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, hindex=datapfAF$hindex[o])
fit.dat.f
y.f2<-predict(loess(fit~hindex, data=fit.dat.f, span=1), data.frame(hindex=seq(0,0.8, by=0.05)), se=TRUE)
y.f2$fit[y.f2$fit > 1] <- 1
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 0.8, by=0.05), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=0.84, y=1.3, legend="D", xpd = NA, bty="n", cex = 1.7)



##FOR qtmean1 ##########

o<-order(datapfAF$qtmean1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$qtmean1[o], datapfAF$pf_pa, xlab="Temperature (°C)",
     ylab= ey, cex.lab = 1.1,
     xlim=c(5,40), ylim=c(0,1),main="GLM fit qtmean1-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qtmean1=datapfAF$qtmean1[o])
y.m2<-predict(loess(fit~qtmean1, data=fit.dat.m, span=1), data.frame(qtmean1=seq(5, 40, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit <0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(5, 40, by=0.1), y.m2$fit, col=2, lwd=3)
legend(15, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)
f.pred
plot(jitter(datapfAF$qtmean1[o]), datapfAF$pf_pa,
     xlab="Temperature (°C)", ylab= ey, cex.lab = 1.1, 
     xlim=c(5,40), ylim=c(0,1),main="CART fit qtmean1-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qtmean1=datapfAF$qtmean1[o])
fit.dat.f
y.f2<-predict(loess(fit~qtmean1, data=fit.dat.f, span=1), data.frame(qtmean1=seq(5, 40, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit <0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(5, 40, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(15, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=44, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)


##FOR qtmean2

o<-order(datapfAF$qtmean2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$qtmean2[o],datapfAF$pf_pa, xlab="Temperature (°C)",
     ylab= ey, cex.lab = 1.1, 
     xlim=c(10,40), ylim=c(0,1),main="GLM fit qtmean2-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qtmean2=datapfAF$qtmean2[o])
y.m2<-predict(loess(fit~qtmean2, data=fit.dat.m, span=1), data.frame(qtmean2=seq(10, 40, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit < 0] <- 0

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(10, 40, by=0.1), y.m2$fit, col=2, lwd=3)
legend(16, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)
f.pred
plot(jitter(datapfAF$qtmean2[o]), datapfAF$pf_pa,
     xlab="Temperature (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,40), ylim=c(0,1),main="CART fit qtmean2-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qtmean2=datapfAF$qtmean2[o])
fit.dat.f
y.f2<-predict(loess(fit~qtmean2, data=fit.dat.f, span=1), data.frame(qtmean2=seq(10, 40, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit < 0] <- 0


##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(10, 40, by=0.1), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(16, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x= 43, y=1.3, legend="A", xpd = NA, bty="n", cex = 1.7)


###FOR qprecp1

o<-order(datapfAF$qprecp1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$qprecp1[o], datapfAF$pf_pa[o], xlab="Precipitation (mm)",
     ylab= ey, cex.lab=1.1,  xlim=c(0,500),ylim=c(0,1),
     main="GLM fit qprecp1-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qprecp1=datapfAF$qprecp1[o])
y.m2<-predict(loess(fit~qprecp1, data=fit.dat.m, span=1), data.frame(qprecp1=seq(0, 500, by=1)), se=TRUE)
y.m2

lines(seq(0, 500, by=1), y.m2$fit, col=2, lwd=3)
legend(100, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$qprecp1[o]), datapfAF$pf_pa[o],
     xlab="Precipitation (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,500), ylim=c(0,1),main="CART fit qprecp1-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qprecp1=datapfAF$qprecp1[o])
fit.dat.f
y.f2<-predict(loess(fit~qprecp1, data=fit.dat.f, span=1), data.frame(qprecp1=seq(0, 500, by=1)), se=TRUE)
y.f2
##lines
lines(seq(0, 500, by=1), y.f2$fit, col=4, lwd=3)
legend(100, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=540, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)

###FOR qprecp2

o<-order(datapfAF$qprecp2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$qprecp2[o], datapfAF$pf_pa[o], xlab="Precipitation (mm)",
     ylab= ey, cex.lab=1.1, xlim=c(0,500),ylim=c(0,1),
     main="GLM fit qprecp2-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, qprecp2=datapfAF$qprecp2[o])
y.m2<-predict(loess(fit~qprecp2, data=fit.dat.m, span=1), data.frame(qprecp2=seq(0, 500, by=1)), se=TRUE)
y.m2$fit[y.m2$fit >1] <- 1

lines(seq(0, 500, by=1), y.m2$fit, col=2, lwd=3)
legend(50, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$qprecp2[o]), datapfAF$pf_pa[o],
     xlab="Precipitation (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,500), ylim=c(0,1),main="CART fit qprecp2-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, qprecp2=datapfAF$qprecp2[o])
fit.dat.f
y.f2<-predict(loess(fit~qprecp2, data=fit.dat.f,span=1), data.frame(qprecp2=seq(0, 500, by=1)), se=TRUE)
y.f2
##lines
lines(seq(0, 500, by=1), y.f2$fit, col=4, lwd=3)
legend(50, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=540, y=1.3, legend="B", xpd = NA, bty="n", cex = 1.7)


### Elevation (DEM)

o<-order(datapfAF$elev)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$elev[o], datapfAF$pf_pa[o], xlab="Elevation (meters)",
     ylab= ey, cex.lab=1.1, xlim=c(0,3500),ylim=c(0,1),
     main="GLM fit elevation-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, elev=datapfAF$elev[o])
y.m2<-predict(loess(fit~elev, data=fit.dat.m, span=1), data.frame(elev=seq(0,3500, by=10)), se=TRUE)
y.m2

lines(seq(0, 3500, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$elev[o]), datapfAF$pf_pa[o],
     xlab="Elevation (meters)", ylab= ey, cex.lab=1.1,
     xlim=c(0,3500), ylim=c(0,1),main="CART fit Elevation-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, elev=datapfAF$elev[o])
fit.dat.f
y.f2<-predict(loess(fit~elev, data=fit.dat.f, span=1), data.frame(elev=seq(0,3500, by=10)), se=TRUE)
y.f2$fit[y.f2$fit <0] <- 0
##lines
lines(seq(0, 3500, by=10), y.f2$fit, col=4, lwd=3)
legend(100, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=3700, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)

###NDVI

o<-order(datapfAF$NDVI)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$NDVI[o], datapfAF$pf_pa[o], xlab="NDVI index",
     ylab= ey, cex.lab=1.1, xlim=c(0,1),ylim=c(0,1),
     main="GLM fit NDVI-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, NDVI=datapfAF$NDVI[o])
y.m2<-predict(loess(fit~NDVI, data=fit.dat.m, span=0.85), data.frame(NDVI=seq(0,1, by=0.05)), se=TRUE)
y.m2

lines(seq(0, 1, by=0.05), y.m2$fit, col=2, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 

## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$NDVI[o]), datapfAF$pf_pa[o],
     xlab="NDVI index", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit NDVI-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, NDVI=datapfAF$NDVI[o])
fit.dat.f
y.f2<-predict(loess(fit~NDVI, data=fit.dat.f, span=0.85), data.frame(NDVI=seq(0,1, by=0.05)), se=TRUE)
y.f2
##lines
lines(seq(0,1, by=0.05), y.f2$fit, col=4, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1.05, y=1.3, legend="D", xpd = NA, bty="n", cex = 1.7)


###FOR R0q1

o<-order(datapfAF$Roq1)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$Roq1[o], datapfAF$pf_pa[o], xlab="R0-Temp 1st quarter",
     ylab= ey, cex.lab=1.1,  xlim=c(0,1),ylim=c(0,1),
     main="GLM fit R0q1-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, Roq1=datapfAF$Roq1[o])
y.m2<-predict(loess(fit~Roq1, data=fit.dat.m, span=1), data.frame(Roq1=seq(0,1, by=0.1)), se=TRUE)
y.m2

lines(seq(0, 1, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$Roq1[o]), datapfAF$pf_pa[o],
     xlab="R0-Temp 1st quarter", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit R0q1-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, Roq1=datapfAF$Roq1[o])
fit.dat.f
y.f2<-predict(loess(fit~Roq1, data=fit.dat.f, span=1), data.frame(Roq1=seq(0,1, by=0.1)), se=TRUE)
y.f2
##lines
lines(seq(0, 1, by=0.1), y.f2$fit, col=4, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

legend(x=1.05, y=1.3, legend="E", xpd = NA, bty="n", cex = 1.7)


##### FOR ROq2

o<-order(datapfAF$Roq2)

op <- par(
  mar=c(4,4.5,2,4),
  oma=c(2,2,2,2),
  xpd = TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(datapfAF$Roq2[o], datapfAF$pf_pa[o], xlab="R0-Temp 2nd quarter",
     ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="GLM fit R0q2-Africa")

fitted<-as.numeric(modelpfAF$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, Roq2=datapfAF$Roq2[o])
y.m2<-predict(loess(fit~Roq2, data=fit.dat.m, span=1), data.frame(Roq2=seq(0,1, by=0.1)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1, by=0.1), y.m2$fit, col=2, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred <- ifelse(f.pred=="1",0,1)

plot(jitter(datapfAF$Roq2[o]), datapfAF$pf_pa[o],
     xlab="R0-Temp 2nd quarter", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1), ylim=c(0,1),main="CART fit R0q2-Africa")
fittedf<- as.numeric(f.pred[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, Roq2=datapfAF$Roq2[o])
fit.dat.f
y.f2<-predict(loess(fit~Roq2, data=fit.dat.f, span=1), data.frame(Roq2=seq(0, 1, by=0.1)), se=TRUE)
y.f2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1, by=0.1), y.f2$fit, col=4, lwd=3)
legend(0.1, 0.28, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=1.05, y=1.3, legend="C", xpd = NA, bty="n", cex = 1.7)

