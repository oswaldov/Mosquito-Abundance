##author: Oswaldo Villena
##Load required packages

library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)
library(quantmod)
library(stargazer)
library(tidyverse)
library(ggplot2)
library(GGally)

##Load malaria-climate data sets
mosqcountcx <-read.csv("data/TrainDataCx.csv", header=TRUE)
head(mosqcountcx,3)

test_Data <- read.csv("data/TestDataCx.csv", header=TRUE)
head(test_Data,3)

##Be sure year and month are factors

mosqcountcx$year <- as.factor(mosqcountcx$year)
mosqcountcx$month <- as.factor(mosqcountcx$month)

test_Data$year <- as.factor(test_Data$year)
test_Data$month <- as.factor(test_Data$month)


##Preliminary visualizations

par(mfrow=c(2,2))
plot(mosqcountcx$site, mosqcountcx$Cxpa)
plot(mosqcountcx$season, mosqcountcx$Cxpa)
plot(mosqcountcx$rainfall_season, mosqcountcx$Cxpa)
plot(mosqcountcx$elevation, mosqcountcx$Cxpa)

par(mfrow=c(2,2))
plot(mosqcountcx$precip, mosqcountcx$Cxpa)
plot(mosqcountcx$tmean, mosqcountcx$Cxpa)
plot(mosqcountcx$tmax, mosqcountcx$Cxpa)
plot(mosqcountcx$tmin, mosqcountcx$Cxpa)

##Check for autocorrelation
##geographic
locat <- c(1,12,13,14)
summary(mosqcount[,locat])
pairs(mosqcount[,locat])

##precip variables
precip <- c(8,15,16,17)
summary(mosqcount[,precip])
pairs(mosqcount[,precip])

## temperature variables 
temp<- c(18,19,20)
summary(mosqcount[,temp])
pairs(mosqcount[,temp])


## Selection of variables of interest

preds<-c("Cxpa")
geo <- c("trap_type","month") ##"rainfall_season") ##"season",
loc <- c("elevation","distance")
precip <- c("precip","preciplag1","preciplag2") ##,"preciplag3"
temp <- c("tmean","tmax","tmin")
tempmean <- c("tmean","tmeanlag1","tmeanlag2") ##,"tmeanlag3"
tempmax <- c("tmax","tmaxlag1","tmaxlag2","tmaxlag3")
tempmin <- c("tmin","tminlag1","tminlag2","tminlag3")

all<- c(preds, geo, loc, precip, tempmean)


dataNS <- mosqcountcx[,all]
head(dataNS)
dim(dataNS)

##Test data
datatestNS <- test_Data[,c(12,2,6,13,14,15,16,17,18,19,20)]
head(datatestNS)


##CART model
## First fit forcing the tree to have lots of branches so we can
## examine it and figure out where to trim
CP=0.0005
MS=100
cnt<-rpart.control(minsplit=MS, cp=CP, xval=100)

f.null<-rpart(Cxpa ~  . , data=dataNS, method="class", control=cnt)
plotcp(f.null) ## use this to decide on a cp/size for trimming the tree
printcp(f.null)
graphics.off()


x11()
par(mfrow=c(1,1))
#plot(f.null, uniform=TRUE, margin=0.0075)
#text(f.null, digits=1, use.n=TRUE, cex=0.5)

rpart.plot(f.null, type = 4, extra = 104, branch.lty = 3)
printcp(f.null)
plotcp(f.null)


## Find the lowest cp value, so the cross-validated error rate is minimum

CP.new <- f.null$cptable[which.min(f.null$cptable[,"xerror"]),"CP"]
CP.new


## now trimming based on the above
cnt<-rpart.control(minsplit=MS, cp=0.0032, xval=100)
f<-rpart(Cxpa ~  . , data=dataNS, method="class", control=cnt)
plotcp(f)
graphics.off()


x11()
par(mfrow=c(1,1))

rpart.plot(f, type = 4, extra = 104, branch.lty = 3, gap = 0.6, space = 0.8, tweak = 1.1)


## Validation

library(ggplot2)
library(gmodels)
library(e1071)
library(gridExtra)
library(randomForest)


preds.rpart = predict(f,newdata = datatestNS,type = "class")
CrossTable(datatestNS$Cxpa,preds.rpart,chisq = F,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)

((2854 + 618)/nrow(datatestNS))*100

PrecisionS <- 618/886 ; PrecisionS

RecallS <- 618/1145 ; RecallS

F1S <- 2 * ((70*54)/(70+54)) ; F1S


## here's what rpart gives for its trimmed tree if I don't include a
## controler

head(dataCX,2)

par(mfrow=c(1,1))
f.d<-rpart(Cxpa ~  . , method="class", data=dataCX)
printcp(f.d)
plotcp(f.d)

par(mfrow=c(1,1))
rpart.plot(f.d,type=4, extra=104, box.palette = "GnBu",branch.lty=3,shadow.col="gray")



##GLM model
##models


m.null<- glm(Cxpa ~ 1, family="binomial", data=dataCX) ## null model
m.null
m.full<- glm(Cxpa ~ .+ tmean * precip, family="binomial", data=dataCX)
m.full


modelcx<- step(m.null, scope=formula(m.full), direction="both", criterion = "BIC")
summary(modelcx)

##Best model

bmodel <- glm(Cxpa ~ site + trap_type + month + tmean + tmeanlag2 + preciplag1 + precip + preciplag2 + 
                tmeanlag1, family = "binomial", data = dataCX)

summary(bmodel)

predict <- predict(bmodel, datatest, type = 'response')

##validation
# confusion matrix
table_mat <- table(datatest$Cxpa, predict > 0.5)
table_mat


accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test

precGLM <- 543/820 ; precGLM
recallGLM <- 543/1145 ; recallGLM

F1GLM <- 2 * ((66*47)/(66+47)) ; F1GLM

## Q-Q plots

qr2<- qresiduals(modelcx)
length(qr2)
summary(qr2)


par(mfrow=c(1,1))
qqnorm(qr2, ylim = c(-6,6), xlim=c(-5,5), main = "Normal Q-Q Plot", las=1, bty="o"); qqline(qr2)


## Plot for the diferent variables using fitted values from the best model

##TEMPERATURE MEAN
##Using GLM

o<-order(dataCX$tmean)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmean[o], dataCX$Cxpa[o], xlab="Temperature (°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmean=dataCX$tmean[o])
y.m2<-predict(loess(fit~tmean, data=fit.dat.m, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmean[o]), dataCX$Cxpa[o],
     xlab="Temperature (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmean=dataCX$tmean[o])
fit.dat.f
y.f2<-predict(loess(fit~tmean, data=fit.dat.f, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=27, y=1.2, legend="F", xpd = NA, bty="n", cex = 1.7)

##Precipitation

o<-order(dataCX$precip)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(dataCX$precip[o], dataCX$Cxpa, xlab="Precipitation (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,1200),main="GLM fit - Precipitation")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, precip=dataCX$precip[o])
y.m2<-predict(loess(fit~precip, data=fit.dat.m, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.m2
##lines
lines(seq(0,1200, by=1), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$precip[o]), dataCX$Cxpa,ylim=c(0,1),xlim=c(0,1200),
     xlab="Precipitation (mm)", ylab = ey, cex.lab = 1.1,
     main="CART fit - Precipitation")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, precip=dataCX$precip[o])
fit.dat.f
y.f2<-predict(loess(fit~precip, data=fit.dat.f, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2


##lines
lines(seq(0, 1200, by=1), y.f2$fit, col=4, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1180, y=1.2, legend="C", xpd = NA, bty="n", cex = 1.7)



###Preciplag1

o<-order(dataCX$preciplag1)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$preciplag1[o],dataCX$Cxpa[o], xlab="Precipitation lag1 (mm)",
     ylab= ey, cex.lab= 1.1, xlim=c(0,1200),ylim=c(0,1),
     main="GLM fit - Precipitation lag1")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag1=dataCX$preciplag1[o])
y.m2<-predict(loess(fit~preciplag1, data=fit.dat.m, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$preciplag1[o]), dataCX$Cxpa[o],
     xlab="Precipitation lag1 (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag1")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag1=dataCX$preciplag1[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag1, data=fit.dat.f, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1200, y=1.2, legend="D", xpd = NA, bty="n", cex = 1.6)



## Precipitation lag2

o<-order(dataCX$preciplag2)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$preciplag2[o],dataCX$Cxpa[o], xlab="Precipitation lag2 (mm)",
     ylab= ey, cex.lab= 1.1, xlim=c(0,1200),ylim=c(0,1),
     main="GLM fit - Precipitation lag2")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag2=dataCX$preciplag2[o])
y.m2<-predict(loess(fit~preciplag2, data=fit.dat.m, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$preciplag2[o]), dataCX$Cxpa[o],
     xlab="Precipitation lag2 (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag2")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag2=dataCX$preciplag2[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag2, data=fit.dat.f, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

legend(x=1200, y=1.2, legend="E", xpd = NA, bty="n", cex = 1.6)


##Distance to anthropogenic features

o<-order(dataCX$distance)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$distance[o], dataCX$Cxpa, xlab="Distance - anthropogenic features (m)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,2500),main="GLM fit - Distance")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, distance=dataCX$distance[o])
y.m2<-predict(loess(fit~distance, data=fit.dat.m, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines
lines(seq(0,2500, by=100), y.m2$fit, col=2, lwd=3)
legend(200, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$distance[o]), dataCX$Cxpa,ylim=c(0,1),xlim=c(0,2500),
     xlab="Distance - anthropogenic features (m)", ylab= ey, cex.lab = 1.1, 
     main="CART fit - Distance")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, distance=dataCX$distance[o])
fit.dat.f
y.f2<-predict(loess(fit~distance, data=fit.dat.f, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit
##lines
lines(seq(0, 2500, by=100), y.f2$fit, col=4, lwd=3)
legend(200, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=2300, y=1.2, legend="B", xpd = NA, bty="n", cex = 1.7)






