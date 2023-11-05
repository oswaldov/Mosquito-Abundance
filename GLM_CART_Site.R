##author: Oswaldo Villena
## Load required packages

library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)
library(quantmod)
library(stargazer)
library(tidyverse)


## Load malaria-climate data sets
mosqcountcx <-read.csv("data/TrainDataCx.csv", header=TRUE)
test_Data <- read.csv("data/TestDataCx.csv", header=TRUE)

## Be sure year and month are factors
mosqcountcx$year <- as.factor(mosqcountcx$year)
mosqcountcx$month <- as.factor(mosqcountcx$month)

test_Data$year <- as.factor(test_Data$year)
test_Data$month <- as.factor(test_Data$month)

## Preliminary visualizations
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


##Check for auto-correlation
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
geo <- c("site","trap_type","month") 
loc <- c("elevation","distance")
precip <- c("precip","preciplag1","preciplag2") 
temp <- c("tmean","tmax","tmin")
tempmean <- c("tmean","tmeanlag1","tmeanlag2") 

all<- c(preds, geo, loc, precip, tempmean)

dataCX <- mosqcountcx[,all]


## Train data
dataTrain <- mosqcountcx[,c(1,2,6,12,14,15,16,17,18)]

## Test data
dataTest <- test_Data[,c(1,2,6,12,14,15,16,17,18)]

## CART model 
## First fit forcing the tree to have lots of branches so we can
## examine it and figure out where to trim

CP=0.0005
MS=150
cnt<-rpart.control(minsplit=MS, cp=CP, xval=100)

f.null<-rpart(Cxpa ~  . , data=dataTrain, method="class", control=cnt)
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
cnt<-rpart.control(minsplit=MS, cp=CP.new, xval=100)
f<-rpart(Cxpa ~  . , data=dataTrain, method="class", control=cnt)
plotcp(f)
graphics.off()

x11()
par(mfrow=c(1,1))

rpart.plot(f, type = 4, extra = 104, branch.lty = 3, gap = 0.6, space = 0.8, tweak = 1.2)


## Validation using data test

library(ggplot2)
library(gmodels)
library(e1071)
library(gridExtra)
library(randomForest)

preds.rpart = predict(f,newdata = dataTest,type = "class")
CrossTable(dataTest$Cxpa,preds.rpart,chisq = F,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)

## here's what rpart gives for its trimmed tree if I don't include a
## controler
par(mfrow=c(1,1))
f.d<-rpart(Cxpa ~  . , method="class", data=dataTrain)
printcp(f.d)
plotcp(f.d)

par(mfrow=c(1,1))
rpart.plot(f.d,type=4, extra=104, box.palette = "GnBu",branch.lty=3,shadow.col="gray")


## GLM model
##models
m.null<- glm(Cxpa ~ 1, family="binomial", data=dataTrain) ## only intercept model
m.null
m.full<- glm(Cxpa ~ . + tmean * precip, family="binomial", data=dataTrain)
m.full

modelcxS<- step(m.null, scope=formula(m.full), direction="both", criterion = "BIC")
summary(modelcxS)

bestmodel <- glm(Cxpa ~ site + trap_type + month + tmean + preciplag2 + precip + preciplag1 + 
                   tmean*precip, family = "binomial", data = dataTrain)

summary(bestmodel)

##validation using data test
predict <- predict(bestmodel, dataTest, type = 'response')
# confusion matrix
table_mat <- table(dataTest$Cxpa, predict > 0.5)
table_mat


accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test

## Q-Q plots for residuals test
qr2<- qresiduals(modelcxS)
length(qr2)
summary(qr2)

par(mfrow=c(1,1))
qqnorm(qr2, ylim = c(-6,6), xlim=c(-5,5), main = "Normal Q-Q Plot", las=1, bty="o"); qqline(qr2)

legend(3, 6.8, "B", cex = 2.5, box.lty = 0, bg = "transparent") 


## Plot for the diferent variables using fitted values from the best model

## TEMPERATURE MEAN

##Using GLM
ey <- expression(italic(Culex) ~ predicted ~ occurrence)

o<-order(dataTrain$tmean)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataTrain$tmean[o], dataTrain$Cxpa[o], xlab="Temperature (°C)",
     ylab= ey, cex.lab=1.4, font=2, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature")

fitted<-as.numeric(modelcxS$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmean=dataTrain$tmean[o])
y.m2<-predict(loess(fit~tmean, data=fit.dat.m, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)


## lets try it with the tree
f.pred<-predict(f, type="vector")
table(f.pred)
f.pred1 <- recode(f.pred, '1' = 0, '2' = 1)
table(f.pred1)


plot(jitter(dataTrain$tmean[o]), dataTrain$Cxpa[o],
     xlab="Temperature (°C)", ylab= ey, cex.lab=1.4, font=2, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmean=dataTrain$tmean[o])
fit.dat.f
y.f2<-predict(loess(fit~tmean, data=fit.dat.f, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)

##legend for the whole graph
legend(x=25, y=1.2, legend="D", xpd = NA, bty="n", cex = 2.2)


##PRECIPITATION

##using GLM
o<-order(dataTrain$precip)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(dataTrain$precip[o], dataTrain$Cxpa, xlab="Precipitation (mm)", font = 2,
     ylab= ey, cex.lab = 1.4, ylim=c(0,1), xlim=c(0,1200),main="GLM fit - Precipitation")
fitted<-as.numeric(modelcxS$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, precip=dataTrain$precip[o])
y.m2<-predict(loess(fit~precip, data=fit.dat.m, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.m2
##lines
lines(seq(0,1200, by=1), y.m2$fit, col=2, lwd=3)


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataTrain$precip[o]), dataTrain$Cxpa,ylim=c(0,1),xlim=c(0,1200),
     xlab="Precipitation (mm)", ylab = ey, cex.lab = 1.4, font=2,
     main="CART fit - Precipitation")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, precip=dataTrain$precip[o])
fit.dat.f
y.f2<-predict(loess(fit~precip, data=fit.dat.f, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2


##lines
lines(seq(0, 1200, by=1), y.f2$fit, col=4, lwd=3)

##legend for the whole graph
legend(x=1050, y=1.2, legend="E", xpd = NA, bty="n", cex = 2.2)


##Precipitation_lag1

##using GLM
o<-order(dataTrain$preciplag1)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataTrain$preciplag1[o],dataTrain$Cxpa[o], xlab="Precipitation lag1 (mm)",
     ylab= ey, cex.lab= 1.4,font = 2, xlim=c(0,1200),ylim=c(0,1), 
     main="GLM fit - Precipitation lag1")

fitted<-as.numeric(modelcxS$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag1=dataTrain$preciplag1[o])
y.m2<-predict(loess(fit~preciplag1, data=fit.dat.m, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataTrain$preciplag1[o]), dataTrain$Cxpa[o],
     xlab="Precipitation lag1 (mm)", ylab= ey, cex.lab=1.4, font=2, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag1")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag1=dataTrain$preciplag1[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag1, data=fit.dat.f, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)

##legend for the whole graph
legend(x=1050, y=1.2, legend="C", xpd = NA, bty="n", cex = 2.2)


## Precipitation_lag2

##using GLM
o<-order(dataTrain$preciplag2)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataTrain$preciplag2[o],dataTrain$Cxpa[o], xlab="Precipitation lag2 (mm)",
     ylab= ey, cex.lab= 1.4, font=2, xlim=c(0,1200),ylim=c(0,1),
     main="GLM fit - Precipitation lag2")

fitted<-as.numeric(modelcxS$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag2=dataTrain$preciplag2[o])
y.m2<-predict(loess(fit~preciplag2, data=fit.dat.m, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataTrain$preciplag2[o]), dataTrain$Cxpa[o],
     xlab="Precipitation lag2 (mm)", ylab= ey, cex.lab=1.4, font=2, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag2")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag2=dataTrain$preciplag2[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag2, data=fit.dat.f, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)

##legend for the whole graph
legend(x=1050, y=1.2, legend="D", xpd = NA, bty="n", cex = 2.2)

##Distance to anthropogenic features

##Using GLM
o<-order(dataTrain$distance)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataTrain$distance[o], dataTrain$Cxpa, xlab="Distance (m)",
     ylab= ey, cex.lab = 1.4, font=2, ylim=c(0,1), xlim=c(0,2500),main="GLM fit - Distance")
fitted<-as.numeric(modelcxS$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, distance=dataTrain$distance[o])
y.m2<-predict(loess(fit~distance, data=fit.dat.m, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,2500, by=100), y.m2$fit, col=2, lwd=3)

## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataTrain$distance[o]), dataTrain$Cxpa,ylim=c(0,1),xlim=c(0,2500),
     xlab="Distance (m)", ylab= ey, cex.lab = 1.4, font=2, 
     main="CART fit - Distance")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, distance=dataTrain$distance[o])
fit.dat.f
y.f2<-predict(loess(fit~distance, data=fit.dat.f, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit
##lines
lines(seq(0, 2500, by=100), y.f2$fit, col=4, lwd=3)

##legend for the whole graph
legend(x=2150, y=1.2, legend="F", xpd = NA, bty="n", cex = 2.2)


















