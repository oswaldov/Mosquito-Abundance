library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)
library(quantmod)
library(stargazer)
library(tidyverse)

getwd()
## cleaning up the malaria-climate data file
mosqcountcx <-read.csv("TrainDataCx.csv", header=TRUE)
head(mosqcountcx,3)
dim(mosqcountcx)
names(mosqcountcx)
test_Data <- read.csv("TestDataCx.csv", header=TRUE)
head(test_Data,3)

test_Data$year <- as.factor(test_Data$year)
test_Data$month <- as.factor(test_Data$month)
##Check for NAs

#mosqcount <- mosqcountcx[complete.cases(mosqcountcx),]
#dim(mosqcount)
#head(mosqcount,3)

summary(mosqcountcx)

mosqcountcx$year <- as.factor(mosqcountcx$year)
mosqcountcx$month <- as.factor(mosqcountcx$month)

## Preliminary visualizations

par(mfrow=c(2,2))
plot(mosqcountcx$site, mosqcountcx$Cxpa)
plot(mosqcountcx$season, mosqcountcx$Cxpa)
plot(mosqcountcx$rainfall_season, mosqcountcx$Cxpa)
plot(mosqcountcx$elevation, mosqcountcx$Cxpa)

par(mfrow=c(2,2))
plot(mosqcount$precip, mosqcount$Cxpa)
plot(mosqcount$tmean, mosqcount$Cxpa)
plot(mosqcount$tmax, mosqcount$Cxpa)
plot(mosqcount$tmin, mosqcount$Cxpa)


names(mosqcountcx)

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


names(mosqcountcx)
dim(mosqcountcx)


preds<-c("Cxpa")
geo <- c("site","trap_type","month") ##"rainfall_season") ##"season",
loc <- c("elevation","distance")
precip <- c("precip","preciplag1","preciplag2") ##,"preciplag3"
temp <- c("tmean","tmax","tmin")
tempmean <- c("tmean","tmeanlag1","tmeanlag2") ##,"tmeanlag3"
tempmax <- c("tmax","tmaxlag1","tmaxlag2","tmaxlag3")
tempmin <- c("tmin","tminlag1","tminlag2","tminlag3")

all<- c(preds, geo, loc, precip, tempmean)


dataCX <- mosqcountcx[,all]
head(dataCX)
dim(dataCX)

head(mosqcountcx,2)

dataTrain <- mosqcountcx[,c(1,2,6,12,14,15,16,17,18)]
head(dataTrain,2)

names(test_Data)

dataTest <- test_Data[,c(1,2,6,12,14,15,16,17,18)]
head(dataTest,2)
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

#, type=2, extra=101, box.palette = "GnBu",branch.lty=3,shadow.col="gray",nn=TRUE)

## Find the lowest cp value, so the cross-validated error rate is minimum

CP.new <- f.null$cptable[which.min(f.null$cptable[,"xerror"]),"CP"]
CP.new


#f.null$cptable


#mm<-min(signif(f.null$cptable[,4],3))+mean(f.null$cptable[,5])
#mm

#w<-which(f.null$cptable[,4]<mm)
#w

#CP.new<-f.null$cptable[min(w)-1,1]
#CP.new

#0.01956152
#0.02066851
#0.02306952
## looking for 1 standard error from 

## now trimming based on the above
cnt<-rpart.control(minsplit=MS, cp=CP.new, xval=100)
f<-rpart(Cxpa ~  . , data=dataTrain, method="class", control=cnt)
plotcp(f)
graphics.off()


x11()
par(mfrow=c(1,1))

rpart.plot(f, type = 4, extra = 104, branch.lty = 3, gap = 0.6, space = 0.8, tweak = 1.2)


## Validation

library(ggplot2)
library(gmodels)
library(e1071)
library(gridExtra)
library(randomForest)


preds.rpart = predict(f,newdata = dataTest,type = "class")
CrossTable(dataTest$Cxpa,preds.rpart,chisq = F,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)

((2800 + 651)/nrow(dataTest))*100

PrecisionS <- 651/973 ; PrecisionS

RecallS <- 651/1145 ; RecallS

F1S <- 2 * ((67*57)/(67+57)) ; F1S

#plot(f, uniform=TRUE, margin=0.0075); text(f, digits=1, use.n=TRUE, cex=0.5)
#rpart.plot(f, type=4, extra=101, box.palette = "GnBu",branch.lty=3,shadow.col="gray",nn=TRUE)


#fancyRpartPlot(f, cex=0.4)

## here's what rpart gives for its trimmed tree if I don't include a
## controler

head(dataTrain,2)

par(mfrow=c(1,1))
f.d<-rpart(Cxpa ~  . , method="class", data=dataTrain)
printcp(f.d)
plotcp(f.d)

par(mfrow=c(1,1))
rpart.plot(f.d,type=4, extra=104, box.palette = "GnBu",branch.lty=3,shadow.col="gray")

#test <- predict(f, type = "vector")
#test
## squared bio variables

#dataCX$tmeansq<- dataCX$tmean^2
#dataCX$precipsq<- dataCX$precip^2
#
#head(dataCX)
#dim(dataCX)
#temps2<-c("tmeansq")
#temps2
#par(mfrow=c(1,1), bty="n")
#for(i in temps2) plot(dataCX[,i], dataCX$Cxpa, main=i)

#datapvasi$pv_pa<-as.integer(datapvasi$pv_pa)
#datapvasi$pv_pa


head(dataTrain)

##models


m.null<- glm(Cxpa ~ 1, family="binomial", data=dataTrain)
m.null
m.full<- glm(Cxpa ~ . + tmean * precip, family="binomial", data=dataTrain)
m.full


modelcxS<- step(m.null, scope=formula(m.full), direction="both", criterion = "BIC")
summary(modelcxS)

##validation

bestmodel <- glm(Cxpa ~ site + trap_type + month + tmean + preciplag2 + precip + preciplag1 + 
                   tmean*precip, family = "binomial", data = dataTrain)

summary(bestmodel)

predict <- predict(bestmodel, dataTest, type = 'response')
# confusion matrix
table_mat <- table(dataTest$Cxpa, predict > 0.5)
table_mat


accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test

precGLM <- 531/805 ; precGLM

recallGLM <- 531/1145 ; recallGLM

F1GLM <- 2 * ((66*46)/(66+46)) ; F1GLM

## Q-Q plots

qr2<- qresiduals(modelcxS)
length(qr2)
summary(qr2)


#y <- qr2[is.finite(qr2) & qr2 >= -5 & qr2 <= 5] 
#summary(y)
#length(y)

par(mfrow=c(1,1))
qqnorm(qr2, ylim = c(-6,6), xlim=c(-5,5), main = "Normal Q-Q Plot", las=1, bty="o"); qqline(qr2)

legend(3, 6.8, "B", cex = 2.5, box.lty = 0, bg = "transparent") 



modelcxS$fitted.values

bestmodel$fitted.values


## TEMPERATURE MEAN

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
#legend(10, 0.95, c("fitted values", "ave fitted"),
#col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



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
#legend(10, 0.95, c("fitted values", "ave fitted"),
# col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


#legend(x=27, y=1.2, legend="A", xpd = NA, bty="n", cex = 1.7)

legend(x=25, y=1.2, legend="D", xpd = NA, bty="n", cex = 2.2)


##PRECIPITATION


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
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,1200, by=1), y.m2$fit, col=2, lwd=3)
#legend(100, 0.95, c("fitted values", "ave fitted"),
#col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



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


##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=1), y.f2$fit, col=4, lwd=3)
#legend(100, 0.95, c("fitted values", "ave fitted"),
#col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=1050, y=1.2, legend="E", xpd = NA, bty="n", cex = 2.2)



names(dataTrain)

#########################
###Preciplag1####

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
#legend(100, 0.95, c("fitted values", "ave fitted"),
#col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



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
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
#legend(100, 0.95, c("fitted values", "ave fitted"),
#col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1050, y=1.2, legend="C", xpd = NA, bty="n", cex = 2.2)



## Precipitation lag2

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
#legend(100, 0.95, c("fitted values", "ave fitted"),
# col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



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
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
#legend(100, 0.95, c("fitted values", "ave fitted"),
# col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1050, y=1.2, legend="D", xpd = NA, bty="n", cex = 2.2)



##Distance
## 

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
#legend(200, 0.95, c("fitted values", "ave fitted"),
# col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



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
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 2500, by=100), y.f2$fit, col=4, lwd=3)
#legend(200, 0.95, c("fitted values", "ave fitted"),
# col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

##extract max value


legend(x=2150, y=1.2, legend="F", xpd = NA, bty="n", cex = 2.2)










## tmean lag2

o<-order(dataCX$tmeanlag2)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmeanlag2[o], dataCX$Cxpa[o], xlab="Temperature lag2 (°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature lag2")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmeanlag2=dataCX$tmeanlag2[o])
y.m2<-predict(loess(fit~tmeanlag2, data=fit.dat.m, span=1), data.frame(tmeanlag2=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmeanlag2[o]), dataCX$Cxpa[o],
     xlab="Temperature lag2 (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature lag2")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmeanlag2=dataCX$tmeanlag2[o])
fit.dat.f
y.f2<-predict(loess(fit~tmeanlag2, data=fit.dat.f, span=1), data.frame(tmeanlag2=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=28, y=1.2, legend="H", xpd = NA, bty="n", cex = 1.7)



## tmean lag3

o<-order(dataCX$tmeanlag3)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmeanlag3[o], dataCX$Cxpa[o], xlab="Temperature lag3 (°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature lag3")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmeanlag3=dataCX$tmeanlag3[o])
y.m2<-predict(loess(fit~tmeanlag3, data=fit.dat.m, span=1), data.frame(tmeanlag3=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 

##extract max value
ep.maxlgdpg <- findPeaks(y.m2$fit)
ep.maxlgdpg


temp <- seq(10, 28, by=0.1)
g <- data.frame(temp, y.m2$fit)
g

max(g$y.m2.fit, na.rm = TRUE)


## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmeanlag3[o]), dataCX$Cxpa[o],
     xlab="Temperature lag3 (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature lag3")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmeanlag3=dataCX$tmeanlag3[o])
fit.dat.f
y.f2<-predict(loess(fit~tmeanlag3, data=fit.dat.f, span=1), data.frame(tmeanlag3=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

##extract max value
ep.maxlgdp <- findPeaks(y.f2$fit)
ep.maxlgdp


temp <- seq(10, 28, by=0.1)
g <- data.frame(temp, y.f2$fit)
g

max(g$y.f2.fit, na.rm = TRUE)

legend(x=28, y=1.2, legend="L", xpd = NA, bty="n", cex = 1.7)













## Year
## lets look at marginal predictions, based on particular predictors

o<-order(dataCX$year)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd = TRUE,
  mfrow=c(1,2), bty="n")

ey <- expression(italic(Cx. ~ quinquefasciatus) ~ predicted ~ presence)

## plot based off predictions from the glm for data pre-2004
plot(dataCX$year[o], dataCX$Cxpa, xlab="Year",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(2002,2004),main="GLM fit - Year")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, year=dataCX$year[o])
fit.dat.m
y.m2<-predict(loess(fit~year, data=fit.dat.m, span = 1), data.frame(year=seq(2002, 2004, by=1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio03[o]/10, y.sym, col=2, lwd=3)
lines(seq(2002, 2004, by=1), y.m2$fit, col=2, lwd=3)
legend(2002, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree

f.pred<-predict(f, type="vector")
f.pred1 <- recode(f.pred, '1' = 0, '2' = 1)
table(f.pred1)


plot(dataCX$year[o], dataCX$Cxpa, ylim=c(0,1),xlim=c(2002,2004),
     xlab="Year", ylab= ey, cex.lab = 1.1,
     main="CART fit - Year")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, year=dataCX$year[o])
fit.dat.f
y.f2<-predict(loess(fit~year, data=fit.dat.f, span = 1), data.frame(year=seq(2002, 2004, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

x<-seq(2002,2004,1)
y<-y.f2$fit

lines(x,y)
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(x,y, col=4, lwd=3)
legend(2002, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=2004, y=1.2, legend="A", xpd = NA, bty="n", cex = 1.7)

names(dataCX)

## Month
## lets look at marginal predictions, based on particular predictors

o<-order(dataCX$month)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$month[o], dataCX$Cxpa, xlab="Month",
     ylab= ey, cex.lab=1.1, ylim=c(0,1), xlim=c(0,12),main="GLM fit - Month")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, month=dataCX$month[o])
y.m2<-predict(loess(fit~month, data=fit.dat.m, span = 1), data.frame(month=seq(0, 12, by=1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

#?loess
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
#y.m2$x <- seq(20, 45, by=0.1)
lines(seq(0, 12, by=1), y.m2$fit, col=2, lwd=3)
legend(1, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
f.pred1 <- recode(f.pred, '1' = 0, '2' = 1)
table(f.pred1)

plot(dataCX$month[o], dataCX$Cxpa,ylim=c(0,1),xlim=c(0,12),
     xlab="Month", ylab= ey, cex.lab=1.1, 
     main="CART fit - Month")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, month=dataCX$month[o])
fit.dat.f
y.f2<-predict(loess(fit~month, data=fit.dat.f, span = 1), data.frame(month=seq(0, 12, by= 1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio05[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 12, by=1), y.f2$fit, col=4, lwd=3)
legend(1, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=12, y=1.2, legend="B", xpd = NA, bty="n", cex = 1.7)



## elevation
## lets look at marginal predictions, based on particular predictors

summary(dataCX$elevation)

o<-order(dataCX$elevation)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(dataCX$elevation[o], dataCX$Cxpa, xlab="Elevation (m)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,2000),main="GLM fit - Elevation")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, elevation=dataCX$elevation[o])
y.m2<-predict(loess(fit~elevation, data=fit.dat.m, span = 1), data.frame(elevation=seq(0, 2000, by=10)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,2000, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
table(f.pred)


#f.pred[f.pred = 1] <- 0 
f.pred1 <- recode(f.pred, '1' = 0, '2' = 1)
table(f.pred1)

plot(jitter(dataCX$elevation[o]), dataCX$Cxpa, ylim=c(0,1),xlim=c(0,2000),
     xlab="Elevation (m)", ylab= ey, cex.lab = 1.1,
     main="CART fit - Elevation")
fittedf<- as.numeric(f.pred1[o])
fit.dat.f<-data.frame(fit=fittedf, elevation=dataCX$elevation[o])
#fit.dat.f
#plot(fit.dat.f$elevation, fit.dat.f$fit)
y.f2<-predict(loess(fit~elevation, data=fit.dat.f, span=1), data.frame(elevation=seq(0, 2000, by=10)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 2000, by=10), y.f2$fit, col=4, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1900, y=1.2, legend="A", xpd = NA, bty="n", cex = 1.7)
























