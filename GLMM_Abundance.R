##author: Oswaldo Villena
##Load package
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)
library(glmmTMB)
library(knitr)
library(bbmle)
library(reshape)
library(dplyr)
library(DHARMa)
library(modEvA)
library(performance)
library(MASS)

getwd()

##Load data
mosquicountcx <-read.csv("mosquitodataculex.csv", header=TRUE)
head(mosquicountcx,4)
dim(mosquicountcx)
names(mosquicountcx)

##Be sure year and month are factors
mosquicountcx$yearf <- as.factor(mosquicountcx$year) 
mosquicountcx$monthf <- as.factor(mosquicountcx$month)

##Add a variable for square temperature
mosquicountcx$sqtmean <- mosquicountcx$tmean^2

##Select your variables of interest
mosquicount <- mosquicountcx[c(1,3,6,7,12,15:18,20,32,33,34)]
head(mosquicount)

## Poisson models

pmd1 <- glmmTMB(count ~ site + (1|yearf), mosquicount, family = poisson)
pmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), mosquicount, family = poisson)
pmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), mosquicount, family = poisson)
pmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|year), mosquicount, family = poisson)
pmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|year), mosquicount, family = poisson)
pmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|year), 
                mosquicount, family = poisson)
pmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|year), 
                mosquicount, family = poisson)
pmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + 
                  (1|year), mosquicount, family = poisson)
pmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 +
                  (1|year), mosquicount, family = poisson)
pmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 +
                   preciplag2 + (1|year), mosquicount, family = poisson)

##Estimate AIC values
AIC(pmd1, pmd2, pmd3, pmd4, pmd5, pmd6, pmd7, pmd8, pmd9, pmd10)


##Negative Binomial models

nbmd1 <- glmmTMB(count ~ site + (1|yearf), mosquicount, family = nbinom2)
nbmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), mosquicount, family = nbinom2)
nbmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), mosquicount, family = nbinom2)
nbmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|year), mosquicount, family = nbinom2)
nbmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|year), mosquicount, family = nbinom2)
nbmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|year), mosquicount, family = nbinom2)
nbmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|year), 
                 mosquicount, family = nbinom2)
nbmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + 
                   (1|year), mosquicount, family = nbinom2)
nbmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 +
                   (1|year), mosquicount, family = nbinom2)
nbmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 +
                    preciplag2 + (1|year), mosquicount, family = nbinom2)

##Estimate AIC values
AIC(nbmd1, nbmd2, nbmd3, nbmd4, nbmd5, nbmd6, nbmd7, nbmd8, nbmd9, nbmd10)


## NEGATIVE BINOMIAL WITH DISPERSION

nbdmd1 <- glmmTMB(count ~ site + (1|yearf), dispformula = ~monthf, mosquicount, family = nbinom2)
nbdmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), dispformula = ~monthf, mosquicount, family = nbinom2)
nbdmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), dispformula = ~monthf, mosquicount, family = nbinom2)
nbdmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|year), dispformula = ~monthf, 
                  mosquicount, family = nbinom2)
nbdmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|year), dispformula = ~monthf,
                  mosquicount, family = nbinom2)
nbdmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|year), dispformula = ~monthf,
                  mosquicount, family = nbinom2)
nbdmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|year), dispformula = ~monthf,
                  mosquicount, family = nbinom2)
nbdmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + 
                    (1|year), dispformula = ~monthf, mosquicount, family = nbinom2)
nbdmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 +
                    (1|year), dispformula = ~monthf, mosquicount, family = nbinom2)
nbdmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 +
                     preciplag2 + (1|year), dispformula = ~monthf, mosquicount, family = nbinom2)


##Estimate AIC values
AIC(nbdmd1, nbdmd2, nbdmd3, nbdmd4, nbdmd5, nbdmd6, nbdmd7, nbdmd8, nbdmd9, nbdmd10)


## ZERO INFLATED POISSON MODELS

zipmd1 <- glmmTMB(count ~ site + (1|year), zi=~site, mosquicount, family = poisson)
zipmd2 <- glmmTMB(count ~ site + month + (1|year), zi=~site + month, mosquicount, family = poisson)
zipmd3 <- glmmTMB(count ~ site + month + trap_type + (1|year), zi=~site + month + trap_type, mosquicount, family = poisson)
zipmd4 <- glmmTMB(count ~ site + month + trap_type + tmean + (1|year), zi=~site + 
                    month + trap_type + tmean, mosquicount, family = poisson)
zipmd5 <- glmmTMB(count ~ site + month + trap_type + tmean + precip + (1|year), zi=~site + 
                    month + trap_type + tmean + precip, mosquicount, family = poisson)
zipmd6 <- glmmTMB(count ~ site + month + trap_type + tmean + precip + distance + (1|year), zi=~site + 
                    month + trap_type + tmean + precip + distance, mosquicount, family = poisson)
zipmd7 <- glmmTMB(count ~ site + month + trap_type + tmean + sqtmean + precip + distance + (1|year), zi=~site + 
                    month + trap_type + tmean + precip + distance, mosquicount, family = poisson)
zipmd8 <- glmmTMB(count ~ site + month + trap_type + tmean + precip + distance + preciplag1 + (1|year), zi=~site + 
                    month + trap_type + tmean + precip + distance + preciplag1, mosquicount, family = poisson)
zipmd9 <- glmmTMB(count ~ site + month + trap_type + tmean + precip + distance + preciplag1 + preciplag2 + (1|year), zi=~site + 
                    month + trap_type + tmean + precip + distance + preciplag1 + preciplag2, mosquicount, family = poisson)
zipmd10 <- glmmTMB(count ~ site + month + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2 + (1|year),
                   zi=~site + month + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2,
                   mosquicount, family = poisson)

##Estimate AIC values
AIC(zipmd1, zipmd2, zipmd3, zipmd4, zipmd5, zipmd6, zipmd7, zipmd8, zipmd9, zipmd10)



## ZERO INFLATED NEGATIVE BINOMIAL MODELS

zinbmd1 <- glmmTMB(count ~ site + (1|yearf), zi=~site, mosquicount, family = nbinom2)
zinbmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), zi=~site + monthf, mosquicount, family = nbinom2)
zinbmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), zi=~site + monthf + trap_type, mosquicount, family = nbinom2)
zinbmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean, mosquicount, family = nbinom2)
zinbmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean + precip, mosquicount, family = nbinom2)
zinbmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean + precip + distance, mosquicount, family = nbinom2)
zinbmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean + sqtmean + precip + distance, mosquicount, family = nbinom2)
zinbmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean + precip + distance + preciplag1, mosquicount, family = nbinom2)
zinbmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 + (1|yearf), zi=~site + 
                     monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2, mosquicount, family = nbinom2)
zinbmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2 + (1|yearf),
                    zi=~site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2,
                    mosquicount, family = nbinom2)

##Estimate AIC values
AIC(zinbmd1, zinbmd2, zinbmd3, zinbmd4, zinbmd5, zinbmd6, zinbmd7, zinbmd8, zinbmd9, zinbmd10)


## HURDLE MODELS POISSON

hpmd1 <- glmmTMB(count ~ site + (1|yearf), zi=~site, mosquicount, family = truncated_poisson)
hpmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), zi=~site + monthf, mosquicount, family = truncated_poisson)
hpmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), zi=~site + monthf + trap_type, mosquicount,
                 family = truncated_poisson)
hpmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean, mosquicount, family = truncated_poisson)
hpmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean + precip, mosquicount, family = truncated_poisson)
hpmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean + precip + distance, mosquicount, family = truncated_poisson)
hpmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean + sqtmean + precip + distance, mosquicount, family = truncated_poisson)
hpmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean + precip + distance + preciplag1, mosquicount, family = truncated_poisson)
hpmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 + (1|yearf), zi=~site + 
                   monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2, mosquicount, 
                 family = truncated_poisson)
hpmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2 + (1|yearf),
                  zi=~site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2,
                  mosquicount, family = truncated_poisson)

##Estimate AIC values
AIC(hpmd1, hpmd2, hpmd3, hpmd4, hpmd5, hpmd6, hpmd7, hpmd8, hpmd9, hpmd10)


## HURDLE MODELS NEGATIVE BINOMIALS

hnbmd1 <- glmmTMB(count ~ site + (1|yearf), zi=~site, mosquicount, family = truncated_nbinom2)
hnbmd2 <- glmmTMB(count ~ site + monthf + (1|yearf), zi=~site + monthf, mosquicount, family = truncated_nbinom2)
hnbmd3 <- glmmTMB(count ~ site + monthf + trap_type + (1|yearf), zi=~site + monthf + trap_type, mosquicount,
                  family = truncated_nbinom2)
hnbmd4 <- glmmTMB(count ~ site + monthf + trap_type + tmean + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean, mosquicount, family = truncated_nbinom2)
hnbmd5 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean + precip, mosquicount, family = truncated_nbinom2)
hnbmd6 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean + precip + distance, mosquicount, family = truncated_nbinom2)
hnbmd7 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean + sqtmean + precip + distance, mosquicount, family = truncated_nbinom2)
hnbmd8 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean + precip + distance + preciplag1, mosquicount, family = truncated_nbinom2)
hnbmd9 <- glmmTMB(count ~ site + monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2 + (1|yearf), zi=~site + 
                    monthf + trap_type + tmean + precip + distance + preciplag1 + preciplag2, mosquicount, 
                  family = truncated_nbinom2)
hnbmd10 <- glmmTMB(count ~ site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2 + (1|yearf),
                   zi=~site + monthf + trap_type + tmean + sqtmean + precip + distance + preciplag1 + preciplag2,
                   mosquicount, family = truncated_nbinom2)

##Estimate AIC values
AIC(hnbmd1, hnbmd2, hnbmd3, hnbmd4, hnbmd5, hnbmd6, hnbmd7, hnbmd8, hnbmd9, hnbmd10)

## Estimate Conditional and Marginal R2 (this piece of code could take between 10-30h)
r2(pmd1);r2(pmd2);r2(pmd3);r2(pmd4);r2(pmd5);r2(pmd6)
r2(pmd7);r2(pmd8);r2(pmd9);r2(pmd10)
r2(nbmd1);r2(nbmd2);r2(nbmd3);r2(nbmd4);r2(nbmd5);r2(nbmd6)
r2(nbmd7);r2(nbmd8);r2(nbmd9);r2(nbmd10)
r2(nbdmd1);r2(nbdmd2);r2(nbdmd3);r2(nbdmd4);r2(nbdmd5);r2(nbdmd6)
r2(nbdmd7);r2(nbdmd8);r2(nbdmd9);r2(nbdmd10)
r2(zipmd1);r2(zipmd2);r2(zipmd3);r2(zipmd4);r2(zipmd5);r2(zipmd6)
r2(zipmd7);r2(zipmd8);r2(zipmd9);r2(zipmd10)
r2(zinbmd1);r2(zinbmd2);r2(zinbmd3);r2(zinbmd4);r2(zinbmd5);r2(zinbmd6)
r2(zinbmd7);r2(zinbmd8);r2(zinbmd9);r2(zinbmd10)
r2(hpmd1);r2(hpmd2);r2(hpmd3);r2(hpmd4);r2(hpmd5);r2(hpmd6)
r2(hpmd7);r2(hpmd8);r2(hpmd9);r2(hpmd10)
r2(hnbmd1);r2(hnbmd2);r2(hnbmd3);r2(hnbmd4);r2(hnbmd5);r2(hnbmd6)
r2(hnbmd7);r2(hnbmd8);r2(hnbmd9);r2(hnbmd10)

## VISUALITATION OF BEST MODEL

##Summary best model
summary(zinbmd7)

##Residuals best model: QQ plot residuals and residuals vs predicted
simulationOutput <- simulateResiduals(fittedModel = zinbmd7, plot = T)
simulationOutput

##Check that the model is able to model the zero count portion of the dataset
testZeroInflation(simulationOutput)

##Generate new data to check and to plot model output
newdata0 = newdata = unique(mosquicount[,c("site","monthf","trap_type","tmean","sqtmean","precip","distance")])
head(newdata0,20)

## Model predictions

X.cond <- model.matrix(lme4::nobars(formula(zinbmd7)[-2]), newdata0)
head(X.cond)
beta.cond <- fixef(zinbmd7)$cond
head(beta.cond)
pred.cond <- X.cond %*% beta.cond
head(pred.cond)


ziformula = zinbmd7$modelInfo$allForm$ziformula
X.zi <- model.matrix(lme4::nobars(ziformula),newdata0)
head(X.zi)
beta.zi <- fixef(zinbmd7)$zi
head(beta.zi)
pred.zi <- X.zi %*% beta.zi
head(pred.zi)


##estimates of the linear predictor

pred.ucount <- exp(pred.cond)*(1-plogis(pred.zi))
head(pred.ucount,10)

##estimates of standard errors and CIs

pred.condpar.psim <- mvrnorm(1000, mu=beta.cond, Sigma = vcov(zinbmd7)$cond) 
head(pred.condpar.psim)
pred.cond.psim <- X.cond %*% t(pred.condpar.psim)
pred.zipar.psim <- mvrnorm(1000, mu=beta.zi, Sigma =vcov(zinbmd7)$zi)
pred.zi.psim <- X.zi %*% t(pred.zipar.psim)
pred.ucount.psim <- exp(pred.cond.psim)*(1-plogis(pred.zi.psim))

ci.ucount <- t(apply(pred.ucount.psim,1,quantile,c(0.025,0.975))) ## alpha = 0.05
head(ci.ucount)
ci.ucount <- data.frame(ci.ucount)
names(ci.ucount) <- c("ucount.low","ucount.high")
head(ci.ucount)
pred.ucount <- data.frame(newdata0, pred.ucount, ci.ucount)
head(pred.ucount,20)

##Mean and median observed mosquito count

library(plyr)

real.count <- ddply(mosquicount, ~site+monthf+trap_type+distance+tmean+precip+sqtmean,
                    summarize, m=median(count), mu=mean(count))

##expression for italic style for plots
ylb <- expression(Predicted~italic(Cx.~quinquefasciatus)~abundance)

##Predicted vs Observed mosquito abundance by month and sites 

##By trap type and month
ggplot(pred.ucount,aes(x=monthf, y=pred.ucount, colour=trap_type)) + geom_point(shape=1, size=2) +
  geom_errorbar(aes(ymin=ucount.low, ymax=ucount.high))+
  geom_point(data=real.count, aes(x=monthf, y=m, colour=trap_type), shape=20, size=2)+
  ylab(ylb)+
  xlab("Months") +
  scale_x_discrete(labels=month.abb)

##By trap type, month, and site
ggplot(pred.ucount,aes(x=site, y=pred.ucount, colour=trap_type)) + geom_point(shape=1, size=2) +
  geom_errorbar(aes(ymin=ucount.low, ymax=ucount.high))+
  geom_point(data=real.count, aes(x=site, y=m, colour=trap_type), shape=20, size=2)+
  facet_wrap(~monthf) +
  ylab(ylb)+
  xlab("Sites") 



## abundamce vs month

meancountS <- pred.ucount %>%
  dplyr::group_by(monthf, trap_type) %>%
  summarise(estimated = mean(pred.ucount, na.rm = T), cilow = mean(ucount.low, na.rm = T),
            cihigh = mean(ucount.high, na.rm = T))

abundco2S <-meancountS[which(meancountS$trap_type=="co2"), ]  ##CO2 traps

abundgravidS <-meancountS[which(meancountS$trap_type=="gravid"), ]  ##Gravid traps


##REAL COUNT

realcountS <- real.count %>%
  dplyr::group_by(monthf, trap_type) %>%
  summarise(observed = mean(m, na.rm = T))

realco2S <-realcountS[which(realcountS$trap_type=="co2"), ]  ## CO2 traps

realgravidS <-realcountS[which(realcountS$trap_type=="gravid"), ]  ##Gravid traps


##Abundance by month plots

library(ggplot2)
library(RCurl)

##For CO2 traps
ggplot(abundco2S,aes(x=monthf, y=estimated)) + geom_point(shape=20, size=4, color = "black") +
  geom_errorbar(aes(ymin=cilow, ymax=cihigh), width = 0.2, size = 1.2, color = "blue") +
  geom_point(data=realco2S, aes(x=monthf, y=observed), shape=20, size=4, color = "red")+
  labs(x="Months", y=ylb, color = "Legend") +
  scale_color_manual("",
                     breaks = c("Estimated","Observed"),
                     values = c("Estimated" = "black", "Observed" = "red"))


##For gravid traps
ggplot(abundgravidS,aes(x=monthf, y=estimated)) + geom_point(shape=20, size=4, color = "black") +
  geom_errorbar(aes(ymin=cilow, ymax=cihigh), width = 0.2, size = 1.2, color = "blue") +
  geom_point(data=realgravidS, aes(x=monthf, y=observed), shape=20, size=4, color = "red")+
  labs(x="Months", y=ylb, color = "Legend") +
  scale_color_manual("",
                     breaks = c("Estimated","Observed"),
                     values = c("Estimated" = "black", "Observed" = "red"))



## BY site 

MeanAbund <- pred.ucount %>%
  dplyr::group_by(site, trap_type) %>%
  summarise(estimated = mean(pred.ucount, na.rm = T), cilow = mean(ucount.low, na.rm = T),
            cihigh = mean(ucount.high, na.rm = T))

##Split data set between CO2 and gravid traps

PAbundCO2 <- MeanAbund[which(MeanAbund$trap_type=="co2"), ]

PAbundGravid <-MeanAbund[which(MeanAbund$trap_type=="gravid"), ]


##Observed counts

Realcount <- real.count %>%
  dplyr::group_by(site, trap_type) %>%
  summarise(observed = mean(m, na.rm = T))

##Split by trap type

RAbundCO2 <- Realcount[which(Realcount$trap_type=="co2"), ]

RAbundGravid <- Realcount[which(Realcount$trap_type=="gravid"), ]


## PLot abundance by site for CO2 traps
ylb <- expression(Average~italic(Cx.~quinquefasciatus)~abundance~-~CO2~traps) ##expression italics for label

ggplot(PAbundCO2,aes(x=site, y=estimated)) + geom_point(shape=20, size=4, color = "black") +
  geom_errorbar(aes(ymin=cilow, ymax=cihigh), width = 0.2, size = 1.2, color = "blue") +
  geom_point(data=RAbundCO2, aes(x=site, y=observed), shape=20, size=4, color = "red")+
  labs(x="Sites", y=ylb, color = "Legend") +
  scale_color_manual("",
                     breaks = c("Estimated","Observed"),
                     values = c("Estimated" = "black", "Observed" = "red"))


## PLot abundance by site for gravid traps
ylb <- expression(Average~italic(Cx.~quinquefasciatus)~abundance~-~gravid~traps) ##expression italics for label

ggplot(PAbundGravid,aes(x=site, y=estimated)) + geom_point(shape=20, size=4, color = "black") +
  geom_errorbar(aes(ymin=cilow, ymax=cihigh), width = 0.2, size = 1.2, color = "blue") +
  geom_point(data=RAbundGravid, aes(x=site, y=observed), shape=20, size=4, color = "red")+
  labs(x="Sites", y=ylb, color = "Legend") +
  scale_color_manual("",
                     breaks = c("Estimated","Observed"),
                     values = c("Estimated" = "black", "Observed" = "red"))

