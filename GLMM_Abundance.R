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
simulateResiduals(fittedModel = zinbmd7, plot = T)


