## R code for "effects of temperature on the life-history traits of Myzus persicae and its efficiency in transmitting PVY virus in potato crops"
## Authors: Bonoukpoe Sokame, Henri Tonnaang, Heidy Gamarra, Pablo carhuapoma, Jan Kreuze, 
## Ali Arab, Peter Armbruster, Leah Johnson, and Oswaldo Villena.


##Load packages

library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(readxl)  
library(ggsurvfit)
library(emmeans)
library(dplyr)

##Load dataset

multiplesheets <- function(fname) { 
  
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  # assigning names to data frames 
  names(data_frame) <- sheets 
  
  # print data frame 
  print(data_frame) 
} 

# specifying the path name 
path <- "./data/Mpersicae data.xls"
dat <- multiplesheets(path)


## 1. Kaplan-Meier non-parametric analysis

## Survival analysis 

## A. Nymph survival
NymphSurv <- dat$NymphSurvival
NymphSurv
attach(NymphSurv)

time<- Time
event<- Status
group<- Treatment
summary(time)
summary(event)
summary(group)
kmsurvival1<- survfit(Surv(time,event)~group, data = NymphSurv)

ggsurvplot(kmsurvival1, xlim=c(0,30), break.x.by=5, ylab="Survival (rate)", xlab="Time (days)", 
           pval = F, risk.table = FALSE, risk.table.title="", legend.labs=c("10\u00B0C","15\u00B0C","20\u00B0C","25\u00B0C","30\u00B0C"), legend=c(0.7,0.8), legend.title="",
           surv.scale="default", break.y.by=0.2, censor.shape=16, censor.size=2,
           palette=c("yellowgreen","black","orange","blue2","red"), title="(A) Nymph Survival")


## B. Adult survival
AdultSurv <- dat$AdultSurvival
attach(AdultSurv)

time<- Time
event<- Status
group<- Treatment
summary(time)
summary(event)
summary(group)
kmsurvival2<- survfit(Surv(time,event)~group, data = AdultSurv)

ggsurvplot(kmsurvival2, xlim=c(0,25), break.x.by=5, ylab="Survival (rate)", xlab="Time (days)",
           pval = F, risk.table = FALSE, risk.table.title="", legend.labs=c("10\u00B0C","15\u00B0C","20\u00B0C","25\u00B0C","30\u00B0C"), legend=c(0.7,0.8), legend.title="",
           surv.scale="default", break.y.by=0.2, censor.shape=16, censor.size=2,
           palette=c("yellowgreen","black","orange","blue","red"), title="(B) Adult Survival")


## C. Nymph life_expected

nymphExp <- dat$ExpectLnymph
nymphExp$Temperature <- as.factor(nymphExp$Temperature)

p <- ggplot(nymphExp, aes(x = Time_days, y = Expectedlife))
p + geom_line() + geom_point() + aes(color=Temperature) +
  scale_color_manual(values = c("yellowgreen", "black","orange","blue","red"))+
  scale_y_continuous(breaks=seq(0, 16, 2)) +
  xlim(0,25)


## D. Adult life_expected

adultexp <- dat$ExpectLadult
adultexp$Temperature <- as.factor(adultexp$Temperature)

p <- ggplot(adultexp, aes(x = Time_days, y = Expectedlife))
p + geom_line() + geom_point() + aes(color=Temperature) +
  scale_color_manual(values = c("yellowgreen", "black","orange","blue","red"))+
  scale_y_continuous(breaks=seq(0, 16, 2)) +
  xlim(0,25)

###############################################################################################################

## 2. Count data (nymph longevity, adult longevity, total longevity, and fecundity expressed as offspring per female) were analysed with generalized linear model (GLM) with negative binomial error distribution considering overdispersion

## Load packages

library(MASS)
library(lme4)
library(car)
library(emmeans)
library(multcompView)
library(lsmeans)


## longevity
longevity <- dat$Longevity

longevity$Treatment <- as.factor(longevity$Treatment)

#mod1<- glm.nb(dependent-variable~Temperature,init.theta = 2, link = log)
mod1 <- glm.nb(Female_longevity~Treatment, data=longevity, init.theta = 2, link = log)
summary(mod1)
anova(mod1, test = "Chisq")
ls_means <- emmeans(mod1, ~ Treatment)
pairs(ls_means, adjust = "tukey")


##fecundity
fecundity <- dat$Fecundity
fecundity$Treatment <- as.factor(fecundity$Treatment)

mod2 <- glm.nb(Fecundity~Treatment, data=fecundity, init.theta = 2, link = log)
summary(mod2)
anova(mod2, test = "Chisq")
ls_means <- emmeans(mod2, ~ Treatment)
pairs(ls_means, adjust = "tukey")


###############################################################################################################

###3. Life table parameter data were analysed using Kruskal-Wallis nonparametric procedure

##Load packages
library(dunn.test)

##Load data
parameters <- dat$`Life table parameters`
parameters
attach(parameters)

## Net reproduction rate (Ro)
shapiro.test(parameters$Ro)
bartlett.test(parameters$Ro ~ Treatment)
dunn.test(parameters$Ro, Treatment, wrap = T)

## Intrinsic rate of increase (rm)
shapiro.test(parameters$rm)
bartlett.test(parameters$rm ~ Treatment)
dunn.test(parameters$rm, Treatment, wrap = T)

## Finite rate of increase (lambda)
shapiro.test(parameters$lambda)
bartlett.test(parameters$lambda ~ Treatment)
dunn.test(parameters$lambda, Treatment, wrap = T)
