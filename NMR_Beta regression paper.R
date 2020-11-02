# Beta Regression Analysis 

# set working directory
setwd("~/Desktop/Beta Regression Paper/Beta Regression R analysis")

# load data
SoilNMR<-read.csv("NMR data.csv")

# load packages
library(betareg)
library(car)
library(gamlss)
library(emmeans)
library(dplyr)
library(lmtest)
library(glmmTMB)
library(DHARMa)
library(rcompanion)

------------------------------------CARBONYL CARBON ANALYSES------------------------------------

# ----- Soil NMR -----
# explore data
head(SoilNMR)
summary(SoilNMR)

# General Linear Model --- Carbonyl Carbon
m1.SoilNMR<-glm(CC~Region+Layer+Region*Layer,data=na.omit(SoilNMR),
             family = gaussian(link="identity"))
plot(m1.SoilNMR)
hist(resid(m1.SoilNMR))

# Plot for manuscript - Carbonyl Carbon GLM
png("CarbonylCarbon_DiagnosticPlotsGLM.png", width = 160, height = 160, units = "mm", res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m1.SoilNMR),y=resid(m1.SoilNMR),main=NULL,xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A", side = 2, line = 2, at = 0.009, col = 'black', font = 2, las = 1)
hist(resid(m1.SoilNMR),main=NULL,xlab="Residuals", cex.lab=1.15)
mtext("B", side = 2, line = 2, at = 10.3, col = 'black', font = 2, las = 1)
qqnorm(resid(m1.SoilNMR),main=NULL, cex.lab=1.15)
mtext("C", side = 2, line = 2, at = 0.009, col = 'black', font = 2, las = 1)
qqline(resid(m1.SoilNMR),col='red')
plot(m1.SoilNMR,which=4,main=NULL, cex.lab=1.15)
mtext("D", side = 2, line = 2, at = .49, col = 'black', font = 2, las = 1)
dev.off()

# Model summary/diagnostics
summary(m1.SoilNMR)
exp(logLik(m1.SoilNMR))
anova(m1.SoilNMR)
#ANODEV
m1.intercept<-glm(CC~1,data=na.omit(SoilNMR),
                  family=gaussian(link = "identity"))
m1.OHCC<-glm(CC~1+Region,data=na.omit(SoilNMR),
              family=gaussian(link = "identity"))
m1.CR<-glm(CC~1+Region+Layer,data=na.omit(SoilNMR),
           family=gaussian(link = "identity"))
m1.OHCCCR<-glm(CC~1+Region+Layer+Region*Layer,data=na.omit(SoilNMR),
              family=gaussian(link = "identity"))
m1ANODEV<-lrtest(m1.intercept,m1.OHCC,m1.CR,m1.OHCCCR)
m1ANODEV

# Regular beta regressions {betareg}
library(betareg)
m3<-betareg(CC~Region+Layer+Region*Layer,data=SoilNMR)
plot(m3)

# Plots for manuscript - Carbonyl Carbon Beta
png("CarbonylCarbon_DiagnosticPlotsBeta.png",  width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(m3,which=1,type="pearson",cex.lab=1.15, caption=NULL)
mtext("A",side=2,line=2,at=2.5,col='black',font=2,las=1)
plot(m3,which=4,type="pearson",cex.lab=1.15, caption=NULL)
mtext("B",side=2,line=2,at=2.5,col='black',font=2,las=1)
plot(m3,which=5,type="deviance",cex.lab=1.15, caption=NULL)
mtext("C",side=2,line=2,at=3.7,col='black',font=2,las=1)
plot(m3,which=2,type="pearson", cex.lab=1.15, caption=NULL)
mtext("D",side=2,line=2,at=0.27,col='black',font=2,las=1)
dev.off()

# model summary
summary(m3)
#ANODEV
m3.intercept<-betareg(CC~1,data=SoilNMR)
m3.OHCC<-betareg(CC~1+Region,data=SoilNMR)
m3.CR<-betareg(CC~1+Region+Layer,data=SoilNMR)
m3.OHCCCR<-betareg(CC~1+Region+Layer+Layer*Region,data=SoilNMR)

install.packages("lrtest")
m3ANODEV<-lrtest(m3.intercept,m3.OHCC,m3.CR,m3.OHCCCR)
m3ANODEV # same result as the GLMMTMB

# ---- Model comparison: SoilNMR ----

glm.model.SoilNMR.CC<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df')

CCglm.SoilNMR<-glm.model.SoilNMR.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

beta.model.SoilNMR.CC<-as.data.frame(m3ANODEV)%>%
  rename(numDf='#Df')

CCbeta.SoilNMR<-beta.model.SoilNMR.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(CCglm.SoilNMR,file ="CCglm.model.csv",row.names=FALSE)
write.csv(CCbeta.SoilNMR,file ="CCbeta.model.csv",row.names = FALSE)

# ---- Model parameters: SoilNMR ----
m1.SoilNMR$coefficients
summary(m1.SoilNMR)
m1.SoilNMR$coefficients

summary(m3)
betacoef<-m3$mu.coefficients
odds<-exp(coef(m3))
beta.prob<-odds/(1+odds)

parameter.summary.CC<-data.frame(glm=m1.SoilNMR$coefficients,
                              beta=beta.prob[-c(13)])

# save table
write.csv(parameter.summary.CC,"CCparameters.csv",row.names = FALSE)

# ---- Fit statistics: SoilNMR -----

#Calculate r-squared for LR using Coz and Snell
library(rcompanion)
nagelkerke(m1.SoilNMR) #R^2 for GLM is the Cox and Snell R-squared value
nagelkerke(m3) #R^2 for Beta is the Cox and Snell R-squared value

fit.SoilNMR<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                     glm=c(AIC(m1.SoilNMR),(1-0.9926)^(-36/2),
                           sum(m1.SoilNMR$residuals^2)/m1.SoilNMR$df.residual),
                     beta=c(AIC(m3),(1-0.9851)^(-36/2),
                            sum(m3$residuals^2)/m3$df.residual))

# save table
write.csv(fit.SoilNMR,"Carbonylcarbonfit.stats.csv",row.names = FALSE)

-------------------------------------ALKYL CARBON ANALYSES------------------------------------

  # ----- Soil NMR -----
# explore data
head(SoilNMR)
summary(SoilNMR)

# General Linear Model --- Alkyl Carbon
m4.SoilNMR<-glm(AC~Region+Layer+Region*Layer,data=na.omit(SoilNMR),
                family = gaussian(link="identity"))
plot(m4.SoilNMR)
hist(resid(m4.SoilNMR))

# model summary/diagnostics
summary(m4.SoilNMR)
exp(logLik(m4.SoilNMR))
anova(m4.SoilNMR)

#ANODEV
m4.intercept<-glm(AC~1,data=na.omit(SoilNMR),
                  family=gaussian(link = "identity"))
m4.OHCC<-glm(AC~1+Region,data=na.omit(SoilNMR),
             family=gaussian(link = "identity"))
m4.CR<-glm(AC~1+Region+Layer,data=na.omit(SoilNMR),
           family=gaussian(link = "identity"))
m4.OHCCCR<-glm(AC~1+Region+Layer+Region*Layer,data=na.omit(SoilNMR),
               family=gaussian(link = "identity"))
m4ANODEV<-lrtest(m4.intercept,m4.OHCC,m4.CR,m4.OHCCCR)
m4ANODEV

# Plot for manuscript - Alkyl Carbon GLM
png("AlkylCarbon_DiagnosticPlotsGLM.png", width = 160, height = 160, units = "mm", res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m4.SoilNMR),y=resid(m1.SoilNMR),main=NULL,xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A", side = 2, line = 2, at = 0.009, col = 'black', font = 2, las = 1)
hist(resid(m4.SoilNMR),main=NULL,xlab="Residuals", cex.lab=1.15)
mtext("B", side = 2, line = 2, at = 14.5, col = 'black', font = 2, las = 1)
qqnorm(resid(m4.SoilNMR),main=NULL, cex.lab=1.15)
mtext("C", side = 2, line = 2, at = 0.06, col = 'black', font = 2, las = 1)
qqline(resid(m4.SoilNMR),col='red')
plot(m4.SoilNMR,which=4,main=NULL, cex.lab=1.15)
mtext("D", side = 2, line = 2, at = .53, col = 'black', font = 2, las = 1)
dev.off()

# Regular beta regressions {betareg}
library(betareg)
m6.SoilNMR<-betareg(AC~Region+Layer+Region*Layer,data=SoilNMR)
plot(m6.SoilNMR)

# model summary
summary(m6.SoilNMR)
#ANODEV
m6.intercept<-betareg(AC~1,data=SoilNMR)
m6.OHCC<-betareg(AC~1+Region,data=SoilNMR)
m6.CR<-betareg(AC~1+Region+Layer,data=SoilNMR)
m6.OHCCCR<-betareg(AC~1+Region+Layer+Region*Layer,data=SoilNMR)

install.packages("lrtest")
m6ANODEV<-lrtest(m6.intercept,m6.OHCC,m6.CR,m6.OHCCCR)
m6ANODEV # same result as the GLMMTMB

# Plots for manuscript - Alkyl Carbon Beta
png("AlkylCarbon_DiagnosticPlotsBeta.png",  width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(m6.SoilNMR,which=1,type="pearson",cex.lab=1.15, caption=NULL)
mtext("A",side=2,line=2,at=3,col='black',font=2,las=1)
plot(m6.SoilNMR,which=4,type="pearson",cex.lab=1.15, caption=NULL)
mtext("B",side=2,line=2,at=3,col='black',font=2,las=1)
plot(m6.SoilNMR,which=5,type="deviance",cex.lab=1.15, caption=NULL)
mtext("C",side=2,line=2,at=3.6,col='black',font=2,las=1)
plot(m6.SoilNMR,which=2,type="pearson", cex.lab=1.15, caption=NULL)
mtext("D",side=2,line=2,at=0.44,col='black',font=2,las=1)
dev.off()

# ---- Model comparison: SoilNMR ----
glm.model.SoilNMR.AC<-as.data.frame(m4ANODEV)%>%
  rename(numDf='#Df')

glm.SoilNMR<-glm.model.SoilNMR.AC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

beta.model.SoilNMR.AC<-as.data.frame(m6ANODEV)%>%
  rename(numDf='#Df')
beta.SoilNMR<-beta.model.SoilNMR.AC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(glm.SoilNMR,file ="AC.glm.model.csv",row.names=FALSE)
write.csv(beta.SoilNMR,file ="AC.beta.model.csv",row.names = FALSE)

# ---- Model parameters: SoilNMR ----
m4.SoilNMR$coefficients
summary(m4.SoilNMR)
m4.SoilNMR$coefficients

summary(m6.SoilNMR)
betacoef.AC<-m6.SoilNMR$mu.coefficients
odds.AC<-exp(coef(m6.SoilNMR))
beta.prob.AC<-odds.AC/(1+odds.AC)

parameter.summary.AC<-data.frame(glm=m4.SoilNMR$coefficients,
                              beta=beta.prob.AC[-c(13)])

# save table
write.csv(parameter.summary.AC,"AC.parameters.csv",row.names = FALSE)

# ---- Fit statistics: SoilNMR -----

#Calculate r-squared for LR using Coz and Snell
library(rcompanion)
nagelkerke(m4.SoilNMR) #R^2 for GLM is the Cox and Snell R-squared value
nagelkerke(m6.SoilNMR) #R^2 for Beta is the Cox and Snell R-squared value

fit.SoilNMR<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                        glm=c(AIC(m4.SoilNMR),(1-0.9853)^(-36/2),
                              sum(m4.SoilNMR$residuals^2)/m4.SoilNMR$df.residual),
                        beta=c(AIC(m6.SoilNMR),(1-0.9830)^(-36/2),
                               sum(m6.SoilNMR$residuals^2)/m6.SoilNMR$df.residual))

# save table
write.csv(fit.SoilNMR,"Alkylcarbonfit.stats.csv",row.names = FALSE)



