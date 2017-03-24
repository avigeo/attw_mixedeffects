# attw_mixedeffects
code used to model ATTW abundance as a function of GIS covariates


setwd("~/Documents/RData/mixedeffects_attwgis")
library(maptools); library(foreign); library(lattice); library(spgwr); library(gam); library(pscl); library(maps); library(splancs); library(spdep); library(graphics); library(MASS); library(grid); library(ggplot2); library(lme4); library(visreg); library(sjPlot); library(sjmisc); library("arm"); library(RColorBrewer); require(gridExtra)



#--- Read Data for Southern Rockies ---#
SRdata <- read.csv("srdata.csv",sep=",",header=T)
names(SRdata)
NRdata <- read.csv("nrdata1.csv",sep=",",header=T)
names(NRdata)

#--- z-score all variables ---#
SRdata.z <- SRdata
NRdata.z <- NRdata

rows <- c("elev","near","per_sf","per_lp","per_mc","per_as","per_pp","per_o")
cols <- c("mean","SD")
varsum <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))

for(i in 1:length(rows)) {
  varsum[i,"mean"] <- mean(c(SRdata[,rows[i]],NRdata[,rows[i]]))
  varsum[i,"SD"] <- sd(c(SRdata[,rows[i]],NRdata[,rows[i]]))
  SRdata.z[,rows[i]] <- (SRdata.z[,rows[i]] - varsum[i,"mean"])/varsum[i,"SD"]
  NRdata.z[,rows[i]] <- (NRdata.z[,rows[i]] - varsum[i,"mean"])/varsum[i,"SD"]
}
head(SRdata.z)
#---generalized poisson mixed effects model with year as a random effect because I lumped all of the data together--#
# Null model, intercept only model
sr0 <- glmer(count1 ~ (1|year), data=SRdata.z, family="poisson")
summary(sr0)
sr1 <- glmer(count1 ~ elev + near + per_sf + per_lp + per_mc + per_as + per_pp + (1|year), data=SRdata.z, family="poisson")
summary(sr1)
#the one below (sr2) is the most supported model!
sr2 <- glmer(count1 ~ elev + near + per_sf + (1|year), data=SRdata.z, family="poisson")
summary(sr2)
sr3 <- glmer(count1 ~ elev + near + per_mc + (1|year), data=SRdata.z, family="poisson")
summary(sr3)
sr4 <- glmer(count1 ~ near + per_sf + (1|year), data=SRdata.z, family="poisson")
summary(sr4)
sr5 <- glmer(count1 ~ elev + near + (1|year), data=SRdata.z, family="poisson")
summary(sr5)
sr6 <- glmer(count1 ~ elev + near + per_sf + per_mc + per_as + (1|year), data=SRdata.z, family="poisson")
summary(sr6)
sr7 <- glmer(count1 ~ elev + near + per_sf + per_mc + per_as + per_pp + (1|year), data=SRdata.z, family="poisson")
summary(sr7)
sr8 <- glmer(count1 ~ elev + per_sf + per_lp + per_mc + per_pp + (1|year), data=SRdata.z, family="poisson")
summary(sr8)


#---BIC table---#
rows <- seq(0,8)
cols <- c("Model","BIC","deltaBIC","no fixed","wtno","wt")
BICout <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
BICout[,"Model"] <- c("Year effect only","elev + near + per_sf + per_lp + per_mc + per_as + per_pp",
                      "elev + near + per_sf","elev + near + per_mc","near + per_sf","elev + near",
                      "elev + near + per_sf + per_mc + per_as","elev + near + per_sf + per_mc + per_as + per_pp", "elev + per_sf + per_lp + per_mc + per_pp")
BICout[,"no fixed"] <- c(0,7,3,3,2,2,5,6,5)
anv <- anova(sr0, sr1, sr2, sr3, sr4, sr5, sr6, sr7, sr8)
BICout[,"BIC"] <- anv$BIC[order(dimnames(anv)[[1]])]
BICout <- data.frame(BICout,stringsAsFactors=F)
for(i in 2:6) BICout[,i] <- as.numeric(BICout[,i])

BICout <- BICout[order(BICout$BIC),] # Orders by BIC
BICout$deltaBIC <- BICout$BIC - min(BICout$BIC)

BICout$wtno <- exp(-0.5*BICout$deltaBIC)
BICout$wt <- round(BICout$wtno/sum(BICout$wtno),digits=3)
BICout <- BICout[,-which(names(BICout)=="wtno")]
write.csv(BICout,"Model_selection_SR2.csv",row.names=F)

#--- Estimates for selected model ---#
select.mod <- sr2
coef <- summary(select.mod)$coefficients
coef <- data.frame(coef)
coef$lower <- coef$Estimate - coef$Std..Error*1.96
coef$upper <- coef$Estimate + coef$Std..Error*1.96
summary(select.mod) # Get random effect from here ("Random effects:")
write.csv(coef,"Coefficients_SR2.csv",row.names=T)

#view estimated fixed effect coefficients
fixef(sr2)
#view predicted random effects
ranef(sr2)
#view coefficients for LMM for each group
coef(sr2)
#compute confidence intervals on the parameters (cutoffs based on liklihood ratio test
confint(sr2)

#---------#----------Northern Rockies-----------#------------#

#first check null model of intercept only model with random effect of year
nr0 <- glmer(count1 ~ (1|year), data=NRdata.z, family="poisson")
summary(nr0)
nr1 <- glmer(count1 ~ elev + near + per_sf + per_lp + per_mc + per_as + per_pp + (1|year), data=NRdata.z, family="poisson")
summary(nr1)
nr2 <- glmer(count1 ~ elev + near + per_sf + (1|year), data=NRdata.z, family="poisson")
summary(nr2)
nr3 <- glmer(count1 ~ elev + near + per_mc + (1|year), data=NRdata.z, family="poisson")
summary(nr3)
nr4 <- glmer(count1 ~ near + per_sf + (1|year), data=NRdata.z, family="poisson")
summary(nr4)
nr5 <- glmer(count1 ~ elev + near + (1|year), data=NRdata.z, family="poisson")
summary(nr5)
nr6 <- glmer(count1 ~ elev + near + per_sf + per_mc + per_as + (1|year), data=NRdata.z, family="poisson")
summary(nr6)
#--the one below (nr8) is the most supported model in the NR--#
nr8 <- glmer(count1 ~ elev + per_sf + per_lp + per_mc + per_pp + (1|year), data=NRdata.z, family="poisson")
summary(nr8)


#---BIC table---#
rows <- seq(0,7)
cols <- c("Model","BIC","deltaBIC","no fixed","wtno","wt")
BICout <- matrix(NA,nrow=length(rows),ncol=length(cols),dimnames=list(rows,cols))
BICout[,"Model"] <- c("Year effect only","elev + near + per_sf + per_lp + per_mc + per_as + per_pp",
                      "elev + near + per_sf","elev + near + per_mc","near + per_sf","elev + near",
                      "elev + near + per_sf + per_mc + per_as", "elev + per_sf + per_lp + per_mc + per_pp")
BICout[,"no fixed"] <- c(0,7,3,3,2,2,5,5)
anv <- anova(nr0, nr1, nr2, nr3, nr4, nr5, nr6, nr8)
BICout[,"BIC"] <- anv$BIC[order(dimnames(anv)[[1]])]
BICout <- data.frame(BICout,stringsAsFactors=F)
for(i in 2:6) BICout[,i] <- as.numeric(BICout[,i])

BICout <- BICout[order(BICout$BIC),] # Orders by BIC
BICout$deltaBIC <- BICout$BIC - min(BICout$BIC)

BICout$wtno <- exp(-0.5*BICout$deltaBIC)
BICout$wt <- round(BICout$wtno/sum(BICout$wtno),digits=3)
BICout <- BICout[,-which(names(BICout)=="wtno")]
write.csv(BICout,"Model_selection_NR3.csv",row.names=F)

#--- Estimates for selected model CHOOSE IT!---#
select.modN <- nr8
coefN <- summary(select.modN)$coefficients
coefN <- data.frame(coefN)
coefN$lower <- coefN$Estimate - coefN$Std..Error*1.96
coefN$upper <- coefN$Estimate + coefN$Std..Error*1.96
summary(select.modN) # Get random effect from here ("Random effects:")
write.csv(coefN,"Coefficients_NR8_2.csv",row.names=T)


tapply(NRdata$count1, NRdata$year, sd)
tapply(SRdata$count1, SRdata$year, mean)

#Feb 9, 2017-- Trying more plotting for PLOS ONE pub ---# #load ggthemes
library(sjPlot); library(sjmisc); library(ggplot2); require(ggthemes)
summary(sr2)
SRdata.z$fit <- predict(sr2)

sjp.setTheme(geom.label.size = 3)

#before plotting, need to change variable names so that they show up on the plot because the plotting function does not support changing axis labesl.
colnames(SRdata.z)[colnames(SRdata.z)=="count1"] <- "Counts"
colnames(SRdata.z)[colnames(SRdata.z)=="elev"] <- "Elevation"
colnames(SRdata.z)[colnames(SRdata.z)=="near"] <- "Distance"
colnames(SRdata.z)[colnames(SRdata.z)=="per_sf"] <- "SpruceFir"
colnames(NRdata.z)[colnames(NRdata.z)=="count1"] <- "Counts"
colnames(NRdata.z)[colnames(NRdata.z)=="per_sf"] <- "SpruceFir"
colnames(NRdata.z)[colnames(NRdata.z)=="per_mc"] <- "MixedConifer"
colnames(NRdata.z)[colnames(NRdata.z)=="per_pp"] <- "Ponderosa"
colnames(NRdata.z)[colnames(NRdata.z)=="per_lp"] <- "Lodgepole"
colnames(NRdata.z)[colnames(NRdata.z)=="elev"] <- "Elevation"

#rerun model with new names for graphing
s2 <- glmer(Counts ~ Elevation + Distance + SpruceFir + (1|year), data=SRdata.z, family="poisson")
summary(s2)

n8 <- glmer(Counts ~ Elevation + SpruceFir + Lodgepole + MixedConifer + Ponderosa + (1|year), data=NRdata.z, family="poisson")
summary(n8) 

#Marginal effects plots of predicted log counts for each fixed term, where remaining co-variates are set to the mean. Use facet.grid to decide whether to plot each coefficient as separate plot or as integrated faceted plot. 
# THESE ARE FIG 4 IN THE PUB #
Fig4 <- sjp.glmer(s2, 
               type = "eff", 
               vars = c("Distance", "Elevation", "SpruceFir"),
               sort.est = "sort.all",
               show.ci = TRUE, 
               title = " ", 
               facet.grid = TRUE,
               free.scale = FALSE,
               hjust = "center"
               )$plot + theme_base() + labs(x="Z-score of Fixed-Effect Parameter", y="Log of Predicted Counts")
quartz()
Fig4
save_plot("Fig4.tif", fig = ggplot2::last_plot(), width = 10, height = 11,
          dpi = 300, label.size = 2.4,
          axis.textsize = 0.8, axis.titlesize = 0.75)

Fig3 <- sjp.glmer(n8, 
                type = "eff", 
                vars = c("Elevation", "SpruceFir", "Lodgepole", "MixedConifer", "Ponderosa"),
                sort.est = "sort.all",
                show.ci = TRUE, 
                title=" ",
                facet.grid = TRUE,
                free.scale = FALSE,
                hjust = "center"
                )$plot + theme_base() + labs(x="Z-score of Fixed-Effect Parameter", y="Log of Predicted Counts")

quartz()
Fig3
save_plot("Fig3.tif", fig = ggplot2::last_plot(), width = 9, height = 12,
          dpi = 300, label.size = 2,
          axis.textsize = 0.8, axis.titlesize = 0.75)
