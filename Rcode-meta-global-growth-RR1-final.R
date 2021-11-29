rm(list = ls()) #clear list

#automatic installation of required packages
packages <- c("xlsx","calibrate","stargazer","sandwich","lmtest","getopt","CausalGAM","ggplot2","reshape2","xts",
              "lattice","gridExtra","gtable","plm","lfe","lmtest","car","tis","foreign","MASS","quantreg","ggrepel",
              "dplyr","stringr","datasets","rio","psych","systemfit","MatchIt","CRTgeeDR","eurostat","plyr","zoo","ggthemes",
              "robumeta","metafor","dplyr","clubSandwich","Hmisc","metafor","pracma","pkgs","broom","sjPlot", "here", "data.table", "pscore", "AER")
ipak(packages)

#load packages
library(xlsx) #Excel-Paket laden
library(calibrate) #Laden des Pakets, das f??r Datenbeschriftung n??tig ist
library (stargazer) #Laden des Pakets, mit dem R-Regressionsoutput in Latex-Tabellen ??bergef??hrt werden kann
library(sandwich)
library(lmtest)
library(getopt)
library(CausalGAM)
library(ggplot2)
library(reshape2)
library(xts)
library(lattice)
library(gridExtra)
library(gtable)
library(plm)
library(lfe)
library(lmtest)
library(car)
library(tis)
library(foreign)
library(MASS)
library(quantreg)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(datasets)
library(rio)
library(psych)
library(systemfit)
library(foreign)
library(MatchIt)
library(CRTgeeDR)
library(eurostat)
library(plyr)
library(zoo)
library(ggthemes)
library("robumeta")
library("metafor")
library("dplyr")
library(clubSandwich)
library(Hmisc)
library(metafor)
library(pracma)
library(broom)
library(sjPlot)
library(here)
library(data.table)
library(pscore)
library(AER)

#Load data
dat <- fread(here("Coding-globalization-growth-RR.csv"))

#exclude estimates with interaction terms
dat <- subset(dat, InteractedGlobalization %in% c('0'))

#calculate the partial correlation coefficient
#formula: see Stanley and Doucouliagos (2012: p. 25)
dat$PartialCorrelationCoefficient <- dat$Tstatistic / (sqrt((dat$Tstatistic^2)+dat$DegreesofFreedom))

dat$PartialCorrelationCoefficient <- dat$PartialCorrelationCoefficient*dat$Transform
dat$Tstatistic <- dat$Tstatistic*dat$Transform

#standard error of the partial correlation coefficient
dat$StandardErrorPartialCorrelation <- sqrt((1-(dat$PartialCorrelationCoefficient)^2)/dat$DegreesofFreedom)

#Precision
dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation
#InverseSE
dat$InverseSE <- 1 / dat$StandardError
#Variance 
dat$Variance <- dat$StandardErrorPartialCorrelation^2
#PrecVariance
dat$PrecVariance <- 1 / dat$Variance

#Proxy for standard error
dat$proxySE <- 1 / (sqrt(dat$Observations))


#transform r to Z and calculate the corresponding sample variances
dat <- escalc(measure="ZCOR", ri=PartialCorrelationCoefficient, ni=Observations, data=dat)

#log transformaitons of variables
dat$YearofPublication <- log(dat$YearofPublication)
dat$NumberofCountries <- log(dat$NumberofCountries)

#transformations of variables
dat$StartYear <- as.numeric(as.integer(dat$StartYear))
dat$EndYear <- as.numeric(as.integer(dat$EndYear))

dat$LengthTimeSpan <- (dat$EndYear - dat$StartYear) + 1
dat$LengthTimeSpan <- log(dat$LengthTimeSpan)
dat$MeanYearData<- (dat$StartYear+dat$EndYear)/2
dat$MeanYearData <- log(dat$MeanYearData)
dat$InstrumentSE <- 1 / (sqrt(dat$DegreesofFreedom))
dat$InstrumentVariance <- dat$InstrumentSE^2
dat$PrecInstrumentVariance <- 1 / dat$InstrumentVariance
dat$Citations <- log(dat$Citations)
dat$PrecNumberRegressions <- 1 / dat$NumberRegressions
#normalize impact factor
dat$JournalImpactFactor <- as.numeric(dat$JournalImpactFactor)
dat$MaxImpactFactor <- max(dat$JournalImpactFactor)
dat$MaxImpactFactor <- as.numeric(dat$MaxImpactFactor)
dat$NormalizedImpactFactor <- dat$JournalImpactFactor / max(dat$JournalImpactFactor)

dat_long <- melt(dat, id=1:116)

#subsets of the data
dat_tradeglobal <- subset(dat_long, TradeGlobalization %in% c('1'))
dat_financialglobal <- subset(dat_long, FinancialGlobalization %in% c('1'))
dat_overallglobal <- subset(dat_long, OverallGlobalization %in% c('1'))

#Descriptive statistics (Table 1)
#all-set
#unweighted average
uwa <- mean(dat$PartialCorrelationCoefficient)
uwa

reguwa <- lm(PartialCorrelationCoefficient~1, data=dat)
summary(reguwa)

coef_test(reguwa, vcov = "CR0", 
          cluster = dat$id, test = "naive-t")

regwa <- lm(PartialCorrelationCoefficient~1, data=dat, weights=PrecVariance)
summary(regwa)

coef_test(regwa, vcov = "CR0", 
          cluster = dat$id, test = "naive-t")
confint(regwa)
#(precision-weighted) average
dat$PreSE <- 1/dat$Variance
numerator <- dat$PartialCorrelationCoefficient*dat$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(dat$PreSE, na.rm=TRUE)
wa

#median
median(dat$PartialCorrelationCoefficient, na.rm=TRUE)

#exclude top and bottom 10%
topbottom <- group_by(dat_long, id) %>%
  mutate(rank = rank(desc(PartialCorrelationCoefficient))) %>%
  filter(rank >=555) %>%
  filter(rank <=4987) %>%
  arrange(rank)

#unweighted average
uwa_topbottom <- sum(topbottom$PartialCorrelationCoefficient, na.rm=TRUE) / 4433
mean(topbottom$PartialCorrelationCoefficient)

#(precision-weighted) average
topbottom$PreSE <- 1/topbottom$Variance
numerator <- topbottom$PartialCorrelationCoefficient*topbottom$PreSE
wa_topbottom <- sum(numerator, na.rm=TRUE)/sum(topbottom$PreSE, na.rm=TRUE)
wa_topbottom

regwa <- lm(PartialCorrelationCoefficient~1, data=topbottom, weights=PrecVariance)
summary(regwa)

coef_test(regwa, vcov = "CR0", 
          cluster = topbottom$id, test = "naive-t")
confint(regwa)

#median
median(topbottom$PartialCorrelationCoefficient, na.rm=TRUE)

#trade globalization
#unweighted average
uwa_trade <- sum(dat_tradeglobal$PartialCorrelationCoefficient, na.rm=TRUE) / 3718
uwa_trade

#(precision-weighted) average
dat_tradeglobal$PreSE <- 1/dat_tradeglobal$Variance
numerator <- dat_tradeglobal$PartialCorrelationCoefficient*dat_tradeglobal$PreSE
wa_trade <- sum(numerator, na.rm=TRUE)/sum(dat_tradeglobal$PreSE, na.rm=TRUE)
wa_trade

regwa <- lm(PartialCorrelationCoefficient~1, data=dat_tradeglobal, weights=PrecVariance)
summary(regwa)

coef_test(regwa, vcov = "CR0", 
          cluster = dat_tradeglobal$id, test = "naive-t")
confint(regwa)

#median
median(dat_tradeglobal$PartialCorrelationCoefficient, na.rm=TRUE)

#financial globalization
#unweighted average
uwa_financial <- sum(dat_financialglobal$PartialCorrelationCoefficient, na.rm=TRUE) / 1788
uwa_financial

#(precision-weighted) average
dat_financialglobal$PreSE <- 1/dat_financialglobal$Variance
numerator <- dat_financialglobal$PartialCorrelationCoefficient*dat_financialglobal$PreSE
wa_financial <- sum(numerator, na.rm=TRUE)/sum(dat_financialglobal$PreSE, na.rm=TRUE)
wa_financial

regwa <- lm(PartialCorrelationCoefficient~1, data=dat_financialglobal, weights=PrecVariance)
summary(regwa)

coef_test(regwa, vcov = "CR0", 
          cluster = dat_financialglobal$id, test = "naive-t")
confint(regwa)

#median
median(dat_financialglobal$PartialCorrelationCoefficient, na.rm=TRUE)

#minimum, maximum, standard deviation
max(dat_long$PartialCorrelationCoefficient)
min(dat_long$PartialCorrelationCoefficient)
sd(dat_long$PartialCorrelationCoefficient)
mean(dat_long$PartialCorrelationCoefficient)
median(dat_long$PartialCorrelationCoefficient)

#write.csv(dat, "test_global.csv")

#Results on publication selection bias (table 3)

#all economic globalization variables
#column (1)
#FAT-PET test
pubbias_1 <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_long)
summary(pubbias_1)

coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (2)
#Fisher's z-transformed
pubbias_3 <- lm(yi ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_long)
summary(pubbias_3)

coef_test(pubbias_3, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (3)
#IV
library(AER)

pubbias_5 <- AER::ivreg(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation | InstrumentSE, weights=PrecInstrumentVariance, data=dat_long)
summary(pubbias_5)

coef_test(pubbias_5, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#trade globalization only
#column (4)
#FAT-PET test
pubbias_1_trade <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_tradeglobal)
summary(pubbias_1_trade)

coef_test(pubbias_1_trade, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_tradeglobal$id, test = "naive-t")

#financial globalization only
#column (5)
pubbias_1_financial <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_financialglobal)
summary(pubbias_1_financial)

coef_test(pubbias_1_financial, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_financialglobal$id, test = "naive-t")

#correct t-values and standard errors for stargazer table
#all globalization variable
ses_pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                  cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#trade globalization
ses_pubbias_1_trade <- list(coef_test(pubbias_1_trade, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                      cluster = dat_tradeglobal$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_trade <- list(coef_test(pubbias_1_trade, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                        cluster = dat_tradeglobal$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#financial globalization
ses_pubbias_1_financial <- list(coef_test(pubbias_1_financial, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                          cluster = dat_financialglobal$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1_financial <- list(coef_test(pubbias_1_financial, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                            cluster = dat_financialglobal$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#all globalization variable
ses_pubbias_3 <- list(coef_test(pubbias_3, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                          cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_3 <- list(coef_test(pubbias_3, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                            cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#all globalization variable
ses_pubbias_5 <- list(coef_test(pubbias_5, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_5 <- list(coef_test(pubbias_5, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                  cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#Visualization of publication selection bias
#Funnel plot (figure 1)

plot_funnel <- ggplot(data=dat,
                      aes(x=PartialCorrelationCoefficient, y=PrecSE)) +
  geom_point(aes(colour=factor(GlobalizationType)), size=0.5) +
  xlab("Partial correlation coefficient") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("Funnel plot of Globalization-growth\n partial correlations (n=5542)")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  geom_vline(xintercept=0, colour="black", linetype=2)+
  theme(legend.text = element_text(colour="black", size = 4))+
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))
plot_funnel

###
#Multivariate meta-regression analysis (table 4)

#column (1)
pubbias_var_est_gts <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + FinancialGlobalization + OverallGlobalization + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + LatinAmericaOnly + Prior + ReviewedJournal + Education + Institutions, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est_gts)

coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (2)
pubbias_var_est_gts_TopTier <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + FinancialGlobalization + OverallGlobalization + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + LatinAmericaOnly + Prior + TopTierJournal + Education + Institutions, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est_gts_TopTier)

coef_test(pubbias_var_est_gts_TopTier, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")


#column (3)
#IV estimation
library(AER)

pubbias_var_gts_IV <- AER::ivreg(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + FinancialGlobalization + OverallGlobalization + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + LatinAmericaOnly + Prior + ReviewedJournal + Education + Institutions | InstrumentSE + FinancialGlobalization + OverallGlobalization + LatinAmericaOnly + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + Prior + ReviewedJournal + Education + Institutions, weights=PrecVariance, data=dat_long)
summary(pubbias_var_gts_IV)

coef_test(pubbias_var_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (4)
#Trade globalization only - no consideration of different trade openness variables
pubbias_trade_var_est_gts <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + LatinAmericaOnly + Prior + ReviewedJournal + Education + Institutions, weights=PrecVariance, data=dat_tradeglobal)
summary(pubbias_trade_var_est_gts)

coef_test(pubbias_trade_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_tradeglobal$id, test = "naive-t")

#column (5)
#Financial globalization only - no consideration of different trade openness variables
pubbias_financial_var_est_gts <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + CrossSection + YearAverage + PerCapita + Dummy1870to1914 + Dummy1915to1945 + Dummy1990s + Dummy2000s + LatinAmericaOnly + Prior + ReviewedJournal + Education + Institutions, weights=PrecVariance, data=dat_financialglobal)
summary(pubbias_financial_var_est_gts)

coef_test(pubbias_financial_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_financialglobal$id, test = "naive-t")

#stargazer table
#WLS with standard errors clustered at the study level
#all globalization variable
ses.WLS.est <- list(coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est <- list(coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS.est
ses.WLS.est

ses.WLS.est_TopTier <- list(coef_test(pubbias_var_est_gts_TopTier, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est_TopTier <- list(coef_test(pubbias_var_est_gts_TopTier, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS.est_TopTier
ses.WLS.est_TopTier

#WLS (IV estimation) with standard errors clustered at the study level
ses.WLS.IV <- list(coef_test(pubbias_var_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                 cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.IV <- list(coef_test(pubbias_var_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                   cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#trade globalization
ses.WLS.est_trade <- list(coef_test(pubbias_trade_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                    cluster = dat_tradeglobal$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est_trade <- list(coef_test(pubbias_trade_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                      cluster = dat_tradeglobal$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#financial globalization
ses.WLS.est_financial <- list(coef_test(pubbias_financial_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                        cluster = dat_financialglobal$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est_financial <- list(coef_test(pubbias_financial_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                          cluster = dat_financialglobal$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#Regression results reported in table 4; IV not included because of stargazer compatibility issues
stargazer(pubbias_var_est_gts, pubbias_var_est_gts_TopTier, pubbias_trade_var_est_gts, pubbias_financial_var_est_gts, t=list(unlist(tvals.WLS.est), unlist(tvals.WLS.est_TopTier), unlist(tvals.WLS.est_trade), unlist(tvals.WLS.est_financial)), se=list(unlist(ses.WLS.est), unlist(ses.WLS.est_TopTier), unlist(ses.WLS.est_trade), unlist(ses.WLS.est_financial)))

#Descriptive statistics (see Table 2)
#mean
mean(dat_long$PartialCorrelationCoefficient)

mean(dat_long$TradeGlobalization)
mean(dat_long$FinancialGlobalization)
mean(dat_long$OverallGlobalization)

mean(dat_long$CrossSection)
mean(dat_long$YearAverage)
mean(dat_long$IntraNationalStudy)
mean(dat_long$PerCapita)
mean(dat_long$Dummy1870to1914)
mean(dat_long$Dummy1915to1945)
mean(dat_long$Dummy1946to1969)
mean(dat_long$Dummy1970s)
mean(dat_long$Dummy1980s)
mean(dat_long$Dummy1990s)
mean(dat_long$Dummy2000s)
mean(dat_long$Dummy2010s)

mean(dat_long$EuropeOnly)
mean(dat_long$NorthAmericaOnly)
mean(dat_long$AfricaOnly)
mean(dat_long$AsiaOnly)
mean(dat_long$LatinAmericaOnly)
mean(dat_long$ChinaOnly)

mean(dat_long$CountryFixedEffects)
mean(dat_long$TimeFixedEffects)
mean(dat_long$NonOLS)
mean(dat_long$LongRun)
mean(dat_long$TacklingEndogeneity)

mean(dat_long$StandardErrorPartialCorrelation)
mean(dat_long$ReviewedJournal)
mean(dat_long$Primary)
mean(dat_long$Prior)
mean(dat_long$Crossauthor)
mean(dat_long$YearofPublication)
mean(dat_long$Citations)
mean(dat_long$NormalizedImpactFactor)

mean(dat_long$Education)
mean(dat_long$LaborForce)
mean(dat_long$GovernmentSpending)
mean(dat_long$Geography)
mean(dat_long$FinancialDevelopment)
mean(dat_long$Institutions)
mean(dat_long$NaturalResources)
mean(dat_long$Democracy)
mean(dat_long$CapitalStock)
mean(dat_long$TopTierJournal)


#standard deviation
sd(dat_long$PartialCorrelationCoefficient)

sd(dat_long$TradeGlobalization)
sd(dat_long$FinancialGlobalization)
sd(dat_long$OverallGlobalization)

sd(dat_long$CrossSection)
sd(dat_long$YearAverage)
sd(dat_long$IntraNationalStudy)
sd(dat_long$PerCapita)
sd(dat_long$Dummy1870to1914)
sd(dat_long$Dummy1915to1945)
sd(dat_long$Dummy1946to1969)
sd(dat_long$Dummy1970s)
sd(dat_long$Dummy1980s)
sd(dat_long$Dummy1990s)
sd(dat_long$Dummy2000s)
sd(dat_long$Dummy2010s)

sd(dat_long$EuropeOnly)
sd(dat_long$NorthAmericaOnly)
sd(dat_long$AfricaOnly)
sd(dat_long$AsiaOnly)
sd(dat_long$LatinAmericaOnly)
sd(dat_long$ChinaOnly)

sd(dat_long$CountryFixedEffects)
sd(dat_long$TimeFixedEffects)
sd(dat_long$NonOLS)
sd(dat_long$LongRun)
sd(dat_long$TacklingEndogeneity)

sd(dat_long$StandardErrorPartialCorrelation)
sd(dat_long$ReviewedJournal)
sd(dat_long$Primary)
sd(dat_long$Prior)
sd(dat_long$Crossauthor)
sd(dat_long$YearofPublication)
sd(dat_long$Citations)
sd(dat_long$NormalizedImpactFactor)

sd(dat_long$Education)
sd(dat_long$LaborForce)
sd(dat_long$GovernmentSpending)
sd(dat_long$Geography)
sd(dat_long$FinancialDevelopment)
sd(dat_long$Institutions)
sd(dat_long$NaturalResources)
sd(dat_long$Democracy)
sd(dat_long$CapitalStock)
sd(dat_long$TopTierJournal)

#Ioannidis et al.2017 (see table A3)
data_new2 <- dat %>%                                      # Top N highest values by group
  arrange(desc(StandardErrorPartialCorrelation)) %>% 
  slice(4987:5542)

regwa_Ioannidis <- lm(PartialCorrelationCoefficient~1, data=data_new2, weights=PrecVariance)
summary(regwa_Ioannidis)
coef_test(regwa_Ioannidis, vcov = "CR0", 
          cluster = data_new2$id, test = "naive-t")
