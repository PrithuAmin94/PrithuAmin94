################################################################
library(pxweb)
#the price index of the swedish housing price.
monthly_prices <- 
  get_pxweb_data(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/BO/BO0501/BO0501B/FastprisManadAr",
                 dims = list(Forvarvsmanad = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'),
                             Fastighetstyp = c('1'),
                             ContentsCode = c('BO0501I7'),
                             Tid = c('*')),
                 clean = TRUE)
values <- monthly_prices[, "values"]
monthly_prices1 <- ts(values, start = c(1992, 1),
                      end = c(2019, 12), frequency = 12) 
quaterly_prices <- aggregate(monthly_prices1, FUN = mean, nfrequency = 4)
quaterly_prices
###############################################################
#inflation:
CPIF_monthly <- 
  get_pxweb_data(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/PR/PR0101/PR0101G/KPIFInd",
                 dims = list(ContentsCode = c('000001R4'),
                             Tid = c('*')),
                 clean = TRUE)
CPIF_monthly_value <- CPIF_monthly[, "values"] 
CPIF <- ts(CPIF_monthly_value, start = c(1992, 1), 
           end = c(2019, 12), frequency = 12)
CPIF <- aggregate(CPIF, FUN = mean, nfrequency = 4)
CPIF

#deflate:
housing_prices <- (quaterly_prices / CPIF) * 100
plot.ts(housing_prices, lty = 4, col = "blue", main = "Housing Prices")
###############################################################
library(urca)
#Function for adf:
ADF.test <- function(data) {
  ADF <- function(type, data) {
    result1 <- ur.df(data, type = type, lags = 3 * frequency(data), 
                     selectlags = "AIC")
    DETERM <- ifelse(type == "trend", 2, ifelse(type == "drift", 
                                                1, 0))
    LAGS <- length(coefficients(result1@testreg)[, 1]) - 
      DETERM - 1
    result2 <- cbind(t(result1@teststat), result1@cval, coefficients(result1@testreg)["z.lag.1", 
                                                                                      1], LAGS)
    round(result2, 2)
  }
  types <- c("trend", "drift", "none")
  result3 <- apply(t(types), 2, ADF, data)
  cat(rep("#", 20), "\n")
  cat(rep("#", 20), "\n")
  cat("Augmented Dickey--Fuller test\n")
  cat(rep("#", 20), "\n")
  cat("type:", "  trend ", "drift ", "none\n")
  cat("AR1:   ", result3[[1]][1, 5] + 1, " ", result3[[2]][1, 
                                                           5] + 1, " ", result3[[3]][5] + 1, "\n")
  cat("lags:  ", result3[[1]][1, 6], "   ", result3[[2]][1, 
                                                         6], "   ", result3[[3]][6], "\n")
  cat(rep("#", 20), "\n")
  result5 <- rbind(result3[[1]][c(1, 3), 1:4], result3[[2]][1:2, 
                                                            1:4], result3[[3]][1:4])
  rownames(result5)[5] <- "tau1"
  result5
}
#################################################################
ADF.test(housing_prices)
#tau3 > cval --> Not stationary
library(urca)
ur.kpss(housing_prices, type = "tau")@teststat 
ur.kpss(housing_prices, type = "tau")@cval
ur.kpss(housing_prices, type = "mu")@teststat
ur.kpss(housing_prices, type = "mu")@cval
#stat > cval --> Not stationary
###############################################################
#needs differencing:
housing_prices_diff <- diff(housing_prices)
ADF.test(housing_prices_diff)
ur.kpss(housing_prices_diff, type = "tau")@teststat 
ur.kpss(housing_prices_diff, type = "tau")@cval
ur.kpss(housing_prices_diff, type = "mu")@teststat
ur.kpss(housing_prices_diff, type = "mu")@cval
asd <- diff(housing_prices_diff)
ADF.test(asd)
ur.kpss(asd, type = "tau")@teststat 
ur.kpss(asd, type = "tau")@cval
ur.kpss(asd, type = "mu")@teststat
ur.kpss(asd, type = "mu")@cval
#it was integrated of order 2 --> I(2)
###############################################################
par(mfrow = c(2,1))
acf(asd, lag.max = 80)
pacf(asd, lag.max = 80)
#shows that i have decaying acf and Pacf of lag order 3
#so this will be an AR(3) model
###############################################################
library(forecast)
auto.arima(asd, ic = "bic", trace = TRUE)
#auto.arima suggested ARIMA(3,0,0)(0,1,2)[4]
###############################################################
par(mfrow = c(1,1))
#deseasonalising:
library(seasonal)
HousePrice <- seas(asd)
plot(HousePrice)
summary(HousePrice)
HousePrice <- HousePrice$data[,"seasonaladj"]
###############################################################
auto.arima(HousePrice, ic = "bic", trace = TRUE)
#now auto.arima suggests ARIMA(0,0,1) with zero mean
###############################################################
par(mfrow = c(2,1))
acf(HousePrice)
pacf(HousePrice)
#decaying pacf and acf of lag order 1. 
#so this will be a ARIMA(0,0,1)
###############################################################
library(portes)
library(tseries)
result.001 <- arima(HousePrice, order = c(0,0,1))
LjungBox(result.001, lags = 1:10)
#suggests autocorrelation in higher lag orders.
#so need to include AR.
result.401<- arima(HousePrice, 
                        order = c(4,0,1))
LjungBox(result.401, lags = 1:30)
result.501<- arima(HousePrice, 
                        order = c(5,0,1))
LjungBox(result.501, lags = 1:30)
result.601<- arima(HousePrice, 
                        order = c(6,0,1))
LjungBox(result.601, lags = 1:30)
result.701 <- arima(HousePrice,
                    order = c(7,0,1))
LjungBox(result.701, lags = 1:30)
BIC(result.401); BIC(result.501)
BIC(result.601); BIC(result.701)
#we choose ARIMA(6,0,1)
################################################################
library(tseries)
jarque.bera.test(resid(result.601))
#we can reject null hypothesis
#our residuals are not normal
################################################################
par(mfrow = c(1,1))
qqnorm(resid(result.601))
qqline(resid(result.601))
#QQ-plot follows with our tests. We see extreme values.
################################################################
#re-estimating model for extreme values:
RES.601 <- residuals(result.601)
MAX.601 <- RES.601 == max(RES.601)
MIN.601 <- RES.601 == min(RES.601)
MAXMIN.601 <- cbind(MAX.601, MIN.601)
################################################################
result.601X <- arima(HousePrice, order = c(6,0,1),
                     xreg = MAXMIN.601)
BIC(result.601); BIC(result.601X)
#even though it seems like our BIC chooses arima(6,0,1),
#we decide to go with result.601X, as it has normal residuals.
################################################################
jarque.bera.test(resid(result.601X))
#residuals are normal.
################################################################
################################################################
# Variable 2:
GDP <- 
  get_pxweb_data(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/NR/NR0103/NR0103B/NR0103ENS2010T10SKv",
                 dims = list(Anvandningstyp = c('BNPM'),
                             ContentsCode = c('NR0103CF'),
                             Tid = c('*')),
                 clean = TRUE)
values <- GDP[, "values"]
GDP1 <- ts(values, start= c(1980, 1), 
           frequency = 4) 
GDP1 <- window(GDP1, start=c(1992, 1), 
               end=c(2019,4))
################################################################
ADF.test(GDP1)#tau 3< cval: Stationary
ur.kpss(GDP1, type = "tau")@teststat #stat < cval: Stationary
ur.kpss(GDP1, type = "tau")@cval
ur.kpss(GDP1, type = "mu")@teststat #stat < cval: Stationary
ur.kpss(GDP1, type = "mu")@cval
#Cointegration not present as GDP1 is integrated to order 0. I(0)
################################################################
library(vars)
var.data <- ts.intersect(HousePrice, GDP1)
colnames(var.data) <- c("HP", "GDP")
ADF.test(var.data[, "HP"])#tau 3< cval: Stationary, trend and drift
ADF.test(var.data[, "GDP"])#tau 3< cval: Stationary, trend and drift
################################################################
VARselect(var.data, type = "both")$selection #BIC suggests 2 lags,
#we chose type = both as we have trend and drift.
var.result <- VAR(var.data, p = 2, type= "both")
################################################################
library(lmtest)
coeftest(var.result$varresult$HP)#to check if the coefficients are statistically significant.
coeftest(var.result$varresult$GDP)#to check if the coefficients are statistically significant.
################################################################
var.resid <- residuals(var.result)
cov(var.resid)
serial.test(var.result)# P-Value is large: No autocorrelation
normality.test(var.result)$jb.mul$JB# small P-Value: Non normal residuals
normality.test(var.result)$jb.mul$Skewness# small P- Value: Skewed.
normality.test(var.result)$jb.mul$kurtosis# Cannot generate value for this command.
#since kurtosis cannot be achieved using the commands above,
#I decided to do a normality.test() without further arguments:
normality.test(var.result)
################################################################
qqnorm(var.resid)#we can see from the qqnorm plot that the extreme values are so far out, 
qqline(var.resid)#that taking extreme values as an external regressor might not help.
################################################################
HP.max <- var.resid[, "HP"] == max(var.resid[, "HP"])
HP.min <- var.resid[, "HP"] == min(var.resid[, "HP"])
GDP.max <- var.resid[, "GDP"] == max(var.resid[, "GDP"])
GDP.min <- var.resid[, "GDP"] == min(var.resid[, "GDP"])
var.data.ext <- cbind(HP.max, HP.min, 
                      GDP.max, GDP.min)
################################################################
var.data2 <- window(var.data, start = c(1993, 1))#we had to decrease the size of the data
#as the row size for the var.data had to match the row size for the dummy variables taken 
#above.
VARselect(var.data2, type = "both", exogen = var.data.ext)$selection
#this time BIC suggests, p=2. 
var.result2 <- VAR(var.data2, p= 2, type = "both", exogen = var.data.ext)
################################################################
qqnorm(resid(var.result2))#we can see that there is a slight improvement in the qqplot,
qqline(resid(var.result2))#however it does not fix it completely. 
################################################################
var.resid2 <- residuals(var.result2)
cov(var.resid2)
serial.test(var.result2)# P-Value is low: Autocorrelation is present 
normality.test(var.result2, multivariate.only = FALSE)$jb.mul$JB
# Large P-Value: normal residuals
normality.test(var.result2, multivariate.only = FALSE)$jb.mul$Skewness
# Large P- Value: Not Skewed
normality.test(var.result2, multivariate.only = FALSE)$jb.mul$kurtosis
# Cannot generate value for this command.
#since kurtosis cannot be achieved using the commands above,
#I decided to do a normality.test() without further arguments:
normality.test(var.result2)
################################################################
#Eigen Values:
roots(var.result2)
#roots are identical and some are unique
#roots are within unit interval. 
#var has stability
################################################################
################################################################
################################################################
# POINT FORECAST:
################################################################
# Point forecast for my arima model, for the next two years -->
library(forecast)
arima.res <- auto.arima(var.data2[, "HP"])
forecast(arima.res)
predict(arima.res, n.ahead = 4)
MAXMIN.601_2 <- matrix(c(rep(0,110), rep(1, 110)), ncol = 2, nrow = 8)
# here we are creating a matrix of 8 by 2, to fit the external regressors.
forecastAR <- predict(result.601X, n.ahead = 8, ci = 0.95, newxreg = MAXMIN.601_2)
forecastAR$pred
# Point forecast for my VAR model, for the next two years -->
var.data.ext2 <- matrix(c(rep(0, 108), rep(1, 108)), ncol = 4, nrow = 8)
# here we are creating a matrix of 8 by 4, to fit the external regressors.
forecastVAR <- forecast(var.result2, h = 8, dumvar = var.data.ext2, fan = TRUE)
forecastVAR$forecast$HP
plot(forecastVAR$forecast$HP)
################################################################
# PSEUDO OUT OF SAMPLE FORECASTS:
################################################################
# Ideally, We need real life data for rolling and recursive Pseudo Out of Sample Forecasts.
# since I had data until 2019, I have instead used that.I have forecasted data for 2018 and 
#2019 and compared with the real data. Also, I have chosen to use ARIMA(6,0,1). Auto arima 
#fails to forecast for 2018 and 2019.This is because the auto.arima suggests a MA(1) model, 
#which means that the shocks die down and come down to zero. 

#Recursive Arima:
ARIMA.RECURSIVE <- ts(matrix(NA,           
                       nrow = 5,     #the number of rows were calculated according to the  
                       ncol = 4),    #formula provided by Andre'
                start = c(2018, 1),
                frequency = 4)
colnames(ARIMA.RECURSIVE) <- paste("Horizon", 1:4, sep = " ")

for (i in 1:5){                                             
  EndDate       <- 2017.75 + (i - 1) / 4                    
  Data          <- window(HousePrice,
                          end = EndDate)                  
  Result        <-  arima(Data, order = c(6,0,1))       
  ARIMA.RECURSIVE[i, ] <- forecast(Result, h = 4)$mean            
}   
ARIMA.RECURSIVE   
#It would be better if I could fit external regressors here, as the model with external 
#regressors have normal residuals which the Arima(6,0,1) does not have.
                                        
#Rolling Arima:
ARIMA.ROLLING  <- ts(matrix(NA, nrow = 5, ncol = 4),
                     start = c(2018, 1), frequency = 4)
colnames(ARIMA.ROLLING) <- paste("Horizon", 1:4, sep = " ")

for (i in 1:5){
  StartDate     <- 1992 + 2/4 + (i - 1) / 4
  EndDate       <- 2017.75 + (i - 1) / 4
  Data          <- window(HousePrice,
                          start = StartDate,
                          end = EndDate)
  Result        <- arima(Data, order = c(6,0,1))
  ARIMA.ROLLING[i, ] <- forecast(Result, h = 4)$mean
}
ARIMA.ROLLING
##############################################################
#Recursive VAR
VAR.RESS <- ts(matrix(NA, 
                      nrow = 5, 
                      ncol = 4),
               start = c(2018, 1), 
               frequency = 4)
colnames(VAR.RESS) <- paste("Horizon", 1:4, sep = " ")
for(i in 1:5){
  EndDate   <- 2017.75 + (i - 1)/ 4
  Data      <- window(var.data2,  end = EndDate)
  Result    <- VAR(Data, p= 2, type = "both")
  VAR.RESS[i,] <- predict(Result, n.ahead =4)$fcst$HP[,"fcst"]
}
VAR.RESS
#It would be better if I could fit external regressors hereas well, since the model external 
#regressors have normal residuals which the this var model does not have.

#ROlling VAR
VAR.RO <- ts(matrix(NA, nrow = 5, ncol = 4),
             start= c(2018, 1), 
             frequency = 4)
colnames(VAR.RO) <- paste("Horizon", 1:4, sep  = " ")

for (i in 1:5){
  StartDate <- 1993 + 1/4 + (i -1) / 4
  EndDate <- 2017.75 + (i -1) / 4
  Data <- window(var.data2, start = StartDate, end = EndDate)
  Result <- VAR(Data, p = 2, type = "both")
  VAR.RO[i, ] <- predict(Result, n.ahead =4)$fcst$HP[, "fcst"]
}
VAR.RO
##########################################################

#Outcome Model:
OUTCOME  <- ts(matrix(NA, nrow = 5, ncol = 4),
               start = c(2018, 1), frequency = 4)

for (i in 1:5) {
  StartDate    <- 2018 + (i - 1) / 4         
  EndDate      <- StartDate + 0.75           
  OUTCOME[i, ] <- window(HousePrice,  
                         start = StartDate,
                         end = EndDate)
}

OUTCOME

#######################################################
#FORECAST EVALUATION:
#######################################################

#Arima Forecast Errors:

ARIMA.RES.FE <- ARIMA.RECURSIVE - OUTCOME
ARIMA.RES.FE
ARIMA.RES.M <- colMeans(ARIMA.RES.FE)
ARIMA.RES.M

ARIMA.RO.FE <- ARIMA.ROLLING - OUTCOME
ARIMA.RO.FE
ARIMA.RO.M  <- colMeans(ARIMA.RO.FE)
ARIMA.RO.M

#Var Forecast Errors:

VAR.RESS.FE <- VAR.RESS - OUTCOME
VAR.RESS.FE
VAR.RESS.M  <- colMeans(VAR.RESS.FE)
VAR.RESS.M

VAR.RO.FE <- VAR.RO - OUTCOME
VAR.RO.FE
VAR.RO.M  <- colMeans(VAR.RO.FE)
VAR.RO.M

#Table of Average Forecasting errors from different model:

Errors <- matrix(c(ARIMA.RES.M, ARIMA.RO.M, VAR.RESS.M, VAR.RO.M), ncol = 4, byrow = TRUE) 
colnames(Errors) <- c("1-step", "2-steps", "3-steps", "4-steps")
rownames(Errors) <- c("ARIMA.RE","ARIMA.RO","VAR.RE", "VAR.RO")
Erorrs <- as.table(Errors)
Errors

rowMeans(Errors)
#we can see that our Arima model is performing far better than our VAR model.this is becasuse
#our VAR model has a lot of extreme values, which could not be adjusted for. 

########################################################
#Function for Bias test:
bias.test <- function(h, FE) {
  require(dynlm); require(lmtest); require(sandwich)
  model <- dynlm(FE[, h] ~ 1) #here we are regressing our forecasts against the intercept. 
  matrix <- NeweyWest(model, lag = h - 1)
  round(coeftest(model, vcov. = matrix)[[4]], 2)
}
#This function takes two arguments, h --> Horizons, FE --> Forecasting Errors
########################################################
#Arima Bias Test:
# Here our command will generate P values for each horizons, For each Model.

ARIMA.RE.p <- apply(t(1:4), 2, bias.test, ARIMA.RES.FE)
names(ARIMA.RE.p) <- 1:4
ARIMA.RE.p
#All the values are higher than 0.05, Null is not rejected. Null Hypothesis: No Bias.

ARIMA.RO.p <- apply(t(1:4), 2, bias.test, ARIMA.RO.FE)
names(ARIMA.RO.p) <- 1:4
ARIMA.RO.p
#All the values are higher than 0.05, Null is not rejected. Null Hypothesis: No Bias.

#VAR Bias Test:

VAR.RE.p <- apply(t(1:4), 2, bias.test, VAR.RESS.FE)
names(VAR.RE.p) <- 1:4
VAR.RE.p
#All the values are higher than 0.05, Null is not rejected. Null Hypothesis: No Bias.

VAR.RO.p <- apply(t(1:4), 2, bias.test, VAR.RO.FE)
names(VAR.RO.p) <- 1:4
VAR.RO.p
#All the values are higher than 0.05, Null is not rejected. Null Hypothesis: No Bias.
########################################################
#FORECAST PRECISION TEST:
########################################################
#we will first create a function of the Diebold - Mariano test, which regresses this LD with
#the intercept.
DM.TEST <- function(ld, h) {
  require(dynlm); require(lmtest); require(sandwich)
  fe1 <- get(ld[1]); fe2 <- get(ld[2])
  res <- dynlm(I(fe1[, h] ^2- fe2[, h]^2) ~1)
  mat <- NeweyWest(res, lag = h - 1)
  round(coeftest(res, vcov. = mat)[4],2)
}
########################################################
#we will now test the hypothesis of equal fore cast precison with the Diebold - Mariano Test.
#We have to calculate the loss difference from any two models: VAR or ARIMA, Recursive and 
#Rolling. 
#creating a list of all pairwise comparisons:
prognoser <- c("ARIMA.RES.FE", "VAR.RESS.FE", "ARIMA.RO.FE", "VAR.RO.FE")
(testlista <- combn(prognoser, 2))
rbind(apply(testlista, 2, DM.TEST, 1), apply(testlista, 2, DM.TEST, 2),
      apply(testlista, 2, DM.TEST, 3), apply(testlista, 2, DM.TEST, 4))
#we can see the columns are each of the possible combinations
#and the rows are the horizons. Except for the combination of VAR rolling and VAR recursive,
#for all of them, we reject null of equal precisions for higher horizons( 3 and 4). 
########################################################


