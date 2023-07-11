# install.packages("xts")
install.packages("BETS")
library(xts)
library(fpp2)
library(dplyr)
library(tidyverse)
library(tsbox)
library(urca)
library(tseries)
library(BETS)
setwd("/Users/camelliasun/Desktop/time series/final project")

# Data from 2017-05-01 to 2022-08-07
litecoin<- read.csv("CBLTCUSD.csv")
litecoin <- litecoin[-1, ]

# Transfer the data from daily to weekly
litecoin$CBLTCUSD <- as.numeric(litecoin$CBLTCUSD)
litecoin <- na.omit(litecoin)
data <- as.xts(litecoin$CBLTCUSD,order.by=as.Date(litecoin$DATE))
weekly <- apply.weekly(data,mean)
colnames(weekly) <- c("Price")
#weekly_test <- as.ts(weekly)
weekly.ts <- ts_ts(weekly)
autoplot(weekly.ts, ylab = "$ Price") +
  ggtitle("Time Series of Weekly Average Price of Litecoin")
length(weekly.ts)

plot(decompose(weekly.ts))

# Train and test data
weekly[220]
train.0 <- window(weekly.ts, start =2017.3446347032, end = 2021.557)
test <- window(weekly.ts, start = 2021.557)

length(test)
length(train.0)
length(weekly)


# Transfer the data from daily to weekly(not used)
#litecoin$CBLTCUSD <- as.numeric(litecoin$CBLTCUSD)
#litecoin$seven_day_index <- c(1, rep(1:(nrow(litecoin)-1)%/%7)+1)
#week_litecoin <- group_by(litecoin, seven_day_index) 
#week_litecoin <- summarize(week_litecoin, mean_prices = mean(CBLTCUSD, na.rm = TRUE))
#week_litecoin <- as.data.frame(week_litecoin)

# Plot of train data
autoplot(train.0, ylab = "$ Weekly average price") +
  ggtitle("Weekly Average Price of Traning Data")
str(train.0)
# Check if Box-cox transformation is necessary  
lambda <- BoxCox.lambda(train.0) 
autoplot(BoxCox(train.0, lambda),ylab = "$ Weekly average price") +
  ggtitle("Weekly Average Price of Traning Data After Transformataion")
train.0 <- BoxCox(train.0, lambda)

lambda

# ACF, ACF and PACF after differencing
Acf(weekly.ts)
train.0  %>%  diff() %>% ggtsdisplay(main="")
ggseasonplot(weekly.ts) +
  ggtitle("Seasonal Plot of Training Data")
# In this case, we expect to make first difference
train.0 %>% adf.test()
train.0 %>% diff() %>% adf.test()
# for loop
n = length(train.0)
max.p = 2
max.d = 1
max.q = 2
max.P = 2
max.D = 2
max.Q = 1

BIC.array =array(NA,dim=c(max.p+1,max.d+1,max.q+1,max.P+1,max.D+1,max.Q+1))
AIC.array =array(NA,dim=c(max.p+1,max.d+1,max.q+1,max.P+1,max.D+1,max.Q+1))

best.bic <- 1e8
x.ts = train.0

for (p in 0:max.p) for(d in 1:max.d) for(q in 0:max.q) 
  for (P in 0:max.P) for(D in 0:max.D) for(Q in 0:max.Q) 
  {
#    cat("p=",p,", d=",d,", q=",q,", P=",P,", D=",D,", Q=",Q,"\n")
    
    fit <- tryCatch(
      {  arima(x.ts, order = c(p,d,q),  
               seas = list(order = c(P,D,Q), 
                           frequency(x.ts)),method="CSS-ML")
      },
      error = function(cond){
        # Choose a return value in case of error
        return(NA)
      }
    )
      if(is.list(fit)){
      number.parameters <- length(fit$coef) + 1
      BIC.array[p+1,d+1,q+1,P+1,D+1,Q+1] = -2*fit$loglik + log(n)*number.parameters
      AIC.array[p+1,d+1,q+1,P+1,D+1,Q+1] = -2*fit$loglik + 2*number.parameters
      
      arima.aic <- AIC.array[p+1,d+1,q+1,P+1,D+1,Q+1]
      arima.bic <- BIC.array[p+1,d+1,q+1,P+1,D+1,Q+1]
      arima.fit <- fit
      arima.model <- c(p,d,q,P,D,Q) 
      
      cat("arima.aic=",arima.aic,", arima.bic=",arima.bic,",arima.model=",arima.model,"\n")
      
    } 
  }

#########



# auto.arima gives model with parameters: 011100
auto.arima(train.0)
checkresiduals(auto.arima(train.0))

# cross validation
# use best BIC: arima.model= 010000 ; 110000 ; 011000 ; 010001 ; 010100 
# where smallest BIC is 1 1 1


library(forecast)

# sarima 010000
farima.010.000 <- function(x, h){
  forecast(Arima(x, order=c(0,1,0),
                 seas=list(order=c(0,0,0),period=52)), h = h)
}

e.010.000 <- tsCV(BoxCox(weekly.ts, lambda), farima.010.000, h=1)
sqrt(mean(e.010.000^2, na.rm=TRUE))

# sarima 110000
farima.110.000 <- function(x, h){
  forecast(Arima(x, order=c(1,1,0),
                 seas=list(order=c(0,0,0),period=52)), h = h)
}

e.110.000 <- tsCV(BoxCox(weekly.ts, lambda), farima.110.000, h=1)
sqrt(mean(e.110.000^2, na.rm=TRUE))

# sarima 011000
farima.011.000 <- function(x, h){
  forecast(Arima(x, order=c(0,1,1),
                 seas=list(order=c(0,0,0),period=52)), h = h)
}

e.011.000 <- tsCV(BoxCox(weekly.ts, lambda), farima.011.000, h=1)
sqrt(mean(e.011.000^2, na.rm=TRUE))

# sarima 010001
farima.010.001 <- function(x, h){
  forecast(Arima(x, order=c(0,1,0),
                 seas=list(order=c(0,0,1),period=52)), h = h)
}

e.010.001 <- tsCV(BoxCox(weekly.ts, lambda), farima.010.001, h=1)
sqrt(mean(e.010.001^2, na.rm=TRUE))

# sarima 010100
farima.010.100 <- function(x, h){
  forecast(Arima(x, order=c(0,1,0),
                 seas=list(order=c(1,0,0),period=52)), h = h)
}

e.010.100 <- tsCV(BoxCox(weekly.ts, lambda), farima.010.100, h=1)
sqrt(mean(e.010.100^2, na.rm=TRUE))

# sarima 011100
farima.011.100 <- function(x, h){
  forecast(Arima(x, order=c(0,1,1),
                 seas=list(order=c(1,0,0),period=52)), h = h)
}

e.011.100 <- tsCV(BoxCox(weekly.ts, lambda), farima.011.100, h=1)
sqrt(mean(e.011.100^2, na.rm=TRUE))

# check WN
sarima.010.100 <- arima(train.0,order=c(0,1,0),
      seas=list(order=c(1,0,0),period=52))

forecast <- predict(sarima.010.100,n.ahead=55)
ggtsdisplay(forecast$pred)


# RMSE = 0.009610696, 0.009768071, 0.009811072, 0.009621089, 0.008935617 respectively.
# Arima(010100) has the smallest cross-validation error when h = 1 (and h = 5). 
# So we would like to choose Arima(010100) as the best fiting model.

# Check residual to see whether the residuals are behaving like white noise
sarima.010.100 <- Arima(train.0,order=c(0,1,0),
                        seas=list(order=c(1,0,0),period=52))
#arima.101 <- Arima(train,order=c(1,0,1))
#arima.112 <- Arima(train,order=c(1,1,2))
#arima.211 <- Arima(train,order=c(2,1,1))
#arima.121 <- Arima(train,order=c(1,2,1))

checkresiduals(sarima.010.100)

sarima.011.100 <- Arima(train.0,order=c(0,1,1),
                        seas=list(order=c(1,0,0),period=52))
checkresiduals(sarima.011.100)

# Check characteristic roots. Points lie within the unit circle indicating stationarity
# of residuals.

# fitted VS train.0(transfer back to original series)

plot(InvBoxCox(train.0, lambda = lambda), col = "red")
lines(InvBoxCox(fitted(sarima.010.100), lambda = lambda), col = "blue")
Arima(train.0, order = c(0, 1, 0), seasonal=c(0, 0, 1))

# forecast vs test
sarima.010.100 <- Arima(train.0,order=c(0,1,0),
                        seas=list(order=c(1,0,0),period=52))

forecast.1 <- predict(sarima.010.100,n.ahead=55)
autoplot(InvBoxCox(forecast.1$pred, lambda))+autolayer(test)

sarima.011.100 <- arima(train.0,order=c(0,1,1),
                        seas=list(order=c(1,0,0),period=52))

forecast.2 <- predict(sarima.011.100,n.ahead=55)
autoplot(InvBoxCox(forecast.2$pred, lambda))+autolayer(test)

plot(InvBoxCox(train.0, lambda = lambda), col = "red")
lines(InvBoxCox(forecast.2$pred, lambda = lambda), col = "green")
lines(weekly.ts)

autoplot(weekly.ts, series = "Original Series") +
  autolayer(InvBoxCox(forecast.2$pred, lambda), series = "Forecasting")


# ARIMA Forecast Accuracy

arima.010.100F <- forecast(sarima.010.100, h = 55)
accuracy(arima.010.100F, BoxCox(test, lambda = lambda))


# ETS model

# ets with parameters A A N with damped 
ets_fit <- ets(train.0)
summary(ets_fit) 

autoplot(ets_fit)

# ets using cross validation to find its cross-validation error
fets <- function(x, h){
  forecast(ets(x, model = "AAN", damped = TRUE), h = h)
}
e_ets <- tsCV(BoxCox(weekly.ts, lambda), fets, h = 1) 
sqrt(mean(e_ets^2, na.rm=TRUE))

# Check whether the fitted ets model has residuals being white noise
# Ljung-Box test rejects residuals being white noise.
checkresiduals(ets_fit)

# ETS Forecast accuracy
atsfa <- forecast(ets_fit,h=55)
accuracy(atsfa,BoxCox(test, lambda = lambda))


# autoarima accuracy

fcauto <- forecast(auto.arima(train.0), h=55)
accuracy(fcauto,BoxCox(test, lambda = lambda))


# confidence interval
forecast.1 <- predict(sarima.010.100,n.ahead=55)
plot(forecast.1$pred,ylab="Predicted processed Litecoin Price",ylim=(c(0,100)))
lines(forecast.1$pred+1.96*forecast.1$se,lty=2,col="red")
lines(forecast.1$pred-1.96*forecast.1$se,lty=2,col="red")



plot(InvBoxCox(forecast.1$pred, lambda),ylab="Predicted Litecoin Price",ylim=(c(-100,1000)))
lines(InvBoxCox(forecast.1$pred+1.96*forecast.1$se, lambda),lty=2,col="red")
lines(InvBoxCox(forecast.1$pred-1.96*forecast.1$se, lambda),lty=2,col="red")


train.0 %>% diff(lag = 12) %>% ggtsdisplay(main="")
summary(sarima.010.100)
coeftest(sarima.010.100)
coeftest(sarima.011.100)
t_test(sarima.010.100)
t_test(sarima.011.100)
