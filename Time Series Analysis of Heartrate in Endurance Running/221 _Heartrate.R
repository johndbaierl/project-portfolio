##########################
####Load 11/20/22 run data
##########################

dat2 = read.csv(file = "Run Data 11_20_22.csv")
colnames(dat2) = cols
#View(dat2)

dat2$Distance = dat2$Distance * 0.6213711922 #Convert km to mi

#trim dist > 1 and dist < 10.81
run2.dat = dat2[dat2$Distance > 1 & dat2$Distance < 10.81, ]

run2.dat$Time = run2.dat$Time - min(run2.dat$Time) + 1

plot(Heartrate ~ Time, data = run2.dat, type = 'l',
     xlab = "Time (s)", ylab = "Heartrate (BPM)",
     main = "Heartrate vs. Time at Constant Percieved Effort")
abline(v = cut, lty = 2)


###75-25 train-test split:

cut = floor(max(run2.dat$Time) * 0.75)

run2.train = run2.dat[1:cut, ]
run2.test = run2.dat[(cut + 1):max(run2.dat$Time), ]


###Detrending
lm.fit = lm(Heartrate ~ Time, data = run2.train)
lm.a = coefficients(lm.fit)[1]
lm.b = coefficients(lm.fit)[2]

plot(Heartrate ~ Time, data = run2.train, type = 'l',
     xlab = "Time (s)", ylab = "Heartrate (BPM)",
     main = "Heartrate vs. Time: Training Data")
abline(lm.fit, col = "red")

hr.res = lm.fit$residuals
plot(lm.fit$residuals, type = 'l', xlab = "Time (s)", ylab = "Residuals (BPM)", 
     main = "Detrended HR vs. Time: Training Data")
abline(h = 0, lty = 2)

library(tseries)
adf.test(hr.res)

###acf/pcf
acf(hr.res, main = "ACF: Training Data")

pacf(hr.res, main = "PACF: Training Data")
#plot(pacf(hr.res, plot = F)[2:35])

###Check for cycles...
spectrum(hr.res, method = "ar", main = "Spectral Density: Training Data") #consistent with pretty smooth data

###Fit ARMA, perform grid search over p's and q's:
library(sarima)
library(astsa)
#test
a1 = arima(hr.res, order = c(4, 0, 0))
a1$coef

aics = matrix(nrow = 10, ncol = 10)
for(p in 0:9) {
  for(q in 0:9) {
    fit = arima(hr.res, order = c(p, 0, q))
    aics[p + 1, q + 1] = fit$aic
  }
}

min(na.omit(aics))
which(aics == min(na.omit(aics)), arr.ind = T) #min at p = 3, q = 2

arma.fit = arima(hr.res, order = c(3, 0, 2))
#arma.fit2 = sarima(hr.res, p = 3, d = 0, q = 2)

ar.coef = arma.fit$coef[1:3]
ma.coef = arma.fit$coef[4:5]
int = arma.fit$coef[6]


###Predictions on training data: One-step-ahead

#ARMA(3, 2)
n.train = length(hr.res)
resid = rep(0, n.train)
x.pred = vector(length = n.train)
for(i in 4:n.train) {
  x.pred[i] = sum(ar.coef * hr.res[(i - 3):(i - 1)]) + 
    sum(ma.coef * resid[(i - 1):(i - 2)]) + int
  resid[i] = x.pred[i] - hr.res[i]
}

#just residuals
plot(x = 4:n.train, y = hr.res[4:n.train], type = 'l', xlim = c(1000, 2000))
lines(x = 4:n.train, y = x.pred[4:n.train], col = "red")

#adding in linear effect
x.lin = lm.a + (1:n.train) * lm.b
x.tot = x.pred + x.lin

plot(x = 4:n.train, y = run2.train$Heartrate[4:n.train], type = 'l',
     xlim = c(400, 500))
lines(x = 4:n.train, y = x.tot[4:n.train], col = "red")

#Try an AR(3)
ar3.fit = arima(hr.res, order = c(3, 0, 0))
ar3.coef = ar3.fit$coef[1:3]
ar3.int = ar3.fit$coef[4]

x.pred2 = vector(length = n.train)
for(i in 4:n.train) {
  x.pred2[i] = sum(ar3.coef * hr.res[(i - 3):(i - 1)]) + ar3.int
  resid2[i] = x.pred2[i] - hr.res[i]
}

x.tot2 = x.pred2 + x.lin
plot(x = 4:n.train, y = run2.train$Heartrate[4:n.train], type = 'l',
     xlim = c(400, 500))
lines(x = 4:n.train, y = x.tot2[4:n.train], col = "red")


###RMSE on test data to compare models: One-step-ahead

#ARMA(3, 2)
hr.test = run2.test$Heartrate
res.test = hr.test - (lm.a + lm.b * run2.test$Time)

n.test = length(hr.test)
resid = rep(0, n.test)
x.pred = vector(length = n.test)
for(i in 4:n.test) {
  x.pred[i] = sum(ar.coef * res.test[(i - 3):(i - 1)]) + 
    sum(ma.coef * resid[(i - 1):(i - 2)]) + int
  resid[i] = x.pred[i] - res.test[i]
}

x.lin = lm.a + run2.test$Time * lm.b
x.tot = x.pred + x.lin

plot(x = run2.test$Time[4:n.test], y = hr.test[4:n.test], type = 'l'
     , xlim = c(3700, 3800))
lines(x = run2.test$Time[4:n.test], y = x.tot[4:n.test], col = "red")

plot(x = run2.test$Time[4:n.test], y = res.test[4:n.test], type = 'l',
     xlim = c(3600, 3700), 
     main = "ARMA(3,2) One-Step-Ahead Prediction: 100 s Test Sample",
     xlab = "Time (s)", ylab = "HR Residual (BPM)")
lines(x = run2.test$Time[4:n.test], y = x.pred[4:n.test], col = "red")

sqrt(sum((hr.test[4:n.test] - x.tot[4:n.test])^2) / (n.test - 3))

#AR(3)
ar3.fit = arima(hr.res, order = c(3, 0, 0))
ar3.coef = ar3.fit$coef[1:3]
ar3.int = ar3.fit$coef[4]

x.pred2 = vector(length = n.test)
resid2 = vector(length = n.test)
for(i in 4:n.test) {
  x.pred2[i] = sum(ar3.coef * res.test[(i - 3):(i - 1)]) + ar3.int
  resid2[i] = x.pred2[i] - res.test[i]
}

x.tot2 = x.pred2 + x.lin

plot(x = run2.test$Time[4:n.test], y = hr.test[4:n.test], type = 'l')
lines(x = run2.test$Time[4:n.test], y = x.tot2[4:n.test], col = "red")

sqrt(sum((hr.test[4:n.test] - x.tot2[4:n.test])^2) / (n.test - 3))

#Linear model only:
sqrt(sum((hr.test - x.lin)^2) / n.test)

#Previous second:
hr.shift = c(0, hr.test)
sqrt(sum((hr.test[2:n.test] - hr.shift[2:n.test])^2) / (n.test - 1))


####Forecasting further ahead (30 s):
forecast_30 = predict(arma.fit, n.ahead = 30)$pred
forecast_30 = forecast_30 + (lm.a + lm.b * (3379:(3379 + 29)))
forecast_30.se = predict(arma.fit, n.ahead = 30)$se

U = forecast_30 + forecast_30.se; L = forecast_30 - forecast_30.se
xx = c(time(U), rev(time(U))); yy = c(L, rev(U))

plot(x = 3300:n.train, y = run2.train$Heartrate[3300:n.train], 
     type = 'l', xlim = c(3300, 3410), ylim = c(163, 170),
     main = "Forecasting: 30 Seconds Ahead", xlab = "Time (s)", 
     ylab = "Heartrate (BPM)")
lines(x = (n.train + 1):(n.train + 30), y = forecast_30, type = 'b')
lines(x = (n.train + 1):(n.train + 30), y = run2.test$Heartrate[1:30],
      col = "red")
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
abline(v = cut, lty = 2)

####Forecasting further ahead (60 s):

#ARMA(3, 2)
forecast_60 = predict(arma.fit, n.ahead = 60)$pred
forecast_60 = forecast_60 + (lm.a + lm.b * (3379:(3379 + 59)))
forecast_60.se = predict(arma.fit, n.ahead = 60)$se

U = forecast_60 + forecast_60.se ; L = forecast_60 - forecast_60.se
xx = c(time(U), rev(time(U))) ; yy = c(L, rev(U))

plot(x = 3300:n.train, y = run2.train$Heartrate[3300:n.train], 
     type = 'l', xlim = c(3300, 3440), ylim = c(163, 170),
     main = "ARMA(3, 2) Forecasting: 60 Seconds Ahead", xlab = "Time (s)", 
     ylab = "Heartrate (BPM)")
lines(x = (n.train + 1):(n.train + 60), y = forecast_60, type = 'b')
lines(x = (n.train + 1):(n.train + 60), y = run2.test$Heartrate[1:60],
      col = "red")
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
abline(v = cut, lty = 2)

#AR(3)
ar.forecast_60 = predict(ar3.fit, n.ahead = 60)$pred
ar.forecast_60 = ar.forecast_60 + (lm.a + lm.b * (3379:(3379 + 59)))
ar.forecast_60.se = predict(arma.fit, n.ahead = 60)$se

####Forecasting: 60-second RMSE

#ARMA(3, 2)
sqrt(sum((run2.test$Heartrate[1:60] - forecast_60)^2) / 60)

#Most recent value:
sqrt(sum((run2.test$Heartrate[1:60] - run2.train$Heartrate[cut])^2) / 60)

#lm
forecast_lm = lm.a + lm.b * ((n.train + 1):(n.train + 60))
sqrt(sum((run2.test$Heartrate[1:60] - forecast_lm)^2) / 60)

#AR(3)
sqrt(sum((run2.test$Heartrate[1:60] - ar.forecast_60)^2) / 60)


####Trying out some more...

#Forecasting further ahead (120 s):

forecast_120 = predict(arma.fit, n.ahead = 120)$pred
forecast_120 = forecast_120 + (lm.a + lm.b * (3379:(3379 + 119)))
forecast_120.se = predict(arma.fit, n.ahead = 120)$se

U = forecast_120 + forecast_120.se; L = forecast_120 - forecast_120.se
xx = c(time(U), rev(time(U))); yy = c(L, rev(U))

plot(x = 3300:n.train, y = run2.train$Heartrate[3300:n.train], 
     type = 'l', xlim = c(3300, 3500), ylim = c(163, 170),
     main = "Forecasting: 60 Seconds Ahead", xlab = "Time (s)", 
     ylab = "Heartrate (BPM)")
lines(x = (n.train + 1):(n.train + 120), y = forecast_120, type = 'b')
lines(x = (n.train + 1):(n.train + 120), y = run2.test$Heartrate[1:120],
      col = "red")
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
abline(v = cut, lty = 2)

