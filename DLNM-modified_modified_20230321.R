############### Distributed Lag Nonlinear Model (DLNM) with example regions
##### After feedback
### 2023-03-21

# 전체기간
# week 변수의 df 값은 2~5로 변경
# temperaure의 lag 12로


library(lubridate)
library(dlnm)
library(foreign)
library(tidyverse)
library(Epi)
library(tsModel)
library(splines)
options(na.action="na.exclude")

load(file = "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\dataset/final_data_전국.rda")


library(readr)
Jeonju <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju.csv")
Jeongeup <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup.csv")
Ulsan <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan.csv")
Jeju <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju.csv")

### 2019년까지만 사용.
Jeonju_2019 <- Jeonju %>% 
  filter(year < 2020)
Jeongeup_2019 <- Jeongeup %>% 
  filter(year < 2020)
Ulsan_2019 <- Ulsan %>% 
  filter(year < 2020)
Jeju_2019 <- Jeju %>% 
  filter(year < 2020)



### 설, 추석이 포함된 주차에 대한 변수 만들기.
Jeonju <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju.csv")

holi <- Jeonju %>% 
  dplyr::select(holi)

write.csv(holi, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi.csv")
holi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi.csv")

Jeonju_2019$holi <- as.factor(holi$holi)
Jeongeup_2019$holi <- as.factor(holi$holi)
Ulsan_2019$holi <- as.factor(holi$holi)
Jeju_2019$holi <- as.factor(holi$holi)




#############################################################
###########  temp*preci / temp lag는 8주까지 / week df = 3, 5 / year trend ns(df = 3)로 
#############################################################


######################### Jeonju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeonju_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeonju_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju_2019$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeonju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeonju_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeonju_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=8, by=1)
pred.preci <- crosspred(cb.preci, model, cen=10, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,40), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeongeup_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeongeup_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup_2019$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeongeup_2019$year,
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeongeup_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeongeup_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeongeup_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=8, by=1)
pred.preci <- crosspred(cb.preci, model, cen=12, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,90), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





######################### Ulsan 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Ulsan_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Ulsan_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan_2019$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Ulsan_2019$year,
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Ulsan_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Ulsan_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Ulsan_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeju_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeju_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju_2019$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeju_2019$year,
                      arglag=list(fun="ns", knots=c(3, 5))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeju_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeju_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeju_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=3, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,60), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 300), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)










#########################################################################################################












#############################################################
###########  temp*preci / temp lag는 12주까지 / week df = 5 / year trend ns(df = 3)로 
#############################################################


######################### Jeonju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeonju_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeonju_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeonju_2019$year,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeonju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeonju_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeonju_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 4) + holi, Jeonju_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=10, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeongeup_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeongeup_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeongeup_2019$year,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeongeup_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeongeup_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeongeup_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=12, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





######################### Ulsan 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Ulsan_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Ulsan_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Ulsan_2019$year,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Ulsan_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Ulsan_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Ulsan_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeju_2019$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeju_2019$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeju_2019$year,
                      arglag=list(fun="ns", knots=c(4, 8))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeju_2019$year,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeju_2019$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeju_2019, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(4, 10, 20, 25), lag=c(0,4,8, 12),
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 300), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju_2019$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)











#########################################################################################################












#############################################################
###########  week는 32~52, 1~4로 / temp*preci / temp lag는 8, 12주까지 / week df = 3, 5 / year trend ns(df = 3)로 
#############################################################
load(file = "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\dataset/final_data_전국.rda")


library(readr)
Jeonju <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju.csv")
Jeongeup <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup.csv")
Ulsan <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan.csv")
Jeju <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju.csv")

### 2019년까지만 사용.
Jeonju_2019 <- Jeonju %>% 
  filter(year < 2020)
Jeongeup_2019 <- Jeongeup %>% 
  filter(year < 2020)
Ulsan_2019 <- Ulsan %>% 
  filter(year < 2020)
Jeju_2019 <- Jeju %>% 
  filter(year < 2020)

### 32주 ~ 52주, 그 다음해 4주까지 사용.
### 53주를 빼는게 에러를 발생시킴. 이유는 진짜 모르겠음.
Jeonju_over32 <- Jeonju_2019 %>% filter(week > 31 | week < 5) 
Jeongeup_over32 <- Jeongeup_2019 %>% filter(week > 31 | week < 5)
Ulsan_over32 <- Ulsan_2019 %>% filter(week > 31 | week < 5)
Jeju_over32 <- Jeju_2019 %>% filter(week > 31 | week < 5)


# 그리고 제일 처음의 4주는 빼자.
Jeonju_over32 <- Jeonju_over32[-c(1:4),]
Jeongeup_over32 <- Jeongeup_over32[-c(1:4),]
Ulsan_over32 <- Ulsan_over32[-c(1:4),]
Jeju_over32 <- Jeju_over32[-c(1:4),]

# 2020년의 1-4주 추가하기.
Jeonju_over32 <- Jeonju_over32 %>% rbind(Jeonju %>% filter(year == 2020 & week < 5))
Jeongeup_over32 <- Jeongeup_over32 %>% rbind(Jeongeup %>% filter(year == 2020 & week < 5))
Ulsan_over32 <- Ulsan_over32 %>% rbind(Ulsan %>% filter(year == 2020 & week < 5))
Jeju_over32 <- Jeju_over32 %>% rbind(Jeju %>% filter(year == 2020 & week < 5))

### 그리고 해당 기간 묶어주는 변수 만들기. 
a <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\a.csv")

Jeonju_over32$s_group <- a$s_group
Jeongeup_over32$s_group <- a$s_group
Ulsan_over32$s_group <- a$s_group
Jeju_over32$s_group <- a$s_group

### 설, 추석이 포함된 주차에 대한 변수 만들기.

holi_over32 <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi_over32.csv")

Jeonju_over32$holi <- as.factor(holi_over32$holi)
Jeongeup_over32$holi <- as.factor(holi_over32$holi)
Ulsan_over32$holi <- as.factor(holi_over32$holi)
Jeju_over32$holi <- as.factor(holi_over32$holi)




######################################## lag 12, df=5

######################### Jeonju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeonju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeonju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju_over32$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeonju_over32$s_group,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeonju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeonju_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeonju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeonju_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,100), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeongeup_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeongeup_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup_over32$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeongeup_over32$s_group,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeongeup_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeongeup_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeongeup_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeongeup_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=8, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





######################### Ulsan 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Ulsan_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Ulsan_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan_over32$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Ulsan_over32$s_group,
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Ulsan_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Ulsan_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Ulsan_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Ulsan_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(0, 10, 20, 25), lag=c(0,4,8, 12), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")


plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,50), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju_over32$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeju_over32$s_group,
                      arglag=list(fun="ns", knots=c(4, 8))) 
cb.preci <- crossbasis(Jeju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeju_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeju_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(4, 10, 20, 25), lag=c(0,4,8, 12),
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,30), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 300), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)


















######################################## lag 8, df=3,5

######################### Jeonju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeonju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeonju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeonju_over32$s_group,
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeonju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeonju_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeonju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeonju_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=8, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,15), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeongeup_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeongeup_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeongeup_over32$s_group,
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeongeup_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeongeup_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeongeup_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeongeup_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=8, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





######################### Ulsan 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Ulsan_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Ulsan_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Ulsan_over32$s_group,
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Ulsan_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Ulsan_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Ulsan_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Ulsan_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")


plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,80), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeju 
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                      group = Jeju_over32$s_group,
                      arglag=list(fun="ns", knots=c(3, 5))) 
cb.preci <- crossbasis(Jeju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       group = Jeju_over32$s_group,
                       arglag=list(fun="ns", knots=c(3,5)))


yr <- length(unique(Jeju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeju_over32, 
             family=quasipoisson)


# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(10, 20, 25), lag=c(0,4,8),
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,30), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 300), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=120, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit[ ,"lag0"] 
pred.temp$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.temp$matRRfit["0" , ] 
pred.temp$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_RR, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju_over32$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)
