############### Distributed Lag Nonlinear Model (DLNM) with example regions


library(dlnm)
library(foreign)
library(tidyverse)
library(Epi)
library(tsModel)
library(splines)
options(na.action="na.exclude")

write.csv(Jeonju, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju.csv")
write.csv(Jeongeup, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup.csv")
write.csv(Ulsan, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan.csv")
write.csv(Jeju, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju.csv")

library(readr)
Jeonju <- read.csv("C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju.csv")
Jeongeup <- read.csv("C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup.csv")
Ulsan <- read.csv("C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan.csv")
Jeju <- read.csv("C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju.csv")

######################### Jeonju 먼저
Jeonju
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeonju$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cbo3int <- crossbasis(Jeonju$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # 이거는 unconstrained
                      arglag=list(fun="integer"))
cb.temp <- crossbasis(Jeonju$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeonju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeonju, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=16, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(-4, 0, 7, 20, 25), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,1.8), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=15, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.temp, "overall", ylim=c(0.8,30.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st


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

write.csv(Jeonju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup
Jeongeup
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeongeup$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeongeup$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeongeup, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=16, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(-4, 0, 7, 20, 25), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,10), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=16, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st


### dcentering 변경해서
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.temp, "overall", ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st


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

write.csv(Jeongeup_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)






#######################     Ulsan  #######################
Ulsan
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Ulsan$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Ulsan$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Ulsan, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=16, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(-1, 0, 7, 20, 25), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,10), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=16, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-1, lty=2, col=4) # 1st


### dcentering 변경해서
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.temp, "overall", ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-1, lty=2, col=4) # 1st


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

write.csv(Ulsan_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




#######################   Jeju   #######################
Jeju
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeju$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeju, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=16, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(4, 7, 12, 20, 25), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,10), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=16, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


### dcentering 변경해서
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.temp, "overall", ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


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

write.csv(Jeju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)









######################### Precipitation

#######################   Jeonju   #######################
Jeonju
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeonju$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeonju$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeonju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeonju, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1) # pct에서 66%인 걸로 적당히


plot(pred.preci, "slices", var=c(0, 10, 20, 100, 200), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall", xlim = c(0, 100), ylim=c(0.9,6), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=100, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st


### centering 변경해서
pred.preci <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.preci, "overall", xlim = c(0, 200), ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_preci_RR <- with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_preci_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju_preci_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)






#######################   Jeongeup   #######################
Jeongeup
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeongeup$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeongeup$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeongeup$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeongeup, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1) # pct에서 66%인 걸로 적당히


plot(pred.preci, "slices", var=c(0, 10, 20, 100, 200), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall", xlim = c(0, 100), ylim=c(0.9,6), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=100, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st


### centering 변경해서
pred.preci <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.preci, "overall", xlim = c(0, 200), ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_preci_RR <- with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_preci_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_preci_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




#######################   Ulsan   #######################
Ulsan
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Ulsan$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Ulsan$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Ulsan$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Ulsan, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1) # pct에서 66%인 걸로 적당히


plot(pred.preci, "slices", var=c(0, 10, 20, 100, 200), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall", xlim = c(0, 100), ylim=c(0.9,6), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=100, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st


### centering 변경해서
pred.preci <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.preci, "overall", xlim = c(0, 200), ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_preci_RR <- with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_preci_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan_preci_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





#######################   Jeju   #######################
Jeju
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeju$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeju$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeju, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=10, by=1) # pct에서 66%인 걸로 적당히


plot(pred.preci, "slices", var=c(0, 10, 20, 100, 200), lag=c(0,2,4,6,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall", xlim = c(0, 200), ylim=c(0.9,6), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=200, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st


### centering 변경해서
pred.preci <- crosspred(cb.temp, model, cen=10, by=1)
plot(pred.preci, "overall", xlim = c(0, 200), ylim=c(0.8,100.0), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=4, lty=2, col=4) # 1st


# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_preci_RR <- with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_preci_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju_preci_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.preci, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)









###################################################
######## Only including months (weeks) with cases
###################################################

Jeonju_weekly <- Jeonju %>% group_by(week) %>% 
  summarise(cases = sum(cases))
Jeongeup_weekly <- Jeongeup %>% group_by(week) %>% 
  summarise(cases = sum(cases))
Ulsan_weekly <- Ulsan %>% group_by(week) %>% 
  summarise(cases = sum(cases))
Jeju_weekly <- Jeju %>% group_by(week) %>% 
  summarise(cases = sum(cases)) 

total_weekly <- cbind(Jeonju_weekly, Jeongeup_weekly$cases, Ulsan_weekly$cases, Jeju_weekly$cases)

total_weekly <- total_weekly %>% 
  mutate(sum_cases = cases + Jeongeup_weekly$cases +
           Ulsan_weekly$cases + Jeju_weekly$cases) %>% 
  dplyr::select(week, sum_cases)

ggplot(aes(x = week, y = sum_cases), data = total_weekly) + geom_point() +
  geom_line() + theme_bw()

write.csv(total_weekly, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\total_weekly.csv")

##################### 32 weeks to 53 weeks. #############################

Jeonju_over32 <- Jeonju %>% filter(week > 31)
Jeongeup_over32 <- Jeongeup %>% filter(week > 31)
Ulsan_over32 <- Ulsan %>% filter(week > 31)
Jeju_over32 <- Jeju %>% filter(week > 31)

Jeonju_over32 <- Jeonju_over32 %>% 
  filter(year < 2020)



######################### Jeonju 
Jeonju_over32
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeonju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju_over32$tmean, lag=12, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeonju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeonju_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=12, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,15), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st



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

write.csv(Jeonju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju_RR.csv")

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
Jeongeup_over32
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeongeup_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeongeup_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeongeup_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=12, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st



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

write.csv(Jeongeup_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_RR.csv")

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
Ulsan_over32
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Ulsan_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Ulsan_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Ulsan_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=12, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,15), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st



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

write.csv(Ulsan_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan_RR.csv")

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
Jeju_over32
######################### Temperature

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeju_over32$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju_over32$tmean, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + ns(week, df = 7*yr) + year, Jeju_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=12, by=1)


d3 <- plot(pred.temp, xlab="Temperature", zlab="RR",
           phi=35, theta=205, ltheta=170, shade=0.4)
lines(trans3d(x=24, y=0:27, z=pred.temp$matRRfit[as.character(24),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=pred.temp$predvar, y=4, z=pred.temp$matRRfit[,"lag4"], pmat=d3), col=3, lwd=2)

plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,15), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=12, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st



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

write.csv(Jeju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju_RR.csv")

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





######################### Jeonju 
Jeonju_over32
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeonju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeonju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeonju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeonju_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=8, by=1)


plot(pred.preci, "slices", var=c(10, 50, 100, 150), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall",xlim = c(0, 200), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=8, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=100, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju_RR.csv")

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
Jeongeup_over32
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeongeup_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeongeup_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeongeup_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeongeup_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.preci, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall",xlim = c(0, 200), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=150, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_RR.csv")

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





######################### Jeongeup 
Jeongeup_over32
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeongeup_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeongeup_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeongeup_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeongeup_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.preci, "slices", var=c(10, 50, 100, 150), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall",xlim = c(0, 200), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=150, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_RR.csv")

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
Ulsan_over32
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Ulsan_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Ulsan_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Ulsan_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Ulsan_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.preci, "slices", var=c(10, 50, 100, 150), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall",xlim = c(0, 200), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=150, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan_RR.csv")

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
Jeju_over32
######################### Precipitation

## internal knots 정해놓고 시작하기.
pct <- quantile(Jeju_over32$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot <- pct[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.preci <- crossbasis(Jeju_over32$total_preci, lag=8, argvar=list(fun="ns", knots=varknot), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...

yr <- length(unique(Jeju_over32$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.preci + ns(week, df = 7*yr) + year, Jeju_over32, 
             family=quasipoisson)

# Prediction and visualization
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.preci, "slices", var=c(10, 50, 100, 150), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.preci, "overall",xlim = c(0, 200), ylim=c(0.9,10), col=2, lwd=2, xlab="Precipitation", ylab="RR",
     main="Overall cumulative association")
abline(v=5, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=150, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st



# Lag-specific RRs at different temperatures 
# lag를 고정하고 temp의 변화에 따른 RR 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit[ ,"lag0"] 
pred.preci$matRRfit[ ,"lag7"]

# Temperature-specific RRs at different lags 
# temp를 고정하고 lag에 따른 RR의 변화를 살펴보기
# Exponentiated point estimates; See also 95% CIs in matRRlow and matRRhigh
# These estiamtes correspond to the 'sliced plots'
pred.preci$matRRfit["0" , ] 
pred.preci$matRRfit["26", ] 

# Overall cummulative (lag-cummulative) association 
# lag는 통합적인 영향으로 고려되고, temp(메인 설명변수)에 따른 RR이 계산된다.
# Exponentiated point estimates & interval estimates 
Jeju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju_RR.csv")

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





#############################################################
###########  temp*preci 해보기. precipitation의 값이 너무 튀어서, 그냥 전체 기간에 대해서 interaction으로.
#############################################################


######################### Jeonju 
Jeonju
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeonju$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeonju$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeonju$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...
cb.preci <- crossbasis(Jeonju$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5)))

yr <- length(unique(Jeonju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 7*yr) + year, Jeonju, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,15), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,5), col=2, lwd=2, xlab="Precipitation", ylab="RR",
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
Jeonju_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeonju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeonju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeonju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Jeongeup 
Jeongeup
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeongeup$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeongeup$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeongeup$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...
cb.preci <- crossbasis(Jeongeup$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5)))

yr <- length(unique(Jeongeup$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 7*yr) + year, Jeongeup, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=-4, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,5), col=2, lwd=2, xlab="Precipitation", ylab="RR",
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
Jeongeup_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Jeongeup_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeongeup_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeongeup$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)




######################### Ulsan 
Ulsan
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Ulsan$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Ulsan$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Ulsan$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...
cb.preci <- crossbasis(Ulsan$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5)))

yr <- length(unique(Ulsan$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 7*yr) + year, Ulsan, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=0, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,5), col=2, lwd=2, xlab="Precipitation", ylab="RR",
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
Ulsan_RR <- with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))

write.csv(Ulsan_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Ulsan_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Ulsan$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)





######################### Jeju 
Jeju
######################### Temperature + Precipitaion

## internal knots 정해놓고 시작하기.
pct_temp <- quantile(Jeju$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th

pct_preci <- quantile(Jeju$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th


## model에 넣을 주요 설명변수인 temperature를 ns로 처리하기 위한 것이다.
# 먼저 위의 varknot에 담아놓은, temperature 자체의 ns 적용을 위한 knot와,
# lag effect의 ns를 설명하기 위해서 설정하는 arglag를 인수로해서 
# crossbasis 객체를 생성한다. 
cb.temp <- crossbasis(Jeju$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), # weekly 자료인데 lag 2달 줄거니까 8로
                      arglag=list(fun="ns", knots=c(3,5))) # lag가 8이니까 knots도 그 안에서 정해야하는데...
cb.preci <- crossbasis(Jeju$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), # weekly 자료인데 lag 2달 줄거니까 8로
                       arglag=list(fun="ns", knots=c(3,5)))

yr <- length(unique(Jeju$year))
# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 7*yr) + year, Jeju, 
             family=quasipoisson)

# Prediction and visualization
pred.temp <- crosspred(cb.temp, model, cen=10, by=1)
pred.preci <- crosspred(cb.preci, model, cen=5, by=1)


plot(pred.temp, "slices", var=c(7, 12, 20, 25), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

plot(pred.preci, "slices", var=c(5, 10, 50, 100), lag=c(0,3,5,8), 
     ci="area", ci.arg=list(col=grey(0.7)), type="l", col=2, ylab="RR")

# 여기서 percentile은 위에서 pct에서 계산된 값 넣는 것.
par(mfrow=c(1,1))
plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
     main="Overall cumulative association")
abline(v=10, lty=2, col="grey") # MMT minimum mortality temperature
abline(v=28, lty=2, col=2) # 99th
abline(v=3, lty=2, col=4) # 1st

plot(pred.preci, "overall", xlim = c(0, 150), ylim=c(0.9,5), col=2, lwd=2, xlab="Precipitation", ylab="RR",
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

write.csv(Jeju_RR, "C:\\Users\\Taehee Chang\\Desktop\\Backup\\Univ Tokyo papers\\Calculated_data\\Jeju_RR.csv")

# minimum mortality temperature (MMT) or the optimal temperature (OT)인 16도를 기준으로
# 각 pct 온도에서 나타나는 RR 값을 계산.
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["19",] 
with(pred.temp, cbind(allRRfit, allRRlow, allRRhigh))["1",]


########## Residuals
res.spline <- residuals(model, type="response")
plot(Jeju$year, res.spline ,ylim=c(-50,150),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/ time adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)


library(lubridate)
holi <- date(c("2019-09-12", "2017-10-01"))
week(holi)