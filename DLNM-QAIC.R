############### Distributed Lag Nonlinear Model (DLNM) with example regions
##### After feedback
### 2023-04-24



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
holi <- holi[1:991,]

Jeonju_2019$holi <- as.factor(holi$holi)
Jeongeup_2019$holi <- as.factor(holi$holi)
Ulsan_2019$holi <- as.factor(holi$holi)
Jeju_2019$holi <- as.factor(holi$holi)


################################################################################
## Model selection
################################################################################

# FUNCTION OF QAIC
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois( model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

lmodels <- list(mod.naive, mod.adjtrend, mod.adjt1, mod.adjtemp)
sapply(lmodels, QAIC)


#############################################################
###########  Jeonju 
#############################################################

######################### lag 8 / df = 3, 4, 5, 6

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


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeonju_2019, 
             family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeonju_2019, 
             family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeonju_2019, 
             family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeonju_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeonju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeonju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeonju_2019, 
                 family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeonju_2019, 
                 family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeonju_2019, 
                 family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeonju_2019, 
                 family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)



#############################################################
###########  Jeongeup 
#############################################################

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeongeup_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeongeup_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeongeup_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeongeup_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeongeup_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeongeup_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeongeup_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeongeup_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeongeup_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)



#############################################################
###########  Ulsan 
#############################################################

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Ulsan_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Ulsan_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Ulsan_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Ulsan_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Ulsan_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Ulsan_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Ulsan_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Ulsan_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Ulsan_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)






#############################################################
###########  Jeju 
#############################################################

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeju_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeju_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeju_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeju_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, Jeju_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, Jeju_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi, Jeju_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi, Jeju_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)









##############################################################################################
### Autocorrelation


#############################################################
###########  Jeonju 
#############################################################
Jeonju_2019$auto <- log(lag(Jeonju_2019$cases, 1) + 0.01)
Jeonju_2019[is.na(Jeonju_2019)] <- log(0.01)

######################### lag 8 / df = 3, 4, 5, 6

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


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeonju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeonju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi + auto, Jeonju_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)



#############################################################
###########  Jeongeup 
#############################################################

Jeongeup_2019$auto <- log(lag(Jeongeup_2019$cases, 1) + 0.01)
Jeongeup_2019[is.na(Jeongeup_2019)] <- log(0.01)

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeongeup_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeongeup_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Jeongeup_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)



#############################################################
###########  Ulsan 
#############################################################

Ulsan_2019$auto <- log(lag(Ulsan_2019$cases, 1) + 0.01)
Ulsan_2019[is.na(Ulsan_2019)] <- log(0.01)

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Ulsan_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Ulsan_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Ulsan_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)






#############################################################
###########  Jeju 
#############################################################

Jeju_2019$auto <- log(lag(Jeju_2019$cases, 1) + 0.01)
Jeju_2019[is.na(Jeju_2019)] <- log(0.01)

######################### lag 8 / df = 3, 4, 5, 6

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
                      arglag=list(fun="ns", knots=c(3,5))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_8_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                 family=quasipoisson)
model_8_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                 family=quasipoisson)
model_8_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                 family=quasipoisson)
model_8_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                 family=quasipoisson)


######################### lag 12 / df = 3, 4, 5, 6

cb.temp <- crossbasis(Jeju_2019$tmean, lag=12, argvar=list(fun="ns", knots=varknot_temp), 
                      arglag=list(fun="ns", knots=c(4,8))) 
cb.preci <- crossbasis(Jeju_2019$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                       arglag=list(fun="ns", knots=c(3,5)))


# 그리고 여기에는 위에서 만든 crossbasis 객체와 covariates를 넣고 모델을 적합한다.
model_12_3 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                  family=quasipoisson)
model_12_4 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                  family=quasipoisson)
model_12_5 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 5) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                  family=quasipoisson)
model_12_6 <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 6) + ns(year, df = 3) + holi+ auto, Jeju_2019, 
                  family=quasipoisson)

lmodels <- list(model_8_3, model_8_4, model_8_5, model_8_6, model_12_3, model_12_4, model_12_5, model_12_6)
sapply(lmodels, QAIC)


