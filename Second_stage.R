############## Second stage
options(na.action="na.exclude")

## ------------------------------------------------------------------------
# Load R packages 
library(dplyr)
library(splines)
library(tsModel)
library(Epi)
library(dlnm)
library(metafor)
library(mixmeta)
library(pscl)

## rda 파일은 그냥 load 써야함. 화살표로 객체에 할당하면 안됨.
load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

final_data_1 <- final_data[c("Jeonju", "Jeongeup", "Suncheon", "Daegu", "Ansoeng", "Hwaseong")]
final_data_2 <- final_data[c("Ulsan", "Busan", "Jeju", "Gimhae", "Haman", "Hapcheon")]

final_data_1_1 <- final_data_1 
final_data_2_1 <- final_data_2 

######################################## final data 조정해서 분석 가능한 형태로 바꾸기
### 2019년까지만 사용.

### 그리고 해당 기간 묶어주는 변수 만들기. 
### 설, 추석이 포함된 주차에 대한 변수 만들기.
Jeonju_over32 <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_over32.csv")

holi <- Jeonju_over32 %>% 
  dplyr::select(holi)

a <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\a.csv")


for (i in 1:length(final_data_1)){
  final_data_1[[i]] <- final_data_1[[i]] %>% filter(year < 2020) %>% filter(week > 31 | week < 5)
  final_data_1[[i]] <- final_data_1[[i]][-c(1:4),] 
  final_data_1[[i]] <- final_data_1[[i]] %>% rbind(final_data_1_1[[i]] %>% filter(year == 2020 & week < 5))
  final_data_1[[i]]$s_group <- a$s_group
  final_data_1[[i]]$holi <- holi$holi
  }

for (i in 1:length(final_data_2)){
  final_data_2[[i]] <- final_data_2[[i]] %>% filter(year < 2020) %>% filter(week > 31 | week < 5)
  final_data_2[[i]] <- final_data_2[[i]][-c(1:4),] 
  final_data_2[[i]] <- final_data_2[[i]] %>% rbind(final_data_2_1[[i]] %>% filter(year == 2020 & week < 5))
  final_data_2[[i]]$s_group <- a$s_group
  final_data_2[[i]]$holi <- holi$holi
  }








################ 12월까지만 쓰는 함수.
######################################## final data 조정해서 분석 가능한 형태로 바꾸기
### 2019년까지만 사용. 

### 그리고 해당 기간 묶어주는 변수 만들기. 
### 설, 추석이 포함된 주차에 대한 변수 만들기.
Jeonju_over32 <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\Jeonju_over32.csv")

holi <- Jeonju_over32 %>% 
  filter(week > 31) %>% 
  dplyr::select(holi)

a <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\a.csv")


for (i in 1:length(final_data_1)){
  final_data_1[[i]] <- final_data_1[[i]] %>% filter(year < 2020) %>% filter(week > 31)
  final_data_1[[i]]$holi <- holi$holi
}

for (i in 1:length(final_data_2)){
  final_data_2[[i]] <- final_data_2[[i]] %>% filter(year < 2020) %>% filter(week > 31)
  final_data_2[[i]]$holi <- holi$holi
}

for (i in 1:length(final_data_1)){
  final_data_1[[i]]$holi <- factor(final_data_1[[i]]$holi, levels = c(0,1))
}
for (i in 1:length(final_data_2)){
  final_data_2[[i]]$holi <- factor(final_data_2[[i]]$holi, levels = c(0,1))
}



#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### 경남 이외 지역 final_data_1
regions_1 <- names(final_data_1)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_1))



par(mfrow=c(2,5))
for(i in seq(regions_1)){
  cat(i,"")
  sub <- final_data_1[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 20) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = round(8) # reference temperature at 33th percentile   
  cen_preci = round(5)
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_1[i])
  plot(pred.preci, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_1[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(19) # 90th pct
  target_preci <- as.character(50) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}



# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions_1, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on tsutsugamushi (5 C° vs. 19 C°)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_1, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (5 mm vs. 50 mm)")






#################### 경남 포함 지역 final_data_2
regions_2 <- names(final_data_2)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_2))




for(i in seq(regions_2)){
  cat(i,"")
  sub <- final_data_2[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 20) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = round(8) # reference temperature at 33th percentile   
  cen_preci = round(5)
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_2[i])
  plot(pred.preci, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_2[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(19) # 90th pct
  target_preci <- as.character(50) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}



# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions_1, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on tsutsugamushi (5 C° vs. 19 C°)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_1, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (5 mm vs. 50 mm)")






################################################################################
# Second-stage modelling - Multivariate meta-analysis
################################################################################

coef_temp <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
               dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp <- vector("list", length(regions_1))

coef_preci <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                    dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci <- vector("list", length(regions_1))


for(i in seq(regions_1)){
  cat(i,"")
  sub <- final_data_1[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 20) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model, cen = 8) 
  coef_temp[i,] <- coef(cr_temp)
  vcov_temp[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model, cen = 5) 
  coef_preci[i,] <- coef(cr_preci)
  vcov_preci[[i]] <- vcov(cr_preci)
}


# Multivariate mixed-effects meta-analysis - Temperature
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound <- rowMeans(sapply(dlist, function(x) range(x$tmean)))
xvar <- seq(bound[1],bound[2],by=0.1)
argvar=list(fun="ns", knots=quantile(xvar,prob=c(.33,.66)))
bvar <- do.call(onebasis, c(list(x=xvar), argvar))

pred.pool <- crosspred(bvar, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=5)
pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar, coef=coef[i,], vcov=vcov[[i]], 
                                         model.link="log", cen=5))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,10), lwd=2,
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_temp <- mixmeta(coef_temp~1, vcov_temp, method="ml")
summary(mv_temp)

# Plot the pooled exposure-response curve
bound <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar <- seq(bound[1],bound[2],by=0.1)
argvar=list(fun="ns", knots=quantile(xvar,prob=c(.33,.66)))
bvar <- do.call(onebasis, c(list(x=xvar), argvar))

pred.pool <- crosspred(bvar, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=8)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar, coef=coef[i,], vcov=vcov[[i]], 
                                         model.link="log", cen=8))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2,
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)




##################################################################################
######### 경남 포함 지역 final_data_2

coef_temp <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                    dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp <- vector("list", length(regions_2))

coef_preci <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                     dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci <- vector("list", length(regions_2))


for(i in seq(regions_2)){
  cat(i,"")
  sub <- final_data_2[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 20) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model, cen = 8) 
  coef_temp[i,] <- coef(cr_temp)
  vcov_temp[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model, cen = 5) 
  coef_preci[i,] <- coef(cr_preci)
  vcov_preci[[i]] <- vcov(cr_preci)
}


# Multivariate mixed-effects meta-analysis - Temperature
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound <- rowMeans(sapply(dlist, function(x) range(x$tmean)))
xvar <- seq(bound[1],bound[2],by=0.1)
argvar=list(fun="ns", knots=quantile(xvar,prob=c(.33,.66)))
bvar <- do.call(onebasis, c(list(x=xvar), argvar))

pred.pool <- crosspred(bvar, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=5)
pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar, coef=coef[i,], vcov=vcov[[i]], 
                                         model.link="log", cen=5))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,10), lwd=2,
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_temp <- mixmeta(coef_temp~1, vcov_temp, method="ml")
summary(mv_temp)

# Plot the pooled exposure-response curve
bound <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar <- seq(bound[1],bound[2],by=0.1)
argvar=list(fun="ns", knots=quantile(xvar,prob=c(.33,.66)))
bvar <- do.call(onebasis, c(list(x=xvar), argvar))

pred.pool <- crosspred(bvar, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=8)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar, coef=coef[i,], vcov=vcov[[i]], 
                                         model.link="log", cen=8))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2,
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)












######################### Zero inflated poisson으로도 한 번 해보기.
regions_1 <- names(final_data_1)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_1))


par(mfrow=c(2,5))
for(i in seq(regions_1)){
  cat(i,"")
  sub <- final_data_1[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        group = sub$s_group,
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         group = sub$s_group,
                         arglag=list(fun="ns", knots=c(3,5)))
  
  model <- zeroinfl(cases ~ cb.temp + cb.preci + ns(week, df = 60) + ns(year, df = 3) | cb.temp, data = sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
                    dist = "poisson")
  
  cen_temp = round(pct_temp[2]) # reference temperature at 33th percentile   
  cen_preci = round(pct_preci[2])
  pred.temp <- crosspred(cb.temp, coef = model$coefficients$count, vcov = vcov(model), cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model$coefficients$count, vcov = vcov(model), cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_1[i])
  plot(pred.preci, "overall", ylim=c(0.9,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_1[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(pct_temp[4]) # 90th pct
  target_preci <- as.character(pct_preci[5]) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}

cbind(logRR, logRRse)


# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions_1, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on mortality (5 vs. 18)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_1, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on mortality (90th vs. 99th)")





















