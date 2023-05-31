############## Second stage
######## 전체 지역 대상으로 
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
library(readr)

## rda 파일은 그냥 load 써야함. 화살표로 객체에 할당하면 안됨.
load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

holi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi.csv")

covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
covariates <- covairate[,-1]

covariates$ndvi_summerGP <-factor(covariates$ndvi_summerGP, 
                                  levels = c(1, 2, 3, 4))
covariates$ndvi_winterGP <-factor(covariates$ndvi_winterGP, 
                                  levels = c(1, 2, 3, 4))
covariates$ndvi_totalGP <-factor(covariates$ndvi_totalGP, 
                                 levels = c(1, 2, 3, 4))
covariates$avgtmeanGP <-factor(covariates$avgtmeanGP, 
                               levels = c(1, 2, 3, 4))
covariates$species <-factor(covariates$species, 
                            levels = c(1, 2))

######################################## final data 조정해서 분석 가능한 형태로 바꾸기
### 2019년까지만 사용. holi 추가.  autocorrelation 항 추가

for (i in 1:length(final_data)){
  final_data[[i]] <- final_data[[i]] %>% filter(year < 2020)
  final_data[[i]]$holi <- holi$holi
  final_data[[i]]$auto <- log(lag(final_data[[i]]$cases, 1) + 0.01)
  final_data[[i]][is.na(final_data[[i]])] <- log(0.01)
}


#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### overall
regions <- names(final_data)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions))



for(i in seq(regions)){
  cat(i,"")
  sub <- final_data[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  

  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}



# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on tsutsugamushi (33th vs. 66th)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (10mm vs. 50mm)")


## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis

ndvi_vector_summer <- covariates %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer <- ndvi_vector_summer[,"ndvi_summer"]

ndvi_vector_winter <- covariates %>% 
  dplyr::select(ndvi_winter)

ndvi_vector_winter <- ndvi_vector_winter[,"ndvi_winter"]

ndvi_vector_total <- covariates %>% 
  dplyr::select(ndvi_total)

ndvi_vector_total <- ndvi_vector_total[,"ndvi_total"]

avgtmean <- covariates %>% 
  dplyr::select(avgtmean)

avgtmean <- avgtmean[,"avgtmean"]

species <- covariates %>% 
  dplyr::select(species)

species <- species[,"species"]


###################### Preci

### Summer
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_summer)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_summer, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), lty="dashed")
abline(h=1, lty = "dotted")


### Winter
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_winter)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_winter, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), lty="dashed")
abline(h=1, lty = "dotted")


### Total
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_total)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_total, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), lty="dashed")
abline(h=1, lty = "dotted")


### Average temperature
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=avgtmean)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(avgtmean, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="Avgtemp", ylab="log(Relative Risk)")

lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), log(preds$pred))
lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), lty="dashed")
abline(h=1, lty = "dotted")


### Species distribution
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=species)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(species, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="species", ylab="log(Relative Risk)")

lines(c(seq(from = min(species), to = max(species), length.out = length(species))), log(preds$pred))
lines(c(seq(from = min(species), to = max(species), length.out = length(species))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(species), to = max(species), length.out = length(species))), lty="dashed")
abline(h=1, lty = "dotted")


###################### Temp 
### Summer
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_summer)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_summer, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_summer), to = max(ndvi_vector_summer), length.out = length(ndvi_vector_summer))), lty="dashed")
abline(h=0, lty = "dotted")


### Winter
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_winter)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_winter, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_winter), to = max(ndvi_vector_winter), length.out = length(ndvi_vector_winter))), lty="dashed")
abline(h=0, lty = "dotted")


### Total
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_total)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_total, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), log(preds$pred))
lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(ndvi_vector_total), to = max(ndvi_vector_total), length.out = length(ndvi_vector_total))), lty="dashed")
abline(h=0, lty = "dotted")


### Average temperature
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=avgtmean)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(avgtmean, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="Avgtemp", ylab="log(Relative Risk)")

lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), log(preds$pred))
lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(avgtmean), to = max(avgtmean), length.out = length(avgtmean))), lty="dashed")
abline(h=0, lty = "dotted")


### Species distribution
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=species)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(species, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="species", ylab="log(Relative Risk)")

lines(c(seq(from = min(species), to = max(species), length.out = length(species))), log(preds$pred))
lines(c(seq(from = min(species), to = max(species), length.out = length(species))), log(preds$ci.lb), lty="dashed")
lines(c(seq(from = min(species), to = max(species), length.out = length(species))), lty="dashed")
abline(h=0, lty = "dotted")



################################################################################
# Second-stage modelling - Multivariate meta-analysis
################################################################################

coef_temp <- matrix(NA, nrow=length(regions), ncol=length(varknot_temp)+1,
                    dimnames=list(regions, paste0("b",seq(3))))
vcov_temp <- vector("list", length(regions))

coef_preci <- matrix(NA, nrow=length(regions), ncol=length(varknot_preci)+1,
                     dimnames=list(regions, paste0("b",seq(3))))
vcov_preci <- vector("list", length(regions))


for(i in seq(regions)){
  cat(i,"")
  sub <- final_data[[i]]
  
  pct_temp <- quantile(sub$tmean, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_temp <- pct_temp[c(3,4)] # knots at 33rd and 66th
  pct_preci <- quantile(sub$total_preci, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  cb.temp <- crossbasis(sub$tmean, lag=8, argvar=list(fun="ns", knots=varknot_temp), 
                        arglag=list(fun="ns", knots=c(3,5))) 
  cb.preci <- crossbasis(sub$total_preci, lag=8, argvar=list(fun="ns", knots=varknot_preci), 
                         arglag=list(fun="ns", knots=c(3,5)))
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp[i,] <- coef(cr_temp)
  vcov_temp[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci[i,] <- coef(cr_preci)
  vcov_preci[[i]] <- vcov(cr_preci)
}


# Multivariate mixed-effects meta-analysis - Temperature
mv_temp <- mixmeta(coef_temp~1, vcov_temp, method="ml")
summary(mv_temp)

# Plot the pooled exposure-response curve
bound_temp <- rowMeans(sapply(final_data, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

pred.pool <- crosspred(bvar_temp, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar_temp, coef=coef_temp[i,], vcov=vcov_temp[[i]], 
                                         model.link="log", cen=10))
par(mfrow=c(1,1))
plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,10), lwd=2, col = "blue", 
     xlab="Temperature (C)", main="Pooled RR")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)

### CI도 음영으로 표현되는 그래프
plot(pred.pool, type="l", ci.arg=list(density=30,col="blue"), ylab="RR", ylim=c(0,15), 
     col="blue", lwd=3, xlab="Temperature (C)", main="pooled RR")




# Multivariate mixed-effects meta-analysis - Precipitation
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound_preci <- rowMeans(sapply(final_data, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

pred.pool <- crosspred(bvar_preci, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=3)
pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar_preci, coef=coef_preci[i,], vcov=vcov_preci[[i]], 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,2), xlim = c(0,250) ,lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled RR")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)

### CI도 음영으로 표현되는 그래프
plot(pred.pool, type="l", ci.arg=list(density=30,col="blue"), ylab="RR", ylim=c(0,3), xlim = c(0,250),
     col="blue", lwd=3, xlab="Precipitation (mm)", main="pooled RR")

?mixmeta

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp <- mixmeta(coef_temp~ndvi_vector_summer + avgtmean + species, vcov_temp, method="ml")
print(summary(mix_temp), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp <- rowMeans(sapply(final_data, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew <- data.frame(ndvi_vector_summer=covariates$ndvi_summer,
                      avgtmean=covariates$avgtmean,
                      species=covariates$species) 

predmix_temp <- predict(mix_temp, datanew,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions)) {
  fit <- fit + predmix_temp[[i]]$fit
  count <- count + 1
}
mix_temp_fit <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions)) {
  vcov <- vcov + predmix_temp[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov <- (vcov / count)


pred.pool <- crosspred(bvar_temp, coef=mix_temp_fit, vcov=mix_temp_vcov, model.link="log", 
                       by=0.1, cen=10)

pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar_temp, coef=predmix_temp[[i]]$fit, vcov=predmix_temp[[i]]$vcov, 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2, col = "blue",
     xlab="Temperature (C)", main="Pooled RR of mixmeta")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

plot(pred.pool, type="l", ci.arg=list(density=30,col="blue"), ylab="RR", ylim=c(0,20), 
     col="blue", lwd=3, xlab="Temperature (C)", main="Pooled RR of mixmeta")







############################# Precipitation
mix_preci <- mixmeta(coef_preci~ndvi_vector_summer + avgtmean + species, vcov_preci, method="ml")
print(summary(mix_preci), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci <- rowMeans(sapply(final_data, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew <- data.frame(ndvi_vector_summer=covariates$ndvi_summer,
                      avgtmean=covariates$avgtmean,
                      species=covariates$species) 

predmix_preci <- predict(mix_preci, datanew,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions)) {
  fit <- fit + predmix_preci[[i]]$fit
  count <- count + 1
}
mix_preci_fit <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions)) {
  vcov <- vcov + predmix_preci[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov <- (vcov / count)


pred.pool <- crosspred(bvar_preci, coef=mix_preci_fit, vcov=mix_preci_vcov, model.link="log", 
                       by=0.1, cen=3)

pred.reg <- lapply(seq(regions), 
                   function(i) crosspred(bvar_preci, coef=predmix_preci[[i]]$fit, vcov=predmix_preci[[i]]$vcov, 
                                         model.link="log", cen=3))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,5), xlim = c(0,250), lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled RR of mixmeta")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

plot(pred.pool, type="l", ci.arg=list(density=30,col="blue"), ylab="RR", ylim=c(0,5), xlim = c(0,250), 
     col="blue", lwd=3, xlab="Precipitation (mm)", main="Pooled RR of mixmeta")



# ################################################################################
# # WALD TEST
# fwald <- function(model,var) {
#   ind <- grep(var,names(coef(model)))
#   coef <- coef(model)[ind]
#   vcov <- vcov(model)[ind,ind]
#   waldstat <- coef%*%solve(vcov)%*%coef
#   df <- length(coef)
#   return(1-pchisq(waldstat,df))
# }
# 
# # TEST THE EFFECTS
# fwald(mvall,"country")
# fwald(mvall,"warmtemp")
# fwald(mvall,"rangetmean")
# 
# ################################################################################
# 
# # INTERCEPT ONLY
# mv.m <- mvmeta(coefall~1,vcovall,cities,control=list(showiter=T))
# summary(mv.m)
# 
# # SINGLE PREDICTOR
# mv.m1 <- mvmeta(coefall~country,vcovall,cities,control=list(showiter=T))
# summary(mv.m1)
# fwald(mv.m1,"country")
# 
# mv.m2 <- mvmeta(coefall~warmtemp,vcovall,cities,control=list(showiter=T))
# summary(mv.m2)
# fwald(mv.m2,"warmtemp")
# 
# mv.m3 <- mvmeta(coefall~rangetmean,vcovall,cities,control=list(showiter=T))
# summary(mv.m3)
# fwald(mv.m3,"rangetmean")
# 
# mv.m4 <- mvmeta(coefall~warmtemp+rangetmean,vcovall,cities,control=list(showiter=T))
# summary(mv.m4)
# fwald(mv.m4,"warmtemp")
# fwald(mv.m4,"rangetmean")
# 
# #

