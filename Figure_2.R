#### Figure 2. 

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

load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

holi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi.csv")

covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
covariates <- covariates[,-1]

summary(covariates)

covariates <- covariates %>% 
  mutate(ndvi_summerGP = ifelse(ndvi_summer < 7000, 1, 
                                ifelse(ndvi_summer  < 7500, 2, 3)),
         avgtmeanGP = ifelse(avgtmean < 12.75, 1, 2),
         avgpreciGP = ifelse(avgpreci < 25.00, 1, 2))

for (i in 1:length(final_data)){
  final_data[[i]] <- final_data[[i]] %>% filter(year < 2020)
  final_data[[i]]$holi <- holi$holi
  final_data[[i]]$auto <- log(lag(final_data[[i]]$cases, 1) + 0.01)
  final_data[[i]][is.na(final_data[[i]])] <- log(0.01)
}

covariates_1 <- covariates %>% 
  dplyr::filter(species == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(species == 2)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']


final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]



################################################################################
# Second-stage modelling - Multivariate meta-analysis
################################################################################
regions <- names(final_data)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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

#####################################################################################

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



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis

ndvi_vector_summer <- covariates %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer <- ndvi_vector_summer[,"ndvi_summer"]


avgtmean <- covariates %>% 
  dplyr::select(avgtmean)

avgtmean <- avgtmean[,"avgtmean"]

avgpreci <- covariates %>% 
  dplyr::select(avgpreci)

avgpreci <- avgpreci[,"avgpreci"]

species <- covariates %>% 
  dplyr::select(species)

species <- species[,"species"]

UrbanRural <- covariates %>% 
  dplyr::select(UrbanRural)

UrbanRural <- UrbanRural[,"UrbanRural"]


## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp <- mixmeta(coef_temp~ndvi_vector_summer + avgtmean + species +UrbanRural +avgpreci, vcov_temp, method="ml")
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
                      species=covariates$species,
                      UrbanRural=covariates$UrbanRural,
                      avgtpreci=covariates$avgpreci) 

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

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,15), lwd=2, col = "#FF9933",
     xlab="Temperature (C)", main="Pooled RR of mixmeta")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

plot(pred.pool, type="l", ci.arg=list(density=30,col="#FF8000"), ylab="RR", ylim=c(0,20), 
     col="#FF8000", lwd=3, xlab="Temperature (C)", main="Pooled RR of mixmeta")







############################# Precipitation
mix_preci <- mixmeta(coef_preci~ndvi_vector_summer + avgtmean + species  +UrbanRural +avgpreci, vcov_preci, method="ml")
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
                      species=covariates$species,
                      UrbanRural=covariates$UrbanRural,
                      avgtpreci=covariates$avgpreci)  

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

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,5), xlim = c(0,250), lwd=2, col = "#0080FF",
     xlab="Precipitation (mm)", main="Pooled RR of mixmeta")
for(i in seq(regions)) lines(pred.reg[[i]], col="grey")

plot(pred.pool, type="l", ci.arg=list(density=30,col="#0080FF"), ylab="RR", ylim=c(0,5), xlim = c(0,250), 
     col="#0080FF", lwd=3, xlab="Precipitation (mm)", main="Pooled RR of mixmeta")









#######################################################################################
########### Species 나눠서 그래프 그리고 한 그래프에 겹쳐서 나타내기
#######################################################################################


regions_1 <- names(final_data_1)
logRR_temp_1 <- logRRse_temp_1 <- logRR_preci_1 <- logRRse_preci_1 <- vector("numeric",length(regions_1))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_1[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_1[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_1[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_1[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
                    dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp_1 <- vector("list", length(regions_1))

coef_preci_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                     dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci_1 <- vector("list", length(regions_1))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_1[i,] <- coef(cr_temp)
  vcov_temp_1[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_1[i,] <- coef(cr_preci)
  vcov_preci_1[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_1 <- covariates_1 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_1 <- ndvi_vector_summer_1[,"ndvi_summer"]

avgtmean_1 <- covariates_1 %>% 
  dplyr::select(avgtmean)

avgtmean_1 <- avgtmean_1[,"avgtmean"]

avgpreci_1 <- covariates_1 %>% 
  dplyr::select(avgpreci)

avgpreci_1 <- avgpreci_1[,"avgpreci"]

UrbanRural_1 <- covariates_1 %>% 
  dplyr::select(UrbanRural)

UrbanRural_1 <- UrbanRural_1[,"UrbanRural"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_1 <- mixmeta(coef_temp_1~ndvi_vector_summer_1 + avgtmean_1 + UrbanRural_1 + avgpreci_1, vcov_temp_1, method="ml") # speceis는 제외.
print(summary(mix_temp_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_1 <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp_1 <- seq(bound_temp_1[1],bound_temp_1[2],by=0.1)
argvar_temp_1=list(fun="ns", knots=quantile(xvar_temp_1,prob=c(.33,.66)))
bvar_temp_1 <- do.call(onebasis, c(list(x=xvar_temp_1), argvar_temp_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_1 <- data.frame(ndvi_vector_summer_1=covariates_1$ndvi_summer,
                      avgtmean_1=covariates_1$avgtmean,
                      UrbanRural_1=covariates_1$UrbanRural,
                      avgpreci_1=covariates_1$avgpreci
                      ) 

predmix_temp_1 <- predict(mix_temp_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_temp_1[[i]]$fit
  count <- count + 1
}
mix_temp_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_temp_1[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_1 <- (vcov / count)


pred.pool_temp_1 <- crosspred(bvar_temp_1, coef=mix_temp_fit_1, vcov=mix_temp_vcov_1, model.link="log", 
                       by=0.1, cen=10)

pred.reg_temp_1 <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_temp_1, coef=predmix_temp_1[[i]]$fit, vcov=predmix_temp_1[[i]]$vcov, 
                                         model.link="log", cen=10))




############################# Precipitation
mix_preci_1 <- mixmeta(coef_preci_1~ndvi_vector_summer_1 + avgtmean_1 + UrbanRural_1 + avgpreci_1, vcov_preci_1, method="ml")
print(summary(mix_preci_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_1 <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci_1 <- seq(bound_preci_1[1],bound_preci_1[2],by=0.1)
argvar_preci_1=list(fun="ns", knots=quantile(xvar_preci_1,prob=c(.33,.66)))
bvar_preci_1 <- do.call(onebasis, c(list(x=xvar_preci_1), argvar_preci_1))


predmix_preci_1 <- predict(mix_preci_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_preci_1[[i]]$fit
  count <- count + 1
}
mix_preci_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_preci_1[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_1 <- (vcov / count)


pred.pool_preci_1 <- crosspred(bvar_preci_1, coef=mix_preci_fit_1, vcov=mix_preci_vcov_1, model.link="log", 
                       by=0.1, cen=3)

pred.reg_preci_1 <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_preci_1, coef=predmix_preci_1[[i]]$fit, vcov=predmix_preci_1[[i]]$vcov, 
                                         model.link="log", cen=3))



#############################################################################



regions_2 <- names(final_data_2)
logRR_temp_2 <- logRRse_temp_2 <- logRR_preci_2 <- logRRse_preci_2 <- vector("numeric",length(regions_2))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_2[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_2[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_2[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_2[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp_2 <- vector("list", length(regions_2))

coef_preci_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci_2 <- vector("list", length(regions_2))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_2[i,] <- coef(cr_temp)
  vcov_temp_2[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_2[i,] <- coef(cr_preci)
  vcov_preci_2[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_2 <- covariates_2 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_2 <- ndvi_vector_summer_2[,"ndvi_summer"]

avgtmean_2 <- covariates_2 %>% 
  dplyr::select(avgtmean)

avgtmean_2 <- avgtmean_2[,"avgtmean"]

avgpreci_2 <- covariates_2 %>% 
  dplyr::select(avgpreci)

avgpreci_2 <- avgpreci_2[,"avgpreci"]

UrbanRural_2 <- covariates_2 %>% 
  dplyr::select(UrbanRural)

UrbanRural_2 <- UrbanRural_2[,"UrbanRural"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_2 <- mixmeta(coef_temp_2~ndvi_vector_summer_2 + avgtmean_2 + UrbanRural_2 + avgpreci_2, vcov_temp_2, method="ml") # speceis는 제외.
print(summary(mix_temp_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_2 <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp_2 <- seq(bound_temp_2[1],bound_temp_2[2],by=0.1)
argvar_temp_2=list(fun="ns", knots=quantile(xvar_temp_2,prob=c(.33,.66)))
bvar_temp_2 <- do.call(onebasis, c(list(x=xvar_temp_2), argvar_temp_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_2 <- data.frame(ndvi_vector_summer_2=covariates_2$ndvi_summer,
                        avgtmean_2=covariates_2$avgtmean,
                        UrbanRural_2=covariates_2$UrbanRural,
                        avgpreci_2=covariates_2$avgpreci) 

predmix_temp_2 <- predict(mix_temp_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_temp_2[[i]]$fit
  count <- count + 1
}
mix_temp_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_temp_2[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_2 <- (vcov / count)


pred.pool_temp_2 <- crosspred(bvar_temp_2, coef=mix_temp_fit_2, vcov=mix_temp_vcov_2, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_2 <- lapply(seq(regions_2), 
                          function(i) crosspred(bvar_temp_2, coef=predmix_temp_2[[i]]$fit, vcov=predmix_temp_2[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_2 <- mixmeta(coef_preci_2~ndvi_vector_summer_2 + avgtmean_2 + UrbanRural_2 + avgpreci_2, vcov_preci_2, method="ml")
print(summary(mix_preci_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_2 <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci_2 <- seq(bound_preci_2[1],bound_preci_2[2],by=0.1)
argvar_preci_2=list(fun="ns", knots=quantile(xvar_preci_2,prob=c(.33,.66)))
bvar_preci_2 <- do.call(onebasis, c(list(x=xvar_preci_2), argvar_preci_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_2 <- predict(mix_preci_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_preci_2[[i]]$fit
  count <- count + 1
}
mix_preci_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_preci_2[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_2 <- (vcov / count)


pred.pool_preci_2 <- crosspred(bvar_preci_2, coef=mix_preci_fit_2, vcov=mix_preci_vcov_2, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_2 <- lapply(seq(regions_2), 
                           function(i) crosspred(bvar_preci_2, coef=predmix_preci_2[[i]]$fit, vcov=predmix_preci_2[[i]]$vcov, 
                                                 model.link="log", cen=3))



############# Species 별로 나눠서 한 그래프에 넣을 때는 CI area를 이렇게 하는게 더 잘 나타날듯

### temp
plot(pred.pool_temp_2, type="l", ci.arg=list("area", col=adjustcolor("#3333FF", alpha = 0.2)), ylab="RR", ylim=c(0,30), 
     col="#0000FF", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_temp_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#FF3333", alpha = 0.2)),  col="#FF0000", lwd=2)



### Preci
plot(pred.pool_preci_2, type="l", ci.arg=list("area", col=adjustcolor("#3333FF", alpha = 0.2)), ylab="RR", ylim=c(0,8), xlim = c(0, 250),
     col="#0000FF", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_preci_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#FF3333", alpha = 0.2)),  col="#FF0000", lwd=2)
















#######################################################################################
########### NDVI 나눠서 그래프 그리고 한 그래프에 겹쳐서 나타내기
#######################################################################################
covariates_1 <- covariates %>% 
  dplyr::filter(ndvi_summerGP == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(ndvi_summerGP == 2)
covariates_3 <- covariates %>% 
  dplyr::filter(ndvi_summerGP == 3)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']
name_3 <- covariates_3[,'NAME_2']


final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]
final_data_3 <- final_data[name_3]


regions_1 <- names(final_data_1)
logRR_temp_1 <- logRRse_temp_1 <- logRR_preci_1 <- logRRse_preci_1 <- vector("numeric",length(regions_1))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_1[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_1[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_1[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_1[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp_1 <- vector("list", length(regions_1))

coef_preci_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci_1 <- vector("list", length(regions_1))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_1[i,] <- coef(cr_temp)
  vcov_temp_1[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_1[i,] <- coef(cr_preci)
  vcov_preci_1[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음음
avgtmean_1 <- covariates_1 %>% 
  dplyr::select(avgtmean)

avgtmean_1 <- avgtmean_1[,"avgtmean"]

species_1 <- covariates_1 %>% 
  dplyr::select(species)

species_1 <- species_1[,"species"]

UrbanRural_1 <- covariates_1 %>% 
  dplyr::select(UrbanRural)

UrbanRural_1 <- UrbanRural_1[,"UrbanRural"]

avgpreci_1 <- covariates_1 %>% 
  dplyr::select(avgpreci)

avgpreci_1 <- avgpreci_1[,"avgpreci"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_1 <- mixmeta(coef_temp_1~ avgtmean_1 + species_1 + UrbanRural_1 + avgpreci_1, vcov_temp_1, method="ml") # NDVI는 제외.
print(summary(mix_temp_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_1 <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp_1 <- seq(bound_temp_1[1],bound_temp_1[2],by=0.1)
argvar_temp_1=list(fun="ns", knots=quantile(xvar_temp_1,prob=c(.33,.66)))
bvar_temp_1 <- do.call(onebasis, c(list(x=xvar_temp_1), argvar_temp_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_1 <- data.frame(avgtmean_1=covariates_1$avgtmean,
                        species_1=covariates_1$species,
                        UrbanRural_1=covariates_1$UrbanRural,
                        avgpreci_1=covariates_1$avgpreci) 

predmix_temp_1 <- predict(mix_temp_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_temp_1[[i]]$fit
  count <- count + 1
}
mix_temp_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_temp_1[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_1 <- (vcov / count)


pred.pool_temp_1 <- crosspred(bvar_temp_1, coef=mix_temp_fit_1, vcov=mix_temp_vcov_1, model.link="log", 
                              by=0.1, cen=10)

pred.reg_temp_1 <- lapply(seq(regions_1), 
                          function(i) crosspred(bvar_temp_1, coef=predmix_temp_1[[i]]$fit, vcov=predmix_temp_1[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_1 <- mixmeta(coef_preci_1~ avgtmean_1 + species_1 + UrbanRural_1 + avgpreci_1, vcov_preci_1, method="ml")
print(summary(mix_preci_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_1 <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci_1 <- seq(bound_preci_1[1],bound_preci_1[2],by=0.1)
argvar_preci_1=list(fun="ns", knots=quantile(xvar_preci_1,prob=c(.33,.66)))
bvar_preci_1 <- do.call(onebasis, c(list(x=xvar_preci_1), argvar_preci_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.


predmix_preci_1 <- predict(mix_preci_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_preci_1[[i]]$fit
  count <- count + 1
}
mix_preci_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_preci_1[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_1 <- (vcov / count)


pred.pool_preci_1 <- crosspred(bvar_preci_1, coef=mix_preci_fit_1, vcov=mix_preci_vcov_1, model.link="log", 
                               by=0.1, cen=3)

pred.reg_preci_1 <- lapply(seq(regions_1), 
                           function(i) crosspred(bvar_preci_1, coef=predmix_preci_1[[i]]$fit, vcov=predmix_preci_1[[i]]$vcov, 
                                                 model.link="log", cen=3))



#############################################################################



regions_2 <- names(final_data_2)
logRR_temp_2 <- logRRse_temp_2 <- logRR_preci_2 <- logRRse_preci_2 <- vector("numeric",length(regions_2))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_2[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_2[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_2[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_2[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp_2 <- vector("list", length(regions_2))

coef_preci_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci_2 <- vector("list", length(regions_2))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_2[i,] <- coef(cr_temp)
  vcov_temp_2[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_2[i,] <- coef(cr_preci)
  vcov_preci_2[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음

avgtmean_2 <- covariates_2 %>% 
  dplyr::select(avgtmean)

avgtmean_2 <- avgtmean_2[,"avgtmean"]

species_2 <- covariates_2 %>% 
  dplyr::select(species)

species_2 <- species_2[,"species"]

avgpreci_2 <- covariates_2 %>% 
  dplyr::select(avgpreci)

avgpreci_2 <- avgpreci_2[,"avgpreci"]

UrbanRural_2 <- covariates_2 %>% 
  dplyr::select(UrbanRural)

UrbanRural_2 <- UrbanRural_2[,"UrbanRural"]


## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_2 <- mixmeta(coef_temp_2~avgtmean_2 + species_2 +UrbanRural_2 + avgpreci_2, vcov_temp_2, method="ml") # speceis는 제외.
print(summary(mix_temp_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_2 <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp_2 <- seq(bound_temp_2[1],bound_temp_2[2],by=0.1)
argvar_temp_2=list(fun="ns", knots=quantile(xvar_temp_2,prob=c(.33,.66)))
bvar_temp_2 <- do.call(onebasis, c(list(x=xvar_temp_2), argvar_temp_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_2 <- data.frame(avgtmean_2=covariates_2$avgtmean,
                        species_2=covariates_2$species,
                        UrbanRural_2=covariates_2$UrbanRural,
                        avgpreci_2=covariates_2$avgpreci) 

predmix_temp_2 <- predict(mix_temp_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_temp_2[[i]]$fit
  count <- count + 1
}
mix_temp_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_temp_2[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_2 <- (vcov / count)


pred.pool_temp_2 <- crosspred(bvar_temp_2, coef=mix_temp_fit_2, vcov=mix_temp_vcov_2, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_2 <- lapply(seq(regions_2), 
                          function(i) crosspred(bvar_temp_2, coef=predmix_temp_2[[i]]$fit, vcov=predmix_temp_2[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_2 <- mixmeta(coef_preci_2~avgtmean_2 + species_2 +UrbanRural_2 + avgpreci_2 , vcov_preci_2, method="ml")
print(summary(mix_preci_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_2 <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci_2 <- seq(bound_preci_2[1],bound_preci_2[2],by=0.1)
argvar_preci_2=list(fun="ns", knots=quantile(xvar_preci_2,prob=c(.33,.66)))
bvar_preci_2 <- do.call(onebasis, c(list(x=xvar_preci_2), argvar_preci_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_2 <- predict(mix_preci_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_preci_2[[i]]$fit
  count <- count + 1
}
mix_preci_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_preci_2[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_2 <- (vcov / count)


pred.pool_preci_2 <- crosspred(bvar_preci_2, coef=mix_preci_fit_2, vcov=mix_preci_vcov_2, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_2 <- lapply(seq(regions_2), 
                           function(i) crosspred(bvar_preci_2, coef=predmix_preci_2[[i]]$fit, vcov=predmix_preci_2[[i]]$vcov, 
                                                 model.link="log", cen=3))



##############################################################################################




regions_3 <- names(final_data_3)
logRR_temp_3 <- logRRse_temp_3 <- logRR_preci_3 <- logRRse_preci_3 <- vector("numeric",length(regions_3))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
for(i in seq(regions_3)){
  cat(i,"")
  sub <- final_data_3[[i]]
  
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
  logRR_temp_3[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_3[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_3[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_3[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_3 <- matrix(NA, nrow=length(regions_3), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_3, paste0("b",seq(3))))
vcov_temp_3 <- vector("list", length(regions_3))

coef_preci_3 <- matrix(NA, nrow=length(regions_3), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_3, paste0("b",seq(3))))
vcov_preci_3 <- vector("list", length(regions_3))


for(i in seq(regions_3)){
  cat(i,"")
  sub <- final_data_3[[i]]
  
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
  coef_temp_3[i,] <- coef(cr_temp)
  vcov_temp_3[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_3[i,] <- coef(cr_preci)
  vcov_preci_3[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음

avgtmean_3 <- covariates_3 %>% 
  dplyr::select(avgtmean)

avgtmean_3 <- avgtmean_3[,"avgtmean"]

species_3 <- covariates_3 %>% 
  dplyr::select(species)

species_3 <- species_3[,"species"]

avgpreci_3 <- covariates_3 %>% 
  dplyr::select(avgpreci)

avgpreci_3 <- avgpreci_3[,"avgpreci"]

UrbanRural_3 <- covariates_3 %>% 
  dplyr::select(UrbanRural)

UrbanRural_3 <- UrbanRural_3[,"UrbanRural"]




## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_3 <- mixmeta(coef_temp_3~avgtmean_3 + species_3 + UrbanRural_3 + avgpreci_3, vcov_temp_3, method="ml") # speceis는 제외.
print(summary(mix_temp_3), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_3 <- rowMeans(sapply(final_data_3, function(x) range(x$tmean)))
xvar_temp_3 <- seq(bound_temp_3[1],bound_temp_3[2],by=0.1)
argvar_temp_3=list(fun="ns", knots=quantile(xvar_temp_3,prob=c(.33,.66)))
bvar_temp_3 <- do.call(onebasis, c(list(x=xvar_temp_3), argvar_temp_3))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_3 <- data.frame(avgtmean_3=covariates_3$avgtmean,
                        species_3=covariates_3$species,
                        UrbanRural_3=covariates_3$UrbanRural,
                        avgpreci_3=covariates_3$avgpreci) 

predmix_temp_3 <- predict(mix_temp_3, datanew_3,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_3)) {
  fit <- fit + predmix_temp_3[[i]]$fit
  count <- count + 1
}
mix_temp_fit_3 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_3)) {
  vcov <- vcov + predmix_temp_3[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_3 <- (vcov / count)


pred.pool_temp_3 <- crosspred(bvar_temp_3, coef=mix_temp_fit_3, vcov=mix_temp_vcov_3, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_3 <- lapply(seq(regions_3), 
                          function(i) crosspred(bvar_temp_3, coef=predmix_temp_3[[i]]$fit, vcov=predmix_temp_3[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_3 <- mixmeta(coef_preci_3~avgtmean_3 + species_3 + UrbanRural_3 + avgpreci_3, vcov_preci_3, method="ml")
print(summary(mix_preci_3), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_3 <- rowMeans(sapply(final_data_3, function(x) range(x$total_preci)))
xvar_preci_3 <- seq(bound_preci_3[1],bound_preci_3[2],by=0.1)
argvar_preci_3=list(fun="ns", knots=quantile(xvar_preci_3,prob=c(.33,.66)))
bvar_preci_3 <- do.call(onebasis, c(list(x=xvar_preci_3), argvar_preci_3))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_3 <- predict(mix_preci_3, datanew_3,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_3)) {
  fit <- fit + predmix_preci_3[[i]]$fit
  count <- count + 1
}
mix_preci_fit_3 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_3)) {
  vcov <- vcov + predmix_preci_3[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_3 <- (vcov / count)


pred.pool_preci_3 <- crosspred(bvar_preci_3, coef=mix_preci_fit_3, vcov=mix_preci_vcov_3, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_3 <- lapply(seq(regions_3), 
                           function(i) crosspred(bvar_preci_3, coef=predmix_preci_3[[i]]$fit, vcov=predmix_preci_3[[i]]$vcov, 
                                                 model.link="log", cen=3))



### temp

plot(pred.pool_temp_3, type="l", ci.arg=list("area", col=adjustcolor("#003300", alpha = 0.2)), ylab="RR", ylim=c(0,20), 
     col="#003300", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_temp_2,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#00CC00", alpha = 0.2)),  col="#00CC00", lwd=2)
lines(pred.pool_temp_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#FF9933", alpha = 0.2)),  col="#FF8000", lwd=2)


### Preci
plot(pred.pool_preci_1, type="l", ci.arg=list("area", col=adjustcolor("#FF9933", alpha = 0.15)), ylab="RR", ylim=c(0,5), xlim = c(0, 250),
     col="#FF8000", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_preci_2,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#80FF00", alpha = 0.15)),  col="#00CC00", lwd=2)
lines(pred.pool_preci_3,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#003300", alpha = 0.3)),  col="#003300", lwd=2)






















covariates_1 <- covariates %>% 
  dplyr::filter(avgpreciGP == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(avgpreciGP == 2)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']

final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]



#######################################################################################
########### Precipitation 나눠서 그래프 그리고 한 그래프에 겹쳐서 나타내기
#######################################################################################


regions_1 <- names(final_data_1)
logRR_temp_1 <- logRRse_temp_1 <- logRR_preci_1 <- logRRse_preci_1 <- vector("numeric",length(regions_1))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_1[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_1[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_1[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_1[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp_1 <- vector("list", length(regions_1))

coef_preci_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci_1 <- vector("list", length(regions_1))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_1[i,] <- coef(cr_temp)
  vcov_temp_1[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_1[i,] <- coef(cr_preci)
  vcov_preci_1[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_1 <- covariates_1 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_1 <- ndvi_vector_summer_1[,"ndvi_summer"]

avgtmean_1 <- covariates_1 %>% 
  dplyr::select(avgtmean)

avgtmean_1 <- avgtmean_1[,"avgtmean"]

species_1 <- covariates_1 %>% 
  dplyr::select(species)

species_1 <- species_1[,"species"]

UrbanRural_1 <- covariates_1 %>% 
  dplyr::select(UrbanRural)

UrbanRural_1 <- UrbanRural_1[,"UrbanRural"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_1 <- mixmeta(coef_temp_1~ndvi_vector_summer_1 +  avgtmean_1 + UrbanRural_1 + species_1, vcov_temp_1, method="ml") # speceis는 제외.
print(summary(mix_temp_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_1 <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp_1 <- seq(bound_temp_1[1],bound_temp_1[2],by=0.1)
argvar_temp_1=list(fun="ns", knots=quantile(xvar_temp_1,prob=c(.33,.66)))
bvar_temp_1 <- do.call(onebasis, c(list(x=xvar_temp_1), argvar_temp_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_1 <- data.frame(ndvi_vector_summer_1=covariates_1$ndvi_summer,
                        avgtmean_1=covariates_1$avgtmean,
                        UrbanRural_1=covariates_1$UrbanRural,
                        species_1=covariates_1$species
) 

predmix_temp_1 <- predict(mix_temp_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_temp_1[[i]]$fit
  count <- count + 1
}
mix_temp_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_temp_1[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_1 <- (vcov / count)


pred.pool_temp_1 <- crosspred(bvar_temp_1, coef=mix_temp_fit_1, vcov=mix_temp_vcov_1, model.link="log", 
                              by=0.1, cen=10)

pred.reg_temp_1 <- lapply(seq(regions_1), 
                          function(i) crosspred(bvar_temp_1, coef=predmix_temp_1[[i]]$fit, vcov=predmix_temp_1[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_1 <- mixmeta(coef_preci_1~ndvi_vector_summer_1 + avgtmean_1 + UrbanRural_1 + species_1, vcov_preci_1, method="ml")
print(summary(mix_preci_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_1 <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci_1 <- seq(bound_preci_1[1],bound_preci_1[2],by=0.1)
argvar_preci_1=list(fun="ns", knots=quantile(xvar_preci_1,prob=c(.33,.66)))
bvar_preci_1 <- do.call(onebasis, c(list(x=xvar_preci_1), argvar_preci_1))


predmix_preci_1 <- predict(mix_preci_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_preci_1[[i]]$fit
  count <- count + 1
}
mix_preci_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_preci_1[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_1 <- (vcov / count)


pred.pool_preci_1 <- crosspred(bvar_preci_1, coef=mix_preci_fit_1, vcov=mix_preci_vcov_1, model.link="log", 
                               by=0.1, cen=3)

pred.reg_preci_1 <- lapply(seq(regions_1), 
                           function(i) crosspred(bvar_preci_1, coef=predmix_preci_1[[i]]$fit, vcov=predmix_preci_1[[i]]$vcov, 
                                                 model.link="log", cen=3))



#############################################################################



regions_2 <- names(final_data_2)
logRR_temp_2 <- logRRse_temp_2 <- logRR_preci_2 <- logRRse_preci_2 <- vector("numeric",length(regions_2))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_2[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_2[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_2[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_2[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp_2 <- vector("list", length(regions_2))

coef_preci_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci_2 <- vector("list", length(regions_2))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_2[i,] <- coef(cr_temp)
  vcov_temp_2[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_2[i,] <- coef(cr_preci)
  vcov_preci_2[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_2 <- covariates_2 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_2 <- ndvi_vector_summer_2[,"ndvi_summer"]

avgtmean_2 <- covariates_2 %>% 
  dplyr::select(avgtmean)

avgtmean_2 <- avgtmean_2[,"avgtmean"]

species_2 <- covariates_2 %>% 
  dplyr::select(species)

species_2 <- species_2[,"species"]

UrbanRural_2 <- covariates_2 %>% 
  dplyr::select(UrbanRural)

UrbanRural_2 <- UrbanRural_2[,"UrbanRural"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_2 <- mixmeta(coef_temp_2~ndvi_vector_summer_2 + avgtmean_2 + UrbanRural_2 + species_2, vcov_temp_2, method="ml") # speceis는 제외.
print(summary(mix_temp_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_2 <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp_2 <- seq(bound_temp_2[1],bound_temp_2[2],by=0.1)
argvar_temp_2=list(fun="ns", knots=quantile(xvar_temp_2,prob=c(.33,.66)))
bvar_temp_2 <- do.call(onebasis, c(list(x=xvar_temp_2), argvar_temp_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_2 <- data.frame(ndvi_vector_summer_2=covariates_2$ndvi_summer,
                        avgtmean_2=covariates_2$avgtmean,
                        UrbanRural_2=covariates_2$UrbanRural,
                        species_2=covariates_2$species) 

predmix_temp_2 <- predict(mix_temp_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_temp_2[[i]]$fit
  count <- count + 1
}
mix_temp_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_temp_2[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_2 <- (vcov / count)


pred.pool_temp_2 <- crosspred(bvar_temp_2, coef=mix_temp_fit_2, vcov=mix_temp_vcov_2, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_2 <- lapply(seq(regions_2), 
                          function(i) crosspred(bvar_temp_2, coef=predmix_temp_2[[i]]$fit, vcov=predmix_temp_2[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_2 <- mixmeta(coef_preci_2~ndvi_vector_summer_2 + avgtmean_2 + UrbanRural_2 + species_2, vcov_preci_2, method="ml")
print(summary(mix_preci_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_2 <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci_2 <- seq(bound_preci_2[1],bound_preci_2[2],by=0.1)
argvar_preci_2=list(fun="ns", knots=quantile(xvar_preci_2,prob=c(.33,.66)))
bvar_preci_2 <- do.call(onebasis, c(list(x=xvar_preci_2), argvar_preci_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_2 <- predict(mix_preci_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_preci_2[[i]]$fit
  count <- count + 1
}
mix_preci_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_preci_2[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_2 <- (vcov / count)


pred.pool_preci_2 <- crosspred(bvar_preci_2, coef=mix_preci_fit_2, vcov=mix_preci_vcov_2, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_2 <- lapply(seq(regions_2), 
                           function(i) crosspred(bvar_preci_2, coef=predmix_preci_2[[i]]$fit, vcov=predmix_preci_2[[i]]$vcov, 
                                                 model.link="log", cen=3))



############# Species 별로 나눠서 한 그래프에 넣을 때는 CI area를 이렇게 하는게 더 잘 나타날듯

### temp
plot(pred.pool_temp_2, type="l", ci.arg=list("area", col=adjustcolor("#0080FF", alpha = 0.2)), ylab="RR", ylim=c(0,30), 
     col="#0080FF", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_temp_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#FF8000", alpha = 0.2)),  col="#FF8000", lwd=2)


### Preci
plot(pred.pool_preci_2, type="l", ci.arg=list("area", col=adjustcolor("#0080FF", alpha = 0.2)), ylab="RR", ylim=c(0,8), xlim = c(0, 250),
     col="#0080FF", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_preci_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#FF8000", alpha = 0.2)),  col="#FF8000", lwd=2)















covariates_1 <- covariates %>% 
  dplyr::filter(UrbanRural == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(UrbanRural == 2)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']

final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]



#######################################################################################
########### Urban Rural 나눠서 그래프 그리고 한 그래프에 겹쳐서 나타내기
#######################################################################################


regions_1 <- names(final_data_1)
logRR_temp_1 <- logRRse_temp_1 <- logRR_preci_1 <- logRRse_preci_1 <- vector("numeric",length(regions_1))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_1[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_1[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_1[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_1[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp_1 <- vector("list", length(regions_1))

coef_preci_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci_1 <- vector("list", length(regions_1))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_1[i,] <- coef(cr_temp)
  vcov_temp_1[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_1[i,] <- coef(cr_preci)
  vcov_preci_1[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_1 <- covariates_1 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_1 <- ndvi_vector_summer_1[,"ndvi_summer"]

avgtmean_1 <- covariates_1 %>% 
  dplyr::select(avgtmean)

avgtmean_1 <- avgtmean_1[,"avgtmean"]

species_1 <- covariates_1 %>% 
  dplyr::select(species)

species_1 <- species_1[,"species"]

avgpreci_1 <- covariates_1 %>% 
  dplyr::select(avgpreci)

avgpreci_1 <- avgpreci_1[,"avgpreci"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_1 <- mixmeta(coef_temp_1~ndvi_vector_summer_1 +  avgtmean_1 + avgpreci_1 + species_1, vcov_temp_1, method="ml") # speceis는 제외.
print(summary(mix_temp_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_1 <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp_1 <- seq(bound_temp_1[1],bound_temp_1[2],by=0.1)
argvar_temp_1=list(fun="ns", knots=quantile(xvar_temp_1,prob=c(.33,.66)))
bvar_temp_1 <- do.call(onebasis, c(list(x=xvar_temp_1), argvar_temp_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_1 <- data.frame(ndvi_vector_summer_1=covariates_1$ndvi_summer,
                        avgtmean_1=covariates_1$avgtmean,
                        avgpreci_1=covariates_1$UrbanRural,
                        species_1=covariates_1$species
) 

predmix_temp_1 <- predict(mix_temp_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_temp_1[[i]]$fit
  count <- count + 1
}
mix_temp_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_temp_1[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_1 <- (vcov / count)


pred.pool_temp_1 <- crosspred(bvar_temp_1, coef=mix_temp_fit_1, vcov=mix_temp_vcov_1, model.link="log", 
                              by=0.1, cen=10)

pred.reg_temp_1 <- lapply(seq(regions_1), 
                          function(i) crosspred(bvar_temp_1, coef=predmix_temp_1[[i]]$fit, vcov=predmix_temp_1[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_1 <- mixmeta(coef_preci_1~ndvi_vector_summer_1 + avgtmean_1 + avgpreci_1 + species_1, vcov_preci_1, method="ml")
print(summary(mix_preci_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_1 <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci_1 <- seq(bound_preci_1[1],bound_preci_1[2],by=0.1)
argvar_preci_1=list(fun="ns", knots=quantile(xvar_preci_1,prob=c(.33,.66)))
bvar_preci_1 <- do.call(onebasis, c(list(x=xvar_preci_1), argvar_preci_1))


predmix_preci_1 <- predict(mix_preci_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_preci_1[[i]]$fit
  count <- count + 1
}
mix_preci_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_preci_1[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_1 <- (vcov / count)


pred.pool_preci_1 <- crosspred(bvar_preci_1, coef=mix_preci_fit_1, vcov=mix_preci_vcov_1, model.link="log", 
                               by=0.1, cen=3)

pred.reg_preci_1 <- lapply(seq(regions_1), 
                           function(i) crosspred(bvar_preci_1, coef=predmix_preci_1[[i]]$fit, vcov=predmix_preci_1[[i]]$vcov, 
                                                 model.link="log", cen=3))



#############################################################################



regions_2 <- names(final_data_2)
logRR_temp_2 <- logRRse_temp_2 <- logRR_preci_2 <- logRRse_preci_2 <- vector("numeric",length(regions_2))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_2[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_2[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_2[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_2[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp_2 <- vector("list", length(regions_2))

coef_preci_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci_2 <- vector("list", length(regions_2))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_2[i,] <- coef(cr_temp)
  vcov_temp_2[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_2[i,] <- coef(cr_preci)
  vcov_preci_2[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_2 <- covariates_2 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_2 <- ndvi_vector_summer_2[,"ndvi_summer"]

avgtmean_2 <- covariates_2 %>% 
  dplyr::select(avgtmean)

avgtmean_2 <- avgtmean_2[,"avgtmean"]

species_2 <- covariates_2 %>% 
  dplyr::select(species)

species_2 <- species_2[,"species"]

avgpreci_2 <- covariates_2 %>% 
  dplyr::select(avgpreci)

avgpreci_2 <- avgpreci_2[,"avgpreci"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_2 <- mixmeta(coef_temp_2~ndvi_vector_summer_2 + avgtmean_2 + avgpreci_2 + species_2, vcov_temp_2, method="ml") # speceis는 제외.
print(summary(mix_temp_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_2 <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp_2 <- seq(bound_temp_2[1],bound_temp_2[2],by=0.1)
argvar_temp_2=list(fun="ns", knots=quantile(xvar_temp_2,prob=c(.33,.66)))
bvar_temp_2 <- do.call(onebasis, c(list(x=xvar_temp_2), argvar_temp_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_2 <- data.frame(ndvi_vector_summer_2=covariates_2$ndvi_summer,
                        avgtmean_2=covariates_2$avgtmean,
                        avgpreci_2=covariates_2$UrbanRural,
                        species_2=covariates_2$species) 

predmix_temp_2 <- predict(mix_temp_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_temp_2[[i]]$fit
  count <- count + 1
}
mix_temp_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_temp_2[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_2 <- (vcov / count)


pred.pool_temp_2 <- crosspred(bvar_temp_2, coef=mix_temp_fit_2, vcov=mix_temp_vcov_2, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_2 <- lapply(seq(regions_2), 
                          function(i) crosspred(bvar_temp_2, coef=predmix_temp_2[[i]]$fit, vcov=predmix_temp_2[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_2 <- mixmeta(coef_preci_2~ndvi_vector_summer_2 + avgtmean_2 + avgpreci_2 + species_2, vcov_preci_2, method="ml")
print(summary(mix_preci_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_2 <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci_2 <- seq(bound_preci_2[1],bound_preci_2[2],by=0.1)
argvar_preci_2=list(fun="ns", knots=quantile(xvar_preci_2,prob=c(.33,.66)))
bvar_preci_2 <- do.call(onebasis, c(list(x=xvar_preci_2), argvar_preci_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_2 <- predict(mix_preci_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_preci_2[[i]]$fit
  count <- count + 1
}
mix_preci_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_preci_2[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_2 <- (vcov / count)


pred.pool_preci_2 <- crosspred(bvar_preci_2, coef=mix_preci_fit_2, vcov=mix_preci_vcov_2, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_2 <- lapply(seq(regions_2), 
                           function(i) crosspred(bvar_preci_2, coef=predmix_preci_2[[i]]$fit, vcov=predmix_preci_2[[i]]$vcov, 
                                                 model.link="log", cen=3))



############# Species 별로 나눠서 한 그래프에 넣을 때는 CI area를 이렇게 하는게 더 잘 나타날듯

### temp
plot(pred.pool_temp_2, type="l", ci.arg=list("area", col=adjustcolor("#990099", alpha = 0.2)), ylab="RR", ylim=c(0,30), 
     col="#990099", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_temp_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#66CC00", alpha = 0.2)),  col="#66CC00", lwd=2)


### Preci
plot(pred.pool_preci_2, type="l", ci.arg=list("area", col=adjustcolor("#990099", alpha = 0.2)), ylab="RR", ylim=c(0,8), xlim = c(0, 250),
     col="#990099", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_preci_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#66CC00", alpha = 0.2)),  col="#66CC00", lwd=2)












covariates_1 <- covariates %>% 
  dplyr::filter(avgtmeanGP == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(avgtmeanGP == 2)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']

final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]



#######################################################################################
########### Avg temp로  나눠서 그래프 그리고 한 그래프에 겹쳐서 나타내기
#######################################################################################


regions_1 <- names(final_data_1)
logRR_temp_1 <- logRRse_temp_1 <- logRR_preci_1 <- logRRse_preci_1 <- vector("numeric",length(regions_1))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_1[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_1[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_1[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_1[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_1, paste0("b",seq(3))))
vcov_temp_1 <- vector("list", length(regions_1))

coef_preci_1 <- matrix(NA, nrow=length(regions_1), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_1, paste0("b",seq(3))))
vcov_preci_1 <- vector("list", length(regions_1))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_1[i,] <- coef(cr_temp)
  vcov_temp_1[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_1[i,] <- coef(cr_preci)
  vcov_preci_1[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_1 <- covariates_1 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_1 <- ndvi_vector_summer_1[,"ndvi_summer"]

UrbanRural_1 <- covariates_1 %>% 
  dplyr::select(UrbanRural)

UrbanRural_1 <- UrbanRural_1[,"UrbanRural"]

species_1 <- covariates_1 %>% 
  dplyr::select(species)

species_1 <- species_1[,"species"]

avgpreci_1 <- covariates_1 %>% 
  dplyr::select(avgpreci)

avgpreci_1 <- avgpreci_1[,"avgpreci"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_1 <- mixmeta(coef_temp_1~ndvi_vector_summer_1 + species_1, vcov_temp_1, method="ml") # speceis는 제외.
print(summary(mix_temp_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_1 <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp_1 <- seq(bound_temp_1[1],bound_temp_1[2],by=0.1)
argvar_temp_1=list(fun="ns", knots=quantile(xvar_temp_1,prob=c(.33,.66)))
bvar_temp_1 <- do.call(onebasis, c(list(x=xvar_temp_1), argvar_temp_1))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_1 <- data.frame(ndvi_vector_summer_1=covariates_1$ndvi_summer,
                        species_1=covariates_1$species
) 

predmix_temp_1 <- predict(mix_temp_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_temp_1[[i]]$fit
  count <- count + 1
}
mix_temp_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_temp_1[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_1 <- (vcov / count)


pred.pool_temp_1 <- crosspred(bvar_temp_1, coef=mix_temp_fit_1, vcov=mix_temp_vcov_1, model.link="log", 
                              by=0.1, cen=10)

pred.reg_temp_1 <- lapply(seq(regions_1), 
                          function(i) crosspred(bvar_temp_1, coef=predmix_temp_1[[i]]$fit, vcov=predmix_temp_1[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_1 <- mixmeta(coef_preci_1~ndvi_vector_summer_1 + species_1, vcov_preci_1, method="ml")
print(summary(mix_preci_1), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_1 <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci_1 <- seq(bound_preci_1[1],bound_preci_1[2],by=0.1)
argvar_preci_1=list(fun="ns", knots=quantile(xvar_preci_1,prob=c(.33,.66)))
bvar_preci_1 <- do.call(onebasis, c(list(x=xvar_preci_1), argvar_preci_1))


predmix_preci_1 <- predict(mix_preci_1, datanew_1,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_1)) {
  fit <- fit + predmix_preci_1[[i]]$fit
  count <- count + 1
}
mix_preci_fit_1 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_1)) {
  vcov <- vcov + predmix_preci_1[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_1 <- (vcov / count)


pred.pool_preci_1 <- crosspred(bvar_preci_1, coef=mix_preci_fit_1, vcov=mix_preci_vcov_1, model.link="log", 
                               by=0.1, cen=3)

pred.reg_preci_1 <- lapply(seq(regions_1), 
                           function(i) crosspred(bvar_preci_1, coef=predmix_preci_1[[i]]$fit, vcov=predmix_preci_1[[i]]$vcov, 
                                                 model.link="log", cen=3))



#############################################################################



regions_2 <- names(final_data_2)
logRR_temp_2 <- logRRse_temp_2 <- logRR_preci_2 <- logRRse_preci_2 <- vector("numeric",length(regions_2))

###################################### Pooled
### varknot 같은 객체때문에 먼저 적합합
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 3
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 66th pct
  target_preci <- as.character(60) 
  logRR_temp_2[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp_2[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci_2[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci_2[i] <- pred.preci$allse[target_preci] # standard error
}

#####################################################################################

coef_temp_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_temp)+1,
                      dimnames=list(regions_2, paste0("b",seq(3))))
vcov_temp_2 <- vector("list", length(regions_2))

coef_preci_2 <- matrix(NA, nrow=length(regions_2), ncol=length(varknot_preci)+1,
                       dimnames=list(regions_2, paste0("b",seq(3))))
vcov_preci_2 <- vector("list", length(regions_2))


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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi + auto, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  # Reduced parameters for the overall cumulative association
  cr_temp <- crossreduce(cb.temp, model) 
  coef_temp_2[i,] <- coef(cr_temp)
  vcov_temp_2[[i]] <- vcov(cr_temp)
  cr_preci <- crossreduce(cb.preci, model) 
  coef_preci_2[i,] <- coef(cr_preci)
  vcov_preci_2[[i]] <- vcov(cr_preci)
}



## ------------------------------------------------------------------------
# Meta-regression by NDVI, avgtemp, speceis
### 여기서, covariates_1이랑 final_data_1이랑 순서 맞추고 가야하는데, 순서 맞음
ndvi_vector_summer_2 <- covariates_2 %>% 
  dplyr::select(ndvi_summer)

ndvi_vector_summer_2 <- ndvi_vector_summer_2[,"ndvi_summer"]

UrbanRural_2 <- covariates_2 %>% 
  dplyr::select(UrbanRural)

UrbanRural_2 <- UrbanRural_2[,"UrbanRural"]

species_2 <- covariates_2 %>% 
  dplyr::select(species)

species_2 <- species_2[,"species"]

avgpreci_2 <- covariates_2 %>% 
  dplyr::select(avgpreci)

avgpreci_2 <- avgpreci_2[,"avgpreci"]

## ------------------------------------------------------------------------
# multivariate meta-regression including NDVI, avgtemp, species distribution
############################# Temp
mix_temp_2 <- mixmeta(coef_temp_2~ndvi_vector_summer_2 +  species_2, vcov_temp_2, method="ml") # speceis는 제외.
print(summary(mix_temp_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_temp를 그대로 사용하게되면  bvar_temp의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_temp_2 <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp_2 <- seq(bound_temp_2[1],bound_temp_2[2],by=0.1)
argvar_temp_2=list(fun="ns", knots=quantile(xvar_temp_2,prob=c(.33,.66)))
bvar_temp_2 <- do.call(onebasis, c(list(x=xvar_temp_2), argvar_temp_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.
datanew_2 <- data.frame(ndvi_vector_summer_2=covariates_2$ndvi_summer,
                        species_2=covariates_2$species) 

predmix_temp_2 <- predict(mix_temp_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_temp_2[[i]]$fit
  count <- count + 1
}
mix_temp_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_temp_2[[i]]$vcov
  count <- count + 1
}
mix_temp_vcov_2 <- (vcov / count)


pred.pool_temp_2 <- crosspred(bvar_temp_2, coef=mix_temp_fit_2, vcov=mix_temp_vcov_2, model.link="log", ci.level = 0.7,
                              by=0.1, cen=10)

pred.reg_temp_2 <- lapply(seq(regions_2), 
                          function(i) crosspred(bvar_temp_2, coef=predmix_temp_2[[i]]$fit, vcov=predmix_temp_2[[i]]$vcov, 
                                                model.link="log", cen=10))




############################# Precipitation
mix_preci_2 <- mixmeta(coef_preci_2~ndvi_vector_summer_2 + species_2, vcov_preci_2, method="ml")
print(summary(mix_preci_2), digits=3)

# Plot the pooled exposure-response curve
# mixmeta로 모든 covariate를 포함한 meta-regression model을 활용해서 대표 그래프를 생성하고 싶음.
# mix_preci를 그대로 사용하게되면  bvar_preci의 basis matrix와 크기가 달라서 불가능.
# 그러므로 먼저 predict로 객체를 생성하고, 해당 객체에 담겨 있는 각 지역 별 값을 다시 평균내서 사용하고자 함.

bound_preci_2 <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci_2 <- seq(bound_preci_2[1],bound_preci_2[2],by=0.1)
argvar_preci_2=list(fun="ns", knots=quantile(xvar_preci_2,prob=c(.33,.66)))
bvar_preci_2 <- do.call(onebasis, c(list(x=xvar_preci_2), argvar_preci_2))

## 그럼 지금 이거로 predict 한거에서 각 지역들의 coef랑 vcov를 꺼내서 평균내면 될듯.

predmix_preci_2 <- predict(mix_preci_2, datanew_2,vcov=T,format="list")

fit <- 0
count <- 0
for (i in seq(regions_2)) {
  fit <- fit + predmix_preci_2[[i]]$fit
  count <- count + 1
}
mix_preci_fit_2 <- (fit / count)

vcov <- 0
count <- 0
for (i in seq(regions_2)) {
  vcov <- vcov + predmix_preci_2[[i]]$vcov
  count <- count + 1
}
mix_preci_vcov_2 <- (vcov / count)


pred.pool_preci_2 <- crosspred(bvar_preci_2, coef=mix_preci_fit_2, vcov=mix_preci_vcov_2, model.link="log", ci.level = 0.9,
                               by=0.1, cen=3)

pred.reg_preci_2 <- lapply(seq(regions_2), 
                           function(i) crosspred(bvar_preci_2, coef=predmix_preci_2[[i]]$fit, vcov=predmix_preci_2[[i]]$vcov, 
                                                 model.link="log", cen=3))



############# Species 별로 나눠서 한 그래프에 넣을 때는 CI area를 이렇게 하는게 더 잘 나타날듯

### temp
plot(pred.pool_temp_2, type="l", ci.arg=list("area", col=adjustcolor("#CC0000", alpha = 0.2)), ylab="RR", ylim=c(0,30), 
     col="#CC0000", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_temp_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#0080FF", alpha = 0.2)),  col="#0080FF", lwd=2)


### Preci
plot(pred.pool_preci_2, type="l", ci.arg=list("area", col=adjustcolor("#CC0000", alpha = 0.2)), ylab="RR", ylim=c(0,8), xlim = c(0, 250),
     col="#CC0000", lwd=2, xlab="Temperature (C)", main="Pooled RR of mixmeta")
lines(pred.pool_preci_1,  type="l", ci = "area", ci.arg=list("area", col=adjustcolor("#0080FF", alpha = 0.2)),  col="#0080FF", lwd=2)

