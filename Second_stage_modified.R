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
library(readr)

## rda 파일은 그냥 load 써야함. 화살표로 객체에 할당하면 안됨.
load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

names(final_data)
final_data <- final_data[order(names(final_data))]
       

final_data_1 <- final_data[c("Ansoeng", "Daegu", "Hwaseong", "Jeongeup", "Jeonju", "Suncheon")]
final_data_2 <- final_data[c("Busan","Gimhae", "Haman", "Hapcheon", "Jeju", "Ulsan")]

ndvi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\ndvi_total_mat.csv")

ndvi_1 <- ndvi %>% 
  mutate(ndvi = (Jul + Aug + Sep + Oct + Nov) / 5 ) %>% 
  dplyr::select(NAME_2, ndvi) %>% 
  arrange(NAME_2)

ndvi_1$NAME_2 %in% names(final_data) ## final_data의 list 요소들 이름과, ndvi의 순서가 같은 것 확인.

final_data_1_1 <- final_data_1 
final_data_2_1 <- final_data_2 

holi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi.csv")


######################################## final data 조정해서 분석 가능한 형태로 바꾸기
### 2019년까지만 사용.

for (i in 1:length(final_data_1)){
  final_data_1[[i]] <- final_data_1[[i]] %>% filter(year < 2020)
  final_data_1[[i]]$holi <- holi$holi
}

for (i in 1:length(final_data_2)){
  final_data_2[[i]] <- final_data_2[[i]] %>% filter(year < 2020)
  final_data_2[[i]]$holi <- holi$holi
}



#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### 경남 이외 지역 final_data_1
regions_1 <- names(final_data_1)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_1))



par(mfrow=c(1,2))
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 10
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_1[i])
  plot(pred.preci, "overall", xlim = c(0,300), ylim=c(0,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_1[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 90th pct
  target_preci <- as.character(60) 
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
       main="Heat effects on tsutsugamushi (33th vs. 66th)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_1, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (10th vs. 90th)")


## ------------------------------------------------------------------------
# Meta-regression by NDVI

ndvi_vector_1 <- ndvi_1 %>% 
  dplyr::filter(NAME_2 == "Jeonju" |
                  NAME_2 == "Jeongeup" |
                  NAME_2 == "Suncheon" | 
                  NAME_2 == "Daegu" |
                  NAME_2 == "Ansoeng" |
                  NAME_2 == "Hwaseong") 

ndvi_vector_1 <- ndvi_vector_1[,"ndvi"]

###################### Temp 
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_1)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_1, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$pred))
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.lb), lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.ub), lty="dashed")
abline(h=0, lty = "dotted")



###################### Preci 
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_1)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_1, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$pred)
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.lb, lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.ub, lty="dashed")
abline(h=1, lty = "dotted")






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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
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
bound_temp <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

pred.pool <- crosspred(bvar_temp, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_temp, coef=coef_temp[i,], vcov=vcov_temp[[i]], 
                                         model.link="log", cen=10))
par(mfrow=c(1,1))
plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2, col = "blue",
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions_1)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound_preci <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

pred.pool <- crosspred(bvar_preci, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_preci, coef=coef_preci[i,], vcov=vcov_preci[[i]], 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,50), lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions_1)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)




## ------------------------------------------------------------------------
# multivariate meta-regression by NDVI
############################# Temp
mixndvi_temp <- mixmeta(coef_temp~ndvi_vector_1, vcov_temp, method="ml")
print(summary(mixndvi_temp), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_temp <- predict(mixndvi_temp, data.frame(ndvi_vector_1=range(ndvi_vector_1)), vcov=T)
predmin_temp <- crosspred(bvar_temp, coef=predndvi_temp[[1]]$fit, vcov=predndvi_temp[[1]]$vcov,
                     model.link="log", by=0.1, cen=10)
predmax_temp <- crosspred(bvar_temp, coef=predndvi_temp[[2]]$fit, vcov=predndvi_temp[[2]]$vcov,
                     model.link="log", by=0.1, cen=10)

plot(predmax_temp, type="l", ci.arg=list(density=50,col="#006600"), ylab="RR", ylim=c(0,20), 
     col="#006600", lwd=3, xlab="Temperature (C)", main="Effect modification by NDVI")
lines(predmin_temp, ci="area", ci.arg=list(col="#FF8000", density=40), col="#FF8000", lwd=3)

legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")


############################# Preci
# multivariate meta-regression by NDVI
mixndvi_preci <- mixmeta(coef_preci~ndvi_vector_1, vcov_preci, method="ml")
print(summary(mixndvi_preci), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_preci <- predict(mixndvi_preci, data.frame(ndvi_vector_1=range(ndvi_vector_1)), vcov=T)
predmin_preci <- crosspred(bvar_preci, coef=predndvi_preci[[1]]$fit, vcov=predndvi_preci[[1]]$vcov,
                     model.link="log", by=0.1, cen=10)
predmax_preci <- crosspred(bvar_preci, coef=predndvi_preci[[2]]$fit, vcov=predndvi_preci[[2]]$vcov,
                     model.link="log", by=0.1, cen=10)

plot(predmax_preci, type="l", ylab="RR", ylim=c(0,50), ci = "lines",
     col="#006600", lwd=3, xlab="Precipitation (mm)", main="Effect modification by NDVI")
lines(predmin_preci, ci = "lines",  col="#FF8000" , lwd=3)


legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")

















#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### 경남 등등 지역 final_data_2
regions_2 <- names(final_data_2)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_2))



par(mfrow=c(1,2))
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 10
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0,30), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_2[i])
  plot(pred.preci, "overall", xlim = c(0,300), ylim=c(0,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_2[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 90th pct
  target_preci <- as.character(60) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}



# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions_2, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on tsutsugamushi (33th vs. 66th)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_2, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (10th vs. 90th)")


## ------------------------------------------------------------------------
# Meta-regression by NDVI

ndvi_vector_2 <- ndvi_1 %>% 
  dplyr::filter(NAME_2 == "Busan" |
                  NAME_2 == "Gimhae" |
                  NAME_2 == "Haman" | 
                  NAME_2 == "Hapcheon" |
                  NAME_2 == "Jeju" |
                  NAME_2 == "Ulsan") 
final_data_2
ndvi_vector_2 <- ndvi_vector_2[,"ndvi"]

###################### Temp 
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_2)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_2, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$pred))
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.lb), lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.ub), lty="dashed")
abline(h=0, lty = "dotted")



###################### Preci 
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_2)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_2, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$pred)
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.lb, lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.ub, lty="dashed")
abline(h=1, lty = "dotted")






################################################################################
# Second-stage modelling - Multivariate meta-analysis
################################################################################

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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
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
bound_temp <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

pred.pool <- crosspred(bvar_temp, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_2), 
                   function(i) crosspred(bvar_temp, coef=coef_temp[i,], vcov=vcov_temp[[i]], 
                                         model.link="log", cen=10))
par(mfrow=c(1,1))
plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2, col = "blue",
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions_2)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound_preci <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

pred.pool <- crosspred(bvar_preci, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_2), 
                   function(i) crosspred(bvar_preci, coef=coef_preci[i,], vcov=vcov_preci[[i]], 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,10), lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions_2)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)




## ------------------------------------------------------------------------
# multivariate meta-regression by NDVI
############################# Temp
mixndvi_temp <- mixmeta(coef_temp~ndvi_vector_2, vcov_temp, method="ml")
print(summary(mixndvi_temp), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_temp <- predict(mixndvi_temp, data.frame(ndvi_vector_2=range(ndvi_vector_2)), vcov=T)
predmin_temp <- crosspred(bvar_temp, coef=predndvi_temp[[1]]$fit, vcov=predndvi_temp[[1]]$vcov,
                          model.link="log", by=0.1, cen=10)
predmax_temp <- crosspred(bvar_temp, coef=predndvi_temp[[2]]$fit, vcov=predndvi_temp[[2]]$vcov,
                          model.link="log", by=0.1, cen=10)

plot(predmax_temp, type="l", ci.arg=list(density=50,col="#006600"), ylab="RR", ylim=c(0,20), 
     col="#006600", lwd=3, xlab="Temperature (C)", main="Effect modification by NDVI")
lines(predmin_temp, ci="area", ci.arg=list(col="#FF8000", density=40), col="#FF8000", lwd=3)

legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")


############################# Preci
# multivariate meta-regression by NDVI
mixndvi_preci <- mixmeta(coef_preci~ndvi_vector_2, vcov_preci, method="ml")
print(summary(mixndvi_preci), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_preci <- predict(mixndvi_preci, data.frame(ndvi_vector_2=range(ndvi_vector_2)), vcov=T)
predmin_preci <- crosspred(bvar_preci, coef=predndvi_preci[[1]]$fit, vcov=predndvi_preci[[1]]$vcov,
                           model.link="log", by=0.1, cen=10)
predmax_preci <- crosspred(bvar_preci, coef=predndvi_preci[[2]]$fit, vcov=predndvi_preci[[2]]$vcov,
                           model.link="log", by=0.1, cen=10)

plot(predmax_preci, type="l", ylab="RR", xlim = c(0,300), ylim=c(0,30),ci.arg=list(density=50,col="#006600"),
     col="#006600", lwd=3, xlab="Precipitation (mm)", main="Effect modification by NDVI")
lines(predmin_preci,ci="lines",  col="#FF8000" , lwd=3)


legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")










##################################################################################################
#### Over_32




######################################## final data 조정해서 분석 가능한 형태로 바꾸기
### 2019년까지만 사용.

### 그리고 해당 기간 묶어주는 변수 만들기. 
### 설, 추석이 포함된 주차에 대한 변수 만들기.

load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

names(final_data)
final_data <- final_data[order(names(final_data))]


final_data_1 <- final_data[c("Ansoeng", "Daegu", "Hwaseong", "Jeongeup", "Jeonju", "Suncheon")]
final_data_2 <- final_data[c("Busan","Gimhae", "Haman", "Hapcheon", "Jeju", "Ulsan")]


holi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\holi_over32.csv")
a <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\a.csv")

final_data_1_1 <- final_data_1 
final_data_2_1 <- final_data_2 


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








#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### 경남 이외 지역 final_data_1
regions_1 <- names(final_data_1)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_1))



par(mfrow=c(1,2))
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 4) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 10
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0,20), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_1[i])
  plot(pred.preci, "overall", xlim = c(0,300), ylim=c(0,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_1[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 90th pct
  target_preci <- as.character(60) 
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
       main="Heat effects on tsutsugamushi (33th vs. 66th)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_1, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (10th vs. 90th)")


## ------------------------------------------------------------------------
# Meta-regression by NDVI

ndvi_vector_1 <- ndvi_1 %>% 
  dplyr::filter(NAME_2 == "Jeonju" |
                  NAME_2 == "Jeongeup" |
                  NAME_2 == "Suncheon" | 
                  NAME_2 == "Daegu" |
                  NAME_2 == "Ansoeng" |
                  NAME_2 == "Hwaseong") 

ndvi_vector_1 <- ndvi_vector_1[,"ndvi"]

###################### Temp 
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_1)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_1, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$pred))
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.lb), lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.ub), lty="dashed")
abline(h=0, lty = "dotted")



###################### Preci 
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_1)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_1, exp(logRR_preci), pch=19, ylim = c(0, 15),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$pred)
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.lb, lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.ub, lty="dashed")
abline(h=1, lty = "dotted")






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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
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
bound_temp <- rowMeans(sapply(final_data_1, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

pred.pool <- crosspred(bvar_temp, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_temp, coef=coef_temp[i,], vcov=vcov_temp[[i]], 
                                         model.link="log", cen=10))
par(mfrow=c(1,1))
plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2, col = "blue",
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions_1)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound_preci <- rowMeans(sapply(final_data_1, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

pred.pool <- crosspred(bvar_preci, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_1), 
                   function(i) crosspred(bvar_preci, coef=coef_preci[i,], vcov=vcov_preci[[i]], 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,50), lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions_1)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)




## ------------------------------------------------------------------------
# multivariate meta-regression by NDVI
############################# Temp
mixndvi_temp <- mixmeta(coef_temp~ndvi_vector_1, vcov_temp, method="ml")
print(summary(mixndvi_temp), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_temp <- predict(mixndvi_temp, data.frame(ndvi_vector_1=range(ndvi_vector_1)), vcov=T)
predmin_temp <- crosspred(bvar_temp, coef=predndvi_temp[[1]]$fit, vcov=predndvi_temp[[1]]$vcov,
                          model.link="log", by=0.1, cen=10)
predmax_temp <- crosspred(bvar_temp, coef=predndvi_temp[[2]]$fit, vcov=predndvi_temp[[2]]$vcov,
                          model.link="log", by=0.1, cen=10)

plot(predmax_temp, type="l",  ylab="RR", ylim=c(0,10), ci = "lines",
     col="#006600", lwd=3, xlab="Temperature (C)", main="Effect modification by NDVI")
lines(predmin_temp, ci = "lines",  col="#FF8000" , lwd=3, col="#FF8000", lwd=3)

legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")


############################# Preci
# multivariate meta-regression by NDVI
mixndvi_preci <- mixmeta(coef_preci~ndvi_vector_1, vcov_preci, method="ml")
print(summary(mixndvi_preci), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_preci <- predict(mixndvi_preci, data.frame(ndvi_vector_1=range(ndvi_vector_1)), vcov=T)
predmin_preci <- crosspred(bvar_preci, coef=predndvi_preci[[1]]$fit, vcov=predndvi_preci[[1]]$vcov,
                           model.link="log", by=0.1, cen=10)
predmax_preci <- crosspred(bvar_preci, coef=predndvi_preci[[2]]$fit, vcov=predndvi_preci[[2]]$vcov,
                           model.link="log", by=0.1, cen=10)

plot(predmax_preci, type="l", ylab="RR", ylim=c(0,50), ci = "lines",
     col="#006600", lwd=3, xlab="Precipitation (mm)", main="Effect modification by NDVI")
lines(predmin_preci, ci = "lines",  col="#FF8000" , lwd=3)


legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")

















#####################################################
# Second-stage modelling - Univariate meta-analysis
#####################################################

#################### 경남 등등 지역 final_data_2
regions_2 <- names(final_data_2)
logRR_temp <- logRRse_temp <- logRR_preci <- logRRse_preci <- vector("numeric",length(regions_2))



par(mfrow=c(1,2))
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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
               family=quasipoisson)
  
  cen_temp = pct_temp[3] # reference temperature at 33th percentile   
  cen_preci = 10
  pred.temp <- crosspred(cb.temp, model, cen=cen_temp, by=1) # default of cen is set as the median 
  pred.preci <- crosspred(cb.preci, model, cen=cen_preci, by=1)
  
  # Plot exposure-response curve
  plot(pred.temp, "overall", ylim=c(0,30), col=2, lwd=2, xlab="Temperature", ylab="RR",
       main=regions_2[i])
  plot(pred.preci, "overall", xlim = c(0,300), ylim=c(0,20), col=2, lwd=2, xlab="Precipitation", ylab="RR",
       main=regions_2[i]) 
  
  # Get risk estimates for heat
  target_temp <- as.character(round(pct_temp[4])) # 90th pct
  target_preci <- as.character(60) 
  logRR_temp[i]   <- pred.temp$allfit[target_temp] # overall cumulative over lag days
  logRRse_temp[i] <- pred.temp$allse[target_temp] # standard error
  logRR_preci[i]   <- pred.preci$allfit[target_preci] # overall cumulative over lag days
  logRRse_preci[i] <- pred.preci$allse[target_preci] # standard error
}



# Random effects meta-analysis
uni_temp <- rma(y=logRR_temp, sei=logRRse_temp, slab=regions_2, measure="RR")
summary(uni_temp)
ci.exp(uni_temp) # combined RR with 95% CI

# Forest plot
par(mfrow=c(1,1))
forest(uni_temp, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Heat effects on tsutsugamushi (33th vs. 66th)")


uni_preci <- rma(y=logRR_preci, sei=logRRse_preci, slab=regions_2, measure="RR")
summary(uni_preci)
ci.exp(uni_preci) # combined RR with 95% CI

# Forest plot
forest(uni_preci, transf=exp, refline=1, pch=23, bg=4, col=2,
       main="Rainfall effects on tsutsugamushi  (10th vs. 90th)")


## ------------------------------------------------------------------------
# Meta-regression by NDVI

ndvi_vector_2 <- ndvi_1 %>% 
  dplyr::filter(NAME_2 == "Busan" |
                  NAME_2 == "Gimhae" |
                  NAME_2 == "Haman" | 
                  NAME_2 == "Hapcheon" |
                  NAME_2 == "Jeju" |
                  NAME_2 == "Ulsan") 
final_data_2
ndvi_vector_2 <- ndvi_vector_2[,"ndvi"]

###################### Temp 
res <- rma(y=logRR_temp, sei=logRRse_temp, mods=ndvi_vector_2)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_temp)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_2, logRR_temp, pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="log(Relative Risk)")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$pred))
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.lb), lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), log(preds$ci.ub), lty="dashed")
abline(h=0, lty = "dotted")



###################### Preci 
res <- rma(y=logRR_preci, sei=logRRse_preci, mods=ndvi_vector_2)
summary(res)

# Bubble plot
preds <- predict(res,  transf = exp)
wi <- 1/sqrt(logRRse_preci)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(ndvi_vector_2, exp(logRR_preci), pch=19, ylim = c(0, 10),
     cex=size, xlab="NDVI", ylab="Relative Risk")

lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$pred)
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.lb, lty="dashed")
lines(c(6000, 6200, 6500, 7000, 7500, 8000), preds$ci.ub, lty="dashed")
abline(h=1, lty = "dotted")






################################################################################
# Second-stage modelling - Multivariate meta-analysis
################################################################################

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
  
  
  model <- glm(cases ~ cb.temp + cb.preci + ns(week, df = 3) + ns(year, df = 3) + holi, sub, ### ppt에 캡쳐한 문제는 여기서 week variable의 df 줄였더니 해결은 됨.
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
bound_temp <- rowMeans(sapply(final_data_2, function(x) range(x$tmean)))
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))

pred.pool <- crosspred(bvar_temp, coef=coef(mv_temp), vcov=vcov(mv_temp), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_2), 
                   function(i) crosspred(bvar_temp, coef=coef_temp[i,], vcov=vcov_temp[[i]], 
                                         model.link="log", cen=10))
par(mfrow=c(1,1))
plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,20), lwd=2, col = "blue",
     xlab="Temperature (C)", main="Pooled and first-stage")
for(i in seq(regions_2)) lines(pred.reg[[i]], col="grey")

lines(pred.pool, lwd=3)



# Multivariate mixed-effects meta-analysis - Precipitation
mv_preci <- mixmeta(coef_preci~1, vcov_preci, method="ml")
summary(mv_preci)

# Plot the pooled exposure-response curve
bound_preci <- rowMeans(sapply(final_data_2, function(x) range(x$total_preci)))
xvar_preci <- seq(bound_preci[1],bound_preci[2],by=0.1)
argvar_preci=list(fun="ns", knots=quantile(xvar_preci,prob=c(.33,.66)))
bvar_preci <- do.call(onebasis, c(list(x=xvar_preci), argvar_preci))

pred.pool <- crosspred(bvar_preci, coef=coef(mv_preci), vcov=vcov(mv_preci), model.link="log", 
                       by=0.1, cen=10)
pred.reg <- lapply(seq(regions_2), 
                   function(i) crosspred(bvar_preci, coef=coef_preci[i,], vcov=vcov_preci[[i]], 
                                         model.link="log", cen=10))

plot(pred.pool, type="l", ci="n", ylab="RR", ylim=c(0,10), lwd=2, col = "blue",
     xlab="Precipitation (mm)", main="Pooled and first-stage")
for(i in seq(regions_2)) lines(pred.reg[[i]], col="grey")
lines(pred.pool, lwd=3)




## ------------------------------------------------------------------------
# multivariate meta-regression by NDVI
############################# Temp
mixndvi_temp <- mixmeta(coef_temp~ndvi_vector_2, vcov_temp, method="ml")
print(summary(mixndvi_temp), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_temp <- predict(mixndvi_temp, data.frame(ndvi_vector_2=range(ndvi_vector_2)), vcov=T)
predmin_temp <- crosspred(bvar_temp, coef=predndvi_temp[[1]]$fit, vcov=predndvi_temp[[1]]$vcov,
                          model.link="log", by=0.1, cen=10)
predmax_temp <- crosspred(bvar_temp, coef=predndvi_temp[[2]]$fit, vcov=predndvi_temp[[2]]$vcov,
                          model.link="log", by=0.1, cen=10)

plot(predmax_temp, type="l", ci.arg=list(density=50,col="#006600"), ylab="RR", ylim=c(0,20), 
     col="#006600", lwd=3, xlab="Temperature (C)", main="Effect modification by NDVI")
lines(predmin_temp, ci="area", ci.arg=list(col="#FF8000", density=40), col="#FF8000", lwd=3)

legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")


############################# Preci
# multivariate meta-regression by NDVI
mixndvi_preci <- mixmeta(coef_preci~ndvi_vector_2, vcov_preci, method="ml")
print(summary(mixndvi_preci), digits=3)

# Plot the pooled exposure-response curve by NDVI
predndvi_preci <- predict(mixndvi_preci, data.frame(ndvi_vector_2=range(ndvi_vector_2)), vcov=T)
predmin_preci <- crosspred(bvar_preci, coef=predndvi_preci[[1]]$fit, vcov=predndvi_preci[[1]]$vcov,
                           model.link="log", by=0.1, cen=10)
predmax_preci <- crosspred(bvar_preci, coef=predndvi_preci[[2]]$fit, vcov=predndvi_preci[[2]]$vcov,
                           model.link="log", by=0.1, cen=10)

plot(predmax_preci, type="l", ylab="RR", xlim = c(0,300), ylim=c(0,30),ci.arg=list(density=50,col="#006600"),
     col="#006600", lwd=3, xlab="Precipitation (mm)", main="Effect modification by NDVI")
lines(predmin_preci,ci="lines",  col="#FF8000" , lwd=3)


legend("topright", c("High","Low"), lty=1, col=c("#006600","#FF8000"), inset=0.05, title="NDVI")









