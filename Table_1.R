##### Table 1
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

for (i in 1:length(final_data)){
  final_data[[i]] <- final_data[[i]] %>% filter(year < 2020)
  final_data[[i]]$holi <- holi$holi
  final_data[[i]]$auto <- log(lag(final_data[[i]]$cases, 1) + 0.01)
  final_data[[i]][is.na(final_data[[i]])] <- log(0.01)
}

covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
covariates <- covariates[,-c(1,8, 9, 10)]

summary(covariates)

covariates <- covariates %>% 
  mutate(ndvi_summerGP = ifelse(ndvi_summer < 7000, 1, 
                                ifelse(ndvi_summer  < 7500, 2, 3)),
         avgtmeanGP = ifelse(avgtmean < 12.75, 1, 2))

covariates$ndvi_summerGP <-factor(covariates$ndvi_summerGP, 
                                  levels = c(1, 2, 3))
covariates$avgtmeanGP <-factor(covariates$avgtmeanGP, 
                               levels = c(1, 2))
covariates$species <-factor(covariates$species, 
                            levels = c(1, 2))

covariates_1 <- covariates %>% 
  dplyr::filter(species == 1)
covariates_2 <- covariates %>% 
  dplyr::filter(species == 2)

name_1 <- covariates_1[,'NAME_2']
name_2 <- covariates_2[,'NAME_2']

final_data_1 <- final_data[name_1]
final_data_2 <- final_data[name_2]

#### Subgroup 1
weekly_avg_cases <- matrix(nrow = length(unlist(final_data_1[[1]][,'cases'])), ncol = length(name_1))
weekly_avg_temp <- matrix(nrow = length(unlist(final_data_1[[1]][,'cases'])), ncol = length(name_1))
weekly_avg_preci <- matrix(nrow = length(unlist(final_data_1[[1]][,'cases'])), ncol = length(name_1))
for (i in 1:length(final_data_1)){
  weekly_avg_cases[,i] <- unlist(final_data_1[[i]][,'cases'])
  weekly_avg_temp[,i] <- unlist(final_data_1[[i]][,'tmean'])
  weekly_avg_preci[,i] <- unlist(final_data_1[[i]][,'total_preci'])
}

weekly_avg_cases_vector <- rowMeans(weekly_avg_cases)
weekly_avg_temp_vector <- rowMeans(weekly_avg_temp)
weekly_avg_preci_vector <- rowMeans(weekly_avg_preci)

summary(weekly_avg_cases_vector)
sd(weekly_avg_cases_vector)
summary(weekly_avg_temp_vector)
sd(weekly_avg_temp_vector)
summary(weekly_avg_preci_vector)
sd(weekly_avg_preci_vector)

summary(covariates_1$ndvi_summer)
sd(covariates_1$ndvi_summer)
summary(covariates_1$avgtmean)
sd(covariates_1$avgtmean)

#### Subgroup 2
weekly_avg_cases <- matrix(nrow = length(unlist(final_data_2[[1]][,'cases'])), ncol = length(name_2))
weekly_avg_temp <- matrix(nrow = length(unlist(final_data_2[[1]][,'cases'])), ncol = length(name_2))
weekly_avg_preci <- matrix(nrow = length(unlist(final_data_2[[1]][,'cases'])), ncol = length(name_2))
for (i in 1:length(final_data_2)){
  weekly_avg_cases[,i] <- unlist(final_data_2[[i]][,'cases'])
  weekly_avg_temp[,i] <- unlist(final_data_2[[i]][,'tmean'])
  weekly_avg_preci[,i] <- unlist(final_data_2[[i]][,'total_preci'])
}

weekly_avg_cases_vector <- rowMeans(weekly_avg_cases)
weekly_avg_temp_vector <- rowMeans(weekly_avg_temp)
weekly_avg_preci_vector <- rowMeans(weekly_avg_preci)

summary(weekly_avg_cases_vector)
sd(weekly_avg_cases_vector)
summary(weekly_avg_temp_vector)
sd(weekly_avg_temp_vector)
summary(weekly_avg_preci_vector)
sd(weekly_avg_preci_vector)

summary(covariates_2$ndvi_summer)
sd(covariates_2$ndvi_summer)
summary(covariates_2$avgtmean)
sd(covariates_2$avgtmean)


#### overall
weekly_avg_cases <- matrix(nrow = length(unlist(final_data[[1]][,'cases'])), ncol = length(final_data))
weekly_avg_temp <- matrix(nrow = length(unlist(final_data[[1]][,'cases'])), ncol = length(final_data))
weekly_avg_preci <- matrix(nrow = length(unlist(final_data[[1]][,'cases'])), ncol = length(final_data))
for (i in 1:length(final_data)){
  weekly_avg_cases[,i] <- unlist(final_data[[i]][,'cases'])
  weekly_avg_temp[,i] <- unlist(final_data[[i]][,'tmean'])
  weekly_avg_preci[,i] <- unlist(final_data[[i]][,'total_preci'])
}

weekly_avg_cases_vector <- rowMeans(weekly_avg_cases)
weekly_avg_temp_vector <- rowMeans(weekly_avg_temp)
weekly_avg_preci_vector <- rowMeans(weekly_avg_preci)

summary(weekly_avg_cases_vector)
sd(weekly_avg_cases_vector)
summary(weekly_avg_temp_vector)
sd(weekly_avg_temp_vector)
summary(weekly_avg_preci_vector)
sd(weekly_avg_preci_vector)

summary(covariates$ndvi_summer)
sd(covariates$ndvi_summer)
summary(covariates$avgtmean)
sd(covariates$avgtmean)








