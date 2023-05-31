#### Calculating metorological factors with buffers

################# 나름대로 QGIS에서 processing 한 csv를 가지고 각 지역 별로 기상 변수를 넣어보기.
library(readr)
library(tidyverse)
Gun_ASOS <- read_csv("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/Scrub_typhus/Meteorology_daily/Gun-ASOS.csv", 
                     col_types = cols(field_1 = col_skip()))
Si_ASOS <- read_csv("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/Scrub_typhus/Meteorology_daily/Si-ASOS.csv")
X2001_2022_import <- read_csv("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/Scrub_typhus/Meteorology_daily/X2001_2022_import.csv")


Si_ASOS <- Si_ASOS %>% 
  arrange(NAME_2)

Si_split <- split(Si_ASOS, Si_ASOS$NAME_2) 

for (i in 1:length(Si_split)) {
  Si_split[[i]] = left_join(Si_split[[i]], X2001_2022_import, by = "ASOS_number")
  Si_split[[i]] = Si_split[[i]] %>% group_by(Date) %>% 
    summarise(tmean = mean(Average_temp, na.rm = TRUE),
              preci = mean(Total_preci, na.rm = TRUE))
  cat(i, "")
}

Busan <- Si_split[["Busan"]]


Gun_ASOS <- Gun_ASOS %>% 
  arrange(NAME_2)

Gun_split <- split(Gun_ASOS, Gun_ASOS$NAME_2) 

for (i in 1:length(Gun_split)) {
  Gun_split[[i]] = left_join(Gun_split[[i]], X2001_2022_import, by = "ASOS_number")
  Gun_split[[i]] = Gun_split[[i]] %>% group_by(Date) %>% 
    summarise(tmean = mean(Average_temp, na.rm = TRUE),
              preci = mean(Total_preci, na.rm = TRUE))
  cat(i, "")
}

Danyang <- Gun_split[["Danyang"]]

##############################################
####### 여기까지 해서 각 list 안에 지역 별로 daily meteorology data를 할당하였음. 
####### 이제 해야할 일은, 이 list 들의 meteorology data를 weekly로 바꾸고, outcome을 추가하는 것.
##############################################
library(lubridate)
library(readxl)
X2001_2022_outcome <- read_excel("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/Scrub_typhus/Disease_outcome_weekly/dataset/2001-2022_import.xlsx")

### 단면 가지고 해보기.
Busan <- Busan %>% 
  group_by(year = year(Busan$Date), week = week(Busan$Date)) %>% 
  summarise(tmean = mean(tmean, na.rm = TRUE),
            total_preci = sum(preci, na.rm = TRUE)) 

X2001_2022_outcome %>% dplyr::select(1, 2, "Busan") %>% rename(cases = Busan) %>% 
  left_join(Busan, by = c("year", "week"))

names(Gun_split)[1]


####### 본격적으로 list에 집어넣기
Si_data <- Si_split

for (i in 1:length(Si_split)) {
  Si_data[[i]] = left_join(Si_split[[i]]%>% 
                             mutate(year = year(Si_split[[names(Si_split)[i]]]$Date), 
                                    week = week(Si_split[[names(Si_split)[i]]]$Date)) %>% 
                             group_by(year, week) %>% 
                             summarise(tmean = mean(tmean, na.rm = TRUE),
                                       total_preci = sum(preci, na.rm = TRUE)), 
                           X2001_2022_outcome %>% 
                             dplyr::select(1, 2, names(Si_split)[i]) %>% 
                             rename(cases = names(Si_split)[i]), by = c("year", "week"))
  cat(i, "")
}

Si_data[[1]] # good


Gun_data <- Gun_split

for (i in 1:length(Gun_split)) {
  Gun_data[[i]] = left_join(Gun_split[[i]]%>% 
                              mutate(year = year(Gun_split[[names(Gun_split)[i]]]$Date), 
                                     week = week(Gun_split[[names(Gun_split)[i]]]$Date)) %>% 
                              group_by(year, week) %>% 
                              summarise(tmean = mean(tmean, na.rm = TRUE),
                                        total_preci = sum(preci, na.rm = TRUE)), 
                            X2001_2022_outcome %>% 
                              dplyr::select(1, 2, names(Gun_split)[i]) %>% 
                              rename(cases = names(Gun_split)[i]), by = c("year", "week"))
  cat(i, "")
}

Gun_data[1] # good


final_data <- c(Si_data, Gun_data)

for (i in 1:length(final_data)) {
  final_data[[i]] = na.omit(final_data[[i]])
}
for (i in 1:length(final_data)) {
  final_data[[i]]$cases = as.numeric(final_data[[i]]$cases)
}


Andong <- final_data[[1]] # good

Jeonju <- final_data[["Jeonju"]]
Jeongeup <- final_data[["Jeongeup"]]
Ulsan <- final_data[["Ulsan"]]
Jeju <- final_data[["Jeju"]]

save(final_data, file = "C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")


######################### 각 지역마다 NDVI, avgtemp, 지역변수 할당하기
load(file = "C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

ndvi <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\ndvi_total_mat.csv")

ndvi_1 <- ndvi %>% 
  mutate(ndvi_summer = (Jul + Aug + Sep + Oct + Nov) / 5,
         ndvi_total = (Jan + Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec) / 12, 
         ndvi_winter = (Dec + Jan + Feb) / 3) %>% 
  dplyr::select(NAME_2, ndvi_summer, ndvi_total, ndvi_winter) %>% 
  arrange(NAME_2)

ndvi_1$NAME_2 %in% names(final_data) ## final_data의 list 요소들 이름과, ndvi의 순서가 같은 것 확인.

### 각 지역 별 평균온도. 
avgtmean <- as.data.frame(sapply(final_data,function(x) mean(x$tmean,na.rm=T)))
class(avgtmean)
avgtmean <- avgtmean %>% 
  rename(avgtmean = `sapply(final_data, function(x) mean(x$tmean, na.rm = T))`)

avgtmean <- cbind(rownames(avgtmean), avgtmean)
avgtmean <- avgtmean %>% 
  rename(NAME_2 = `rownames(avgtmean)`)

### 각 지역 별 평균 강수량. 
avgpreci <- as.data.frame(sapply(final_data,function(x) mean(x$total_preci,na.rm=T)))
class(avgpreci)
avgpreci <- avgpreci %>% 
  rename(avgpreci = `sapply(final_data, function(x) mean(x$total_preci, na.rm = T))`)

avgpreci <- cbind(rownames(avgpreci), avgpreci)
avgpreci <- avgpreci %>% 
  rename(NAME_2 = `rownames(avgpreci)`)

### 각 지역의 우점종에 따른 구분.
species <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\region.csv")

### 각 지역을 시,군으로 구분. Urban, Rural 구분 하는거. 시를 1로, 군을 2로
covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
UrbanRural <- covariates %>% 
  dplyr::select(NAME_2, UrbanRural)

### 합쳐서 하나로 
covariates <- ndvi_1 %>% 
  left_join(species, by = "NAME_2") %>% 
  left_join(avgtmean, by = "NAME_2") %>% 
  left_join(avgpreci, by = "NAME_2") %>% 
  left_join(UrbanRural, by = "NAME_2")

  

write.csv(covariates, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")

















