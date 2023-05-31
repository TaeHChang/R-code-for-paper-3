### Figure 1-1 map
library(readxl)
library(dplyr)
library(tmap)
library(ggplot2)
library(sf)

load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

for (i in 1:length(final_data)){
  final_data[[i]] <- final_data[[i]] %>% filter(year < 2020)
  final_data[[i]]$holi <- holi$holi
  final_data[[i]]$auto <- log(lag(final_data[[i]]$cases, 1) + 0.01)
  final_data[[i]][is.na(final_data[[i]])] <- log(0.01)
}

covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
covariates <- covariates[,-c(1,8, 9, 10)]




### Tsutsugamushi map annual 
Korea <- st_read("D:\\Environmental data\\Scrub_typhus\\KOR_adm\\Si-gun\\Korea_Si-gun.shp")
X2001_2022 <- read_excel("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/2001-2022_total_import.xlsx") 

#### 100000명 당 발생률
X2001_2022 <- X2001_2022 %>% 
  mutate(per2001 = (case2001/population) * 100000,
         per2002 = (case2002/population) * 100000,
         per2003 = (case2003/population) * 100000,
         per2004 = (case2004/population) * 100000,
         per2005 = (case2005/population) * 100000,
         per2006 = (case2006/population) * 100000,
         per2007 = (case2007/population) * 100000,
         per2008 = (case2008/population) * 100000,
         per2009 = (case2009/population) * 100000,
         per2010 = (case2010/population) * 100000,
         per2011 = (case2011/population) * 100000,
         per2012 = (case2012/population) * 100000,
         per2013 = (case2013/population) * 100000,
         per2014 = (case2014/population) * 100000,
         per2015 = (case2015/population) * 100000,
         per2016 = (case2016/population) * 100000,
         per2017 = (case2017/population) * 100000,
         per2018 = (case2018/population) * 100000,
         per2019 = (case2019/population) * 100000) %>% 
  mutate(percases = (per2001+ per2002+ per2003+ per2004+ per2005+
         per2006+ per2007+ per2008+ per2009+ per2010+
         per2011+ per2012+ per2013+ per2014+ per2015+
         per2016+ per2017+ per2018+ per2019) / 20) %>% 
  dplyr::select(NAME_2, percases)

Korea <- Korea %>% 
  left_join(X2001_2022, by = "NAME_2")


library("RColorBrewer")
display.brewer.all()

map1 <- tm_shape(Korea) + tm_polygons(col = "percases", palette = "OrRd", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map1



### Species distribution, NDVI, avg temp
Korea <- st_read("D:\\Environmental data\\Scrub_typhus\\KOR_adm\\Si-gun\\Korea_Si-gun.shp")
covariates <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Calculated_data\\covariates.csv")
covariates <- covariates[,-c(1,8, 9, 10)]

Korea_covariates <- Korea %>% 
  left_join(covariates, by = "NAME_2")

## Species distribution
map2 <- tm_shape(Korea_covariates) + tm_polygons(col = "speceis", palette = "Blues", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map2

## NDVI_ summer
map3 <- tm_shape(Korea_covariates) + tm_polygons(col = "ndvi_summer", palette = "Greens", style = "cont", alpha = 0.80) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map3

## Avg_temp
map4 <- tm_shape(Korea_covariates) + tm_polygons(col = "avgtmean", palette = "Oranges", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map4


## Avg_preci
map5 <- tm_shape(Korea_covariates) + tm_polygons(col = "avgpreci", palette = "Blues", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map5


## UrbanRural
map6 <- tm_shape(Korea_covariates) + tm_polygons(col = "UrbanRural", palette = "-Greys", alpha = 0.7) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map6


