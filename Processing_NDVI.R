############################################
############ Processing NDVI value
############################################

library(sf)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(tidyverse)

library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
detectCores()
registerDoParallel(cores=12)
getDoParWorkers()

setwd("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\KOR_adm\\Si-gun")
Korea <- readOGR(dsn = "Korea_Si-gun.shp")
plot(Korea)
crs(Korea) #4326

ndvi <- raster("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\VI_Monthly_1Km_v6\\NDVI\\MOD13A3_NDVI_2001_001.tif")
ndvi <- projectRaster(ndvi, crs = "+proj=longlat +datum=WGS84 +no_defs")



ndvi_list <- list.files(path="C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\VI_Monthly_1Km_v6\\NDVI", pattern = "tif$", full.names = TRUE)

ndvi_list_raster<-list()
for(i in 1:length(ndvi_list)){
  ndvi_list_raster[[i]]<-raster(ndvi_list[i])
  ndvi_list_raster[[i]] <- projectRaster(ndvi_list_raster[[i]], crs = "+proj=longlat +datum=WGS84 +no_defs")
}
plot(ndvi_list_raster[[20]])



ndvi_mat <- matrix(nrow = 161, ncol = length(ndvi_list))
a<-Sys.time()
for(i in 1:length(ndvi_list)){
ndvi_crop <- crop(ndvi_list_raster[[i]], Korea)
ndvi_mat[,i] <- as.vector(raster::extract(ndvi_crop, Korea, fun = mean, na.rm = TRUE))
}
b<-Sys.time()
b-a

ndvi_mat <- as.data.frame(ndvi_mat)
Korea@data$NAME_2
ndvi_mat_1 <- ndvi_mat %>% cbind(Korea@data$NAME_2)
  
write.csv(ndvi_mat_1, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\ndvi_mat.csv" )
library(readr)
ndvi_mat <- read.csv("C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\ndvi_mat.csv")

#################################### Calculate total monthly average / 이미 monthly로 다운로드 받았음
#### 12번째마다 되풀이 된다. 


ndvi_total <- ndvi_mat %>% 
  mutate(Jan = (V1 + V13 + V25 + V37 + V49 + V61 + V73 +
           V85 + V97 + V109 + V121 + V133) /11,
         Feb = (V2 + V14 + V26 + V38 + V50 + V62 + V74 +
                  V86 + V98 + V110 + V122 + V134) /11,
         Mar = (V3 + V15 + V27 + V39 + V51 + V63 + V75 +
                  V87 + V99 + V111 + V123 + V135) /11,
         Apr = (V4 + V16 + V28 + V40 + V52 + V64 + V76 +
                  V88 + V100 + V112 + V124 + V136) /11,
         May = (V5 + V17 + V29 + V41 + V53 + V65 + V77 +
                  V89 + V101 + V113 + V125 + V137) /11,
         Jun = (V6 + V18 + V30 + V42 + V54 + V66 + V78 +
                  V90 + V102 + V114 + V126 + V138) /11,
         Jul = (V7 + V19 + V31 + V43 + V55 + V67 + V79 +
                  V91 + V103 + V115 + V127 + V139) /11,
         Aug = (V8 + V20 + V32 + V44 + V56 + V68 + V80 +
                  V92 + V104 + V116 + V128 + V140) /11,
         Sep = (V9 + V21 + V33 + V45 + V57 + V69 + V81 +
                  V93 + V105 + V117 + V129 + V141) /11,
         Oct = (V10 + V22 + V34 + V46 + V58 + V70 + V82 +
                  V94 + V106 + V118 + V130 + V142) /11,
         Nov = (V11 + V23 + V35 + V47 + V59 + V71 + V83 +
                  V95 + V107 + V119 + V131 + V143) /11,
         Dec = (V12 + V24 + V36 + V48 + V60 + V72 + V84 +
                  V96 + V108 + V120 + V132 + V144) /11,) %>% 
  rename(NAME_2 = V253) %>% 
  dplyr::select(NAME_2, Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)


write.csv(ndvi_total, "C:\\Users\\Taehee\\OneDrive\\바탕 화면\\My papers\\Univ Tokyo papers_작업중\\Scrub_typhus\\NDVI\\ndvi_total_mat.csv" )


#################################### Visualization
library(tmap)

Korea <- st_read("D:\\Environmental data\\Scrub_typhus\\KOR_adm\\Si-gun\\Korea_Si-gun.shp")

Korea <- Korea %>% 
  left_join(ndvi_total, by = "NAME_2")

library("RColorBrewer")
display.brewer.all()

map1 <- tm_shape(Korea) + tm_polygons(col = "Jan", palette = "Greens", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map1

map2 <- tm_shape(Korea) + tm_polygons(col = "Jul", palette = "Greens", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map2

map3 <- tm_shape(Korea) + tm_polygons(col = "Sep", palette = "Greens", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map3

map4 <- tm_shape(Korea) + tm_polygons(col = "Nov", palette = "Greens", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map4



