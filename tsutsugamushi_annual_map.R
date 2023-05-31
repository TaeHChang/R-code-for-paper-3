
### Tsutsugamushi map annual 
### By Si-gun

library(readxl)
library(dplyr)
library(tmap)
library(ggplot2)
library(sf)

Korea <- st_read("D:\\Environmental data\\Scrub_typhus\\KOR_adm\\Si-gun\\Korea_Si-gun.shp")
X2001_2022 <- read_excel("C:/Users/Taehee/OneDrive/바탕 화면/Univ Tokyo papers/dataset/2001-2022_total_import.xlsx") 

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
         per2019 = (case2019/population) * 100000,
         per2020 = (case2020/population) * 100000,
         per2021 = (case2021/population) * 100000,
         per2022 = (case2022/population) * 100000) %>% 
  dplyr::select(NAME_2, per2001, per2002, per2003, per2004, per2005,
                per2006, per2007, per2008, per2009, per2010,
                per2011, per2012, per2013, per2014, per2015,
                per2016, per2017, per2018, per2019, per2020,
                per2021, per2022)

Korea <- Korea %>% 
  left_join(X2001_2022, by = "NAME_2")

###data <- Korea %>% dplyr::select(NAME_2, 12)
library("RColorBrewer")
display.brewer.all()

map1 <- tm_shape(Korea) + tm_polygons(col = "per2022", palette = "OrRd", style = "cont", alpha = 0.75) +
  tm_scale_bar(size = 1, position = "right") +
  tmap_options(check.and.fix = TRUE)
map1

