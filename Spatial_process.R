############ Spatial processing

library(tidyverse)
library(readr)

setwd("D:\\Environmental data\\Scrub_typhus\\Meteorology_daily")

X2001_2010 <- read_csv("2001_2010.csv", col_types = cols(Name = col_skip()))
X2011_2020 <- read_csv("2011_2020.csv", col_types = cols(Name = col_skip()))
X2021_2022 <- read_csv("2021_2022.csv", col_types = cols(Name = col_skip()))
meta <- read_csv("META.csv", col_types = cols(Name = col_skip()))

X2001_2022 <- rbind(X2001_2010, X2011_2020, X2021_2022 )
X2001_2022 <- X2001_2022 %>% 
  arrange(Number)

meta_1 <- meta[!(duplicated(meta$Number)),]

total <- left_join(X2001_2022, meta_1, by = "Number")

total[is.na(total)] <- 0

write.csv(total, "X2001_2022.csv")

latlong <- total %>% 
  dplyr::select(Number, lat, long)

latlong <- latlong[!(duplicated(latlong$Number)),]
write.csv(latlong, "latlong.csv")

