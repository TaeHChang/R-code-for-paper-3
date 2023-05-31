############ Line and dot plots
library(readr)
library(tidyverse)
library(readxl)

X2001_2022 <- read_excel("C:/Users/Taehee/OneDrive/바탕 화면/Univ Tokyo papers/dataset/2001-2022_modified_import.xlsx")

######### Jeonju
Jeonju <- X2001_2022 %>% 
  dplyr::select(year, week, number, Jeonju)

ggplot(Jeonju, aes(x = number, y = Jeonju)) + geom_line(size = 0.9, alpha= 0.8) + 
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()

######### Namwon
Namwon <- X2001_2022 %>% 
  dplyr::select(year, week, number, Namwon)

Namwon$Namwon <- as.numeric(Namwon$Namwon) 

ggplot(Namwon, aes(x = number, y = Namwon)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()

######### Jeongeup
Jeongeup <- X2001_2022 %>% 
  dplyr::select(year, week, number, Jeongeup)

Jeongeup$Jeongeup <- as.numeric(Jeongeup$Jeongeup) 

ggplot(Jeongeup, aes(x = number, y = Jeongeup)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()

######### Yesan
Yesan <- X2001_2022 %>% 
  dplyr::select(year, week, number, Yesan)

Yesan$Yesan <- as.numeric(Yesan$Yesan) 

ggplot(Yesan, aes(x = number, y = Yesan)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()

######### Okcheon
Okcheon <- X2001_2022 %>% 
  dplyr::select(year, week, number, Okcheon)

Okcheon$Okcheon <- as.numeric(Okcheon$Okcheon) 

ggplot(Okcheon, aes(x = number, y = Okcheon)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()

######### Hamyang
Hamyang <- X2001_2022 %>% 
  dplyr::select(year, week, number, Hamyang)

Hamyang$Hamyang <- as.numeric(Hamyang$Hamyang) 

ggplot(Hamyang, aes(x = number, y = Hamyang)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()


######### Jeju
Jeju <- X2001_2022 %>% 
  dplyr::select(year, week, number, Jeju)

Jeju$Jeju <- as.numeric(Jeju$Jeju) 

ggplot(Jeju, aes(x = number, y = Jeju)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()


######### Hapcheon
Hapcheon <- X2001_2022 %>% 
  dplyr::select(year, week, number, Hapcheon)

Hapcheon$Hapcheon <- as.numeric(Hapcheon$Hapcheon) 

ggplot(Hapcheon, aes(x = number, y = Hapcheon)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()


######### Yangsan
Yangsan <- X2001_2022 %>% 
  dplyr::select(year, week, number, Yangsan)

Yangsan$Yangsan <- as.numeric(Yangsan$Yangsan) 

ggplot(Yangsan, aes(x = number, y = Yangsan)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()


######### Ulsan
Ulsan <- X2001_2022 %>% 
  dplyr::select(year, week, number, Ulsan)

Ulsan$Ulsan <- as.numeric(Ulsan$Ulsan) 

ggplot(Ulsan, aes(x = number, y =Ulsan)) + geom_line(size = 0.9, alpha= 0.8) +
  geom_smooth(method = 'gam', aes(color="gam")) + geom_smooth(method = 'loess', aes(color="loess")) +
  theme_bw()



