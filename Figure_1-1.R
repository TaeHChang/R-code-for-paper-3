######## Figure 1-1 line plots
library(readr)
library(tidyverse)
library(readxl)

line_plot <- read_excel("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/2001-2022_modified_import.xlsx") 

### Cases
line_plot <- line_plot %>% 
  dplyr::filter(year < 2020)
for (i in 4:165) { 
  line_plot[,i] <- as.numeric(unlist(line_plot[,i]))
}

line_plot$average <- rowMeans(line_plot[,c(4:165)])
line_plot$sum <- rowSums(line_plot[,c(4:165)])

ggplot(line_plot, aes(x = number, y = sum)) + geom_line(size = 0.9, alpha= 0.8) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


################### Temp, Preci
load("C:/Users/Taehee/OneDrive/바탕 화면/My papers/Univ Tokyo papers_작업중/dataset/final_data_전국.rda")

for (i in 1:length(final_data)){
  final_data[[i]] <- final_data[[i]] %>% filter(year < 2020)
  final_data[[i]]$number <- 1:length(final_data[[i]]$year)
}

### Temp는 point plot으로로
plot(x=final_data[[1]]$number, y = final_data[[1]]$tmean, type = "p", ylim = c(-10,35), col = "#FF9933") 
for (i in 2:40) {
  lines(x=final_data[[i]]$number, y = final_data[[i]]$tmean, type = "p", ylim = c(-10,35), col = "#FF9933")
}

### Preci는 histogram으로
weekly_avg_preci <- matrix(nrow = length(unlist(final_data[[1]][,'cases'])), ncol = length(final_data))
for (i in 1:length(final_data)){
  weekly_avg_preci[,i] <- unlist(final_data[[i]][,'total_preci'])
}
weekly_avg_preci_vector <- rowMeans(weekly_avg_preci)
number <- as.vector(1:length(weekly_avg_preci_vector))
weekly_preci <- as.data.frame(cbind(number, weekly_avg_preci_vector))

ggplot(data = weekly_preci, aes(x = number, y = weekly_avg_preci_vector)) + 
  geom_bar(stat = "identity", fill = "#3399FF") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


