# bike data management
rm(list = ls())
library(readr)
library(dplyr)
library(stringr)
library(readxl)
result <- read_excel("~/DCCOX/bikedata/result.xlsx")
result[, -1]/60 -> distance
dt201801 <- read.csv("/data/zm/bikedata/201801-capitalbikeshare-tripdata.csv")
dt201802 <- read.csv("/data/zm/bikedata/201802-capitalbikeshare-tripdata.csv")
dt201803 <- read.csv("/data/zm/bikedata/201803-capitalbikeshare-tripdata.csv")
dt201804 <- read.csv("/data/zm/bikedata/201804-capitalbikeshare-tripdata.csv")
dt201805 <- read.csv("/data/zm/bikedata/201805-capitalbikeshare-tripdata.csv")
dt201806 <- read.csv("/data/zm/bikedata/201806-capitalbikeshare-tripdata.csv")
dt201807 <- read.csv("/data/zm/bikedata/201807-capitalbikeshare-tripdata.csv")
dt201808 <- read.csv("/data/zm/bikedata/201808-capitalbikeshare-tripdata.csv")
dt201809 <- read.csv("/data/zm/bikedata/201809-capitalbikeshare-tripdata.csv")
dt201810 <- read.csv("/data/zm/bikedata/201810-capitalbikeshare-tripdata.csv")
dt201811 <- read.csv("/data/zm/bikedata/201811-capitalbikeshare-tripdata.csv")
dt201812 <- read.csv("/data/zm/bikedata/201812-capitalbikeshare-tripdata.csv")

clean_station_spaces <- function(df) {
  target_cols <- c("Start.station", "End.station")
  valid_cols <- intersect(target_cols, names(df))
  df %>%
    mutate(across(all_of(valid_cols), 
                  ~ {
                    if (!is.character(.)) return(.)  # 跳过非字符类型
                    str_trim(., side = "right") %>%    # 1. 去除右侧所有空格
                      str_replace_all("\\s{2,}", " ")  # 2. 合并中间连续空格
                  }))
}


locations <- read.csv("~/DCCOX/bikedata/locations.csv")
locations$stations_id <- locations$X
dt <- rbind(dt201801, dt201802,dt201803, dt201804,
            dt201805, dt201806,dt201807, dt201808,
            dt201809, dt201810,dt201811, dt201812) %>% clean_station_spaces
unique(c(dt$Start.station, dt$End.station)) %in% locations$address %>% table

merge(dt, locations, by.x = "Start.station", by.y = "address", all.x = T) %>%
  select(c("Start.station", "Duration", "Start.date", "End.date", "Start.station.number", 
           "End.station.number", "End.station", "Bike.number", "Member.type", 
           "addr", "LATITUDE", "LONGITUDE", "matched", "Errors", "stations_id")) -> aa
names(aa) <- c("Start.station", "Duration", "Start.date", "End.date", "Start.station.number", 
               "End.station.number", "End.station", "Bike.number", "Member.type", 
               "start_station", "start_lat", "start_long", "start_matched", "start_errors", "start_stations_id")
merge(aa, locations, by.x = "End.station", by.y = "address", all.x = T)%>%
  select(c("Start.station", "Duration", "Start.date", "End.date", "Start.station.number", 
           "End.station.number", "End.station", "Bike.number", "Member.type", 
           "start_station", "start_lat", "start_long", "start_matched", "start_errors", "start_stations_id",
           "addr", "LATITUDE", "LONGITUDE", "matched", "Errors", "stations_id")) -> bb
names(bb) <- c("Start.station", "Duration", "Start.date", "End.date", "Start.station.number", 
               "End.station.number", "End.station", "Bike.number", "Member.type", 
               "start_station", "start_lat", "start_long", "start_matched", "start_errors", "start_stations_id",
               "end_station", "end_lat", "end_long", "end_matched", "end_errors", "end_stations_id")
bb %>% select(c("start_stations_id", "Start.station", "Start.date", 
                "end_stations_id", "End.station", "End.date", "Duration", 
                "Bike.number", "Member.type", 
                "start_lat", "start_long", "start_matched", "start_errors", 
                 "end_lat", "end_long", "end_matched", "end_errors")) -> ready_dt
c(ready_dt$start_stations_id, ready_dt$end_stations_id) %>% unique() %>% sort() -> ids
distance[ids, ids] -> distance
colnames(distance) <- ids
rownames(distance) <- ids

# 假设 time_matrix 是预定义的非对称时间矩阵
compute_station_density <- function(time_matrix, threshold = 3) {
  n <- nrow(time_matrix)
  non_diag <- matrix(TRUE, n, n)
  diag(non_diag) <- FALSE
  # 计算可达站点数量
  out_counts <- rowSums((time_matrix <= threshold) & non_diag, na.rm = TRUE)  # 从i出发可达
  in_counts <- colSums((time_matrix <= threshold) & non_diag, na.rm = TRUE)    # 到达i的可达
  n_i <- (out_counts + in_counts) / 2
  names(n_i) <- rownames(time_matrix)
  return(n_i)
}
n_i <- compute_station_density(distance, 3)

n = length(ids)
new_stations <- 1:n
names(new_stations) <- ids
zij <- array(0, c(n, n, 4))
C_matrix <- array(0, c(n, n))
# ready_dt$x <- ready_dt$y <- c()

for(e in 1:nrow(ready_dt)){
  x_old <- ready_dt[e, "start_stations_id"]
  y_old <- ready_dt[e, "end_stations_id"]
  dxy <- distance[as.character(x_old), as.character(y_old)]
  x = new_stations[as.character(x_old)]
  y = new_stations[as.character(y_old)]
  # ready_dt[e, 'x'] = x
  # ready_dt[e, 'y'] = y
  zij[x, y, 1] = log(max(dxy, 1))
  zij[x, y, 2] = log(max(dxy, 1))^2
  zij[x, y, 3] = log(max(n_i[as.character(x_old)], 1))
  zij[x, y, 4] = log(max(n_i[as.character(y_old)], 1))
  C_matrix[x, y] = 1
  if(e%%10000==1)cat(x, y, "\n")
}


ready_dt$x = new_stations[as.character(ready_dt$start_stations_id)]
ready_dt$y = new_stations[as.character(ready_dt$end_stations_id)]
ready_dt$Start.date <- as.Date(ready_dt$Start.date)
ready_dt$time = ready_dt$Start.date - min(ready_dt$Start.date)
ready_dt$time = as.numeric(ready_dt$time)/max(as.numeric(ready_dt$time))

trail <- ready_dt %>% select(x, y, time)

bike_data <- list(trail = trail, zij = zij, C_matrix = C_matrix, 
                  new_stations = new_stations, ready_dt = ready_dt)

save(bike_data, file = "~/DCCOX/bikedata/bike_data3.rdata")
