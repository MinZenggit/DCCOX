library(readxl)
library(dplyr)
haversine_distance <- function(x) {
  lat1=x[1]; lon1=x[2]; lat2=x[3]; lon2=x[4]
  # 将经纬度从度数转换为弧度
  radian_conv <- function(degree) degree * pi / 180
  lat1 <- radian_conv(lat1)
  lon1 <- radian_conv(lon1)
  lat2 <- radian_conv(lat2)
  lon2 <- radian_conv(lon2)
  # 经纬度差值
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  # Haversine 公式
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  # 地球平均半径（单位：米）
  R <- 6371e3
  # 返回距离（单位：米）
  return(R * c)
}

haversine_distance(c(38.7661887,-77.1688199, 39.123513,-77.168902))

dt1 <- read_excel("output3.xlsx")
dt2 <- read.csv("Capital_Bikeshare_Locations.csv")
dt3 <- read.csv("bikedata/201805-capitalbikeshare-tripdata.csv")

merge(dt1, dt2[, c("NAME","LATITUDE", "LONGITUDE")], 
      by.x = "address", by.y = "NAME", all.x = T) -> dd
which(is.na(dd$LATITUDE)) -> NA_ind
dd$matched = !is.na(dd$LATITUDE)
dd$LATITUDE[NA_ind] = dd$lat[NA_ind]
dd$LONGITUDE[NA_ind] = dd$lon[NA_ind]
dd$Errors = apply(dd[, c("lat", "lon", "LATITUDE", "LONGITUDE")], 1, haversine_distance)
write.csv(dd, file = "locations.csv")

View(dd[, c("address", "name", "lat", "lon", "LATITUDE", "LONGITUDE", "Errors", "matched")])

View(dd[(dd$address %in% dt3$Start.station),
        c("address", "name","lat", "lon",
          "LATITUDE", "LONGITUDE", "Errors", "matched")])

haversine_distance <- function(x) {
  lat1=x[1]; lon1=x[2]; lat2=x[3]; lon2=x[4]
  # 将经纬度从度数转换为弧度
  radian_conv <- function(degree) degree * pi / 180
  lat1 <- radian_conv(lat1)
  lon1 <- radian_conv(lon1)
  lat2 <- radian_conv(lat2)
  lon2 <- radian_conv(lon2)
  # 经纬度差值
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  # Haversine 公式
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  # 地球平均半径（单位：米）
  R <- 6371e3
  # 返回距离（单位：米）
  return(R * c)
}


