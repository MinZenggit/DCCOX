
# 设置工作目录到bikedata/rawdata
setwd("bikedata/rawdata")

# 获取所有.zip文件的路径
zip_files <- list.files(pattern = "*.zip", full.names = TRUE)

# 解压缩所有.zip文件到当前目录
for (zip_file in zip_files) {
  unzip(zip_file, exdir = ".")
}

# 获取所有.csv文件的路径
csv_files <- list.files(pattern = "*.csv", full.names = TRUE)

# 读取所有.csv文件到列表中
data_list <- lapply(csv_files, read.csv)

# 提取所有文件的Start station和End station
all_stations <- unlist(lapply(data_list, function(x) c(x$`Start station`, x$`End station`)))

# 获取唯一的站点名称
unique_stations <- unique(all_stations)

# 查看结果
print(unique_stations)



# 安装并加载必要的包
if (!require(hash)) install.packages("hash")
library(hash)

# 设置工作目录到包含数据的文件夹（根据实际路径调整）
setwd("bikedata/rawdata")

# 获取所有.csv文件的路径
csv_files <- list.files(pattern = "*.csv", full.names = TRUE)

# 初始化一个hash集合存储唯一的站点名称
Stations <- c()

# 遍历每个CSV文件
for (csv_file in csv_files) {
  # 打开文件连接
  con <- file(csv_file, "r")
  read.csv(con) -> aa
  unique(c(aa$Start.station,  aa$End.station, 
           aa$start_station_name, aa$end_station_name)) -> dd
  Stations <- c(Stations, dd)
  print(csv_file);  print(length(dd))
}

# 获取所有唯一的站点名称
unique_stations <- unique(Stations)
write.table(unique_stations, file = "all_stations.txt")
# 输出结果




