# bike data management
library(readr)
dt201805 <- read_csv("bikedata/201805-capitalbikeshare-tripdata.csv")
dt201804 <- read_csv("bikedata/201804-capitalbikeshare-tripdata.csv")
# View(dt201804)
dt201804$`Start station number` %>% unique() -> aa
dt201804$`End station number` %>% unique() -> bb

# compute the distance
dt201804 %>%
  group_by(`Start station number`, `End station number`) %>%
  summarise(connections = n(), dij = mean(Duration),.groups = 'drop') %>% 
  filter(connections >= 10) -> dd

# compute the number of 3-mins-reached stations
n_i1 <- dd %>%
  group_by(`Start station number`) %>% summarise(reached_i = sum(dij < 3*60))
n_i2 <- dd %>%
  group_by(`End station number`) %>% summarise(reach_i = sum(dij < 3*60))

merge(n_i1, n_i2, by.x = "Start station number", by.y = "End station number") -> n_i
n_i$n_i = (n_i$reached_i + n_i$reach_i)/2
names(n_i) = c("Station_number", "reached_i", "reach_i", "n_i")
n_i <- n_i %>% select("Station_number", "n_i")

merge(dd, n_i, by.x = "Start station number", by.y = "Station_number") -> aa
merge(aa, n_i, by.x = "End station number", by.y = "Station_number") -> bb

bb$Z1 = log(pmax(bb$dij/60, 1))
bb$Z2 = log(pmax(bb$dij/60, 1))^2

names(bb) <-c("End station number", "Start station number", "connections", 
              "dij", "Z3", "Z4", "Z1", "Z2")
all_stations <- unique(c(bb$`End station number`, bb$`Start station number`)) %>% sort()
n = length(all_stations)
new_stations <- 1:n
names(new_stations) <- all_stations

zij <- array(0, c(n, n, 4))
C_matrix <- array(0, c(n, n))
for(e in 1:nrow(bb)){
  x = new_stations[as.character(bb[e, "Start station number"])]
  y = new_stations[as.character(bb[e, "End station number"])]
  z1 = bb[e, "Z1"]
  z2 = bb[e, "Z2"]
  z3 = bb[e, "Z3"]
  z4 = bb[e, "Z4"]
  zij[x, y, 1] = z1
  zij[x, y, 2] = z2
  zij[x, y, 3] = z3
  zij[x, y, 4] = z4
  C_matrix[x, y] = 1
  cat(x, y, "\n")
}

dt201805 %>% filter(`Start station number` %in% all_stations, 
                    `End station number` %in% all_stations) %>% 
  select(`Start station number`, `End station number`, `Start date`)-> cc

cc$x = new_stations[as.character(cc$`Start station number`)]
cc$y = new_stations[as.character(cc$`End station number`)]
cc$time = cc$`Start date` - cc$`Start date`[1]
cc$time = as.numeric(cc$time)/max(as.numeric(cc$time))

trail <- cc %>% select(x, y, time)

bike_data <- list(trail = trail, zij = zij, C_matrix = C_matrix)
save(bike_data, file = "bikedata/bike_data.rdata")

###
# 只研究四月份超过十次以上connection的那些边

