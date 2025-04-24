
library(tidyverse)
load(file = "mit_plots/xkk.rdata")
load(file = "mit_plots/xkkCI.rdata")
load(file = "mit_plots/xkkHomo.rdata")
ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)
xkk[2*n, ] <- 0
apply(xkk[1:n, ], 1, function(x)(x - ave1)) %>% t()-> xkk_a
apply(xkk[(n+1):(2*n), ], 1, function(x)(x - ave2)) %>% t() -> xkk_b

data_matrix <- xkk_b
colnames(data_matrix) <- ftrans(seq(0.24,0.9,0.01))
df <- as.data.frame(data_matrix) %>%
  mutate(
    sample_id = row_number(),
    year_school = proxyData$sector$year_school,  # 替换为实际数据
    new_id = proxyData$sector$new_id,
    floor = proxyData$sector$floor
  )

df_long <- df %>%
  pivot_longer(
    cols = -c(sample_id, year_school, new_id, floor),
    names_to = "time",
    values_to = "value"
  ) %>%
  mutate(time = as.numeric(gsub("t", "", time)))  # 提取时间数值

df_long <- df_long %>%
  mutate(year_school = factor(year_school))

# ggplot(df_long, aes(x = time, y = value, group = sample_id)) +
#   geom_line(
#     aes(color = floor),
#     alpha = 0.6,
#     linewidth = 0.5
#   ) + ylim(-5, 5) +
#   scale_color_brewer(palette = "Set1", name = "Floor") +  # 使用离散调色板
#   labs(x = "Time", y = "Parameters", title = "") +
#   theme_minimal()

# 将数据转换为时间序列列表
ts_list <- split(df_long$value, df_long$sample_id)
dtw_dist <- proxy::dist(ts_list, method = "DTW")
hc <- hclust(dtw_dist, method = "ward.D2")

# 划分簇（假设分为3类）
clusters <- cutree(hc, k = 5)
# 将聚类结果合并到原始数据
df_clustered <- df_long %>%
  mutate(cluster = factor(clusters[sample_id]))

ggplot(df_clustered, aes(x = time, y = value)) +
  geom_line(
    aes(group = sample_id, color = cluster),
    alpha = 0.5,
    linewidth = 0.5
  ) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(x = "Time", y = "Value", title = "Time Series Clustering via DTW") +
  theme_minimal()
