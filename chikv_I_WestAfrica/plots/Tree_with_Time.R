##############################
## 0. 加载 R 包
##############################
library(readr)
library(dplyr)
library(stringr)
library(ggtree)
library(ggplot2)
library(ape)
library(lubridate)
library(Biostrings)

##############################
## 1. 读入 TreeTime 输出的“修剪后时间树”
##############################
tree <- ape::read.nexus(
  "/scr/u/dongw21/Chikungunya/treetime_out_pruned/timetree.nexus"
)
tree

##############################
## 2. 读入 metadata 并整理
##############################
## 假设列名：Accession_ID, Collection_date, location, host, ...
meta <- read_tsv(
  "/scr/u/dongw21/Chikungunya/metadata_chikv_rr.tsv",
  show_col_types = FALSE
)

meta <- meta %>%
  rename(
    Accession_ID    = Accession_ID,
    Collection_date = Collection_date
  ) %>%
  mutate(
    Collection_date = ymd(Collection_date)
  )

##############################
## 3. 从 tip label 中提取 Accession_ID
##############################
tip_df <- tibble(label = tree$tip.label) %>%
  mutate(
    Accession_ID = str_extract(label, "EPI_ISL_[0-9]+")
  )

##############################
## 4. 合并 metadata
##############################
tip_df2 <- tip_df %>%
  left_join(meta, by = "Accession_ID")

# 检查是否有没匹配上的
cat("未匹配 metadata 的 tip 数：",
    sum(is.na(tip_df2$Collection_date)), "\n")

##############################
## 5. 时间树：按 location 着色
##############################
p_time_loc <- ggtree(tree) %<+% tip_df2 +
  geom_tippoint(
    aes(color = location),
    size = 1.5,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  theme_tree2() +                         # 带时间轴
  scale_x_continuous(name = "Time (year)") +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 6. 时间树：按 host 着色
##############################
p_time_host <- ggtree(tree) %<+% tip_df2 +
  geom_tippoint(
    aes(color = host),
    size = 1.5,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Time (year)") +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 7. 时间树：按采集日期渐变着色
##############################
p_time_date <- ggtree(tree) %<+% tip_df2 +
  geom_tippoint(
    aes(color = Collection_date),
    size = 1.5,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  scale_color_viridis_c(
    option   = "plasma",
    na.value = "grey80",
    name     = "Collection date",
    breaks   = scales::pretty_breaks(5),
    labels   = function(x) format(as.Date(x), "%Y-%m-%d")
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Time (year)") +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 8. 同一棵时间树的圆形布局（可选）
##############################
p_circ_loc <- ggtree(tree, layout = "circular") %<+% tip_df2 +
  geom_tippoint(
    aes(color = location),
    size = 1.3,
    na.rm = TRUE
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 9. 显示和保存
##############################
print(p_time_loc)
print(p_time_host)
print(p_time_date)
print(p_circ_loc)

outdir <- "/scr/u/dongw21/Chikungunya/plots"
dir.create(outdir, showWarnings = FALSE)

ggsave(file.path(outdir, "timetree_location.pdf"),
       p_time_loc, width = 10, height = 6)
ggsave(file.path(outdir, "timetree_host.pdf"),
       p_time_host, width = 10, height = 6)
ggsave(file.path(outdir, "timetree_date.pdf"),
       p_time_date, width = 10, height = 6)
ggsave(file.path(outdir, "timetree_circular_location.pdf"),
       p_circ_loc, width = 10, height = 10)

## 如需 PNG：
ggsave(file.path(outdir, "timetree_location.png"),
       p_time_loc, width = 10, height = 6, dpi = 300)
