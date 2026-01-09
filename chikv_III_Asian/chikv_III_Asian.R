##############################
## chikv_III_Asian.R
##############################

## 0. 加载 R 包
library(readr)
library(dplyr)
library(stringr)
library(ggtree)
library(ggplot2)
library(ape)
library(lubridate)  # 处理日期（可选）

##############################
## 1. 读入系统发育树
##############################

tree <- ape::read.tree("/scr/u/dongw21/Chikungunya/chikv_III_Asian/chikv_III_Asian.treefile")
tree
# Phylogenetic tree with 101 tips and 99 internal nodes.

##############################
## 2. 读入元数据并整理
##############################
## metadata_chikv_rr： Accession_ID ID  date  location  genotype  host  gender  age

meta <- read_tsv("/scr/u/dongw21/Chikungunya/metadata_chikv_rr.tsv",
                 show_col_types = FALSE)

meta <- read_tsv("/scr/u/dongw21/Chikungunya/metadata_chikv_rr.tsv",
                 show_col_types = FALSE)

meta <- meta %>%
  dplyr::rename(
    Accession_ID       = Accession_ID,
    Collection_date = Collection_date  )%>%
  mutate(
    Collection_date = ymd(Collection_date))
  
str(meta$Collection_date)
##############################
## 3. 从 tip label 中提取 Accession_ID ID
##############################
tip_df <- tibble(
  label = tree$tip.label
) %>%
  mutate(
    Accession_ID = stringr::str_extract(label, "EPI_ISL_[0-9]+")
  )

##############################
## 4. 把 metadata 合并到 tip 信息里
##############################
## 这次把 metadata 的所有列都带上，方便之后按不同变量着色
tip_df2 <- tip_df %>%
  dplyr::left_join(meta, by = "Accession_ID") %>%
  mutate(
    ## 从 location 中提取国家名（格式：大洲 / 国家 / 省 / 市，取第二部分）
    location_country = ifelse(
      !is.na(location) & location != "unknown",
      stringr::str_trim(stringr::word(location, 2, sep = " / ")),
      location
    )
  )

## 检查匹配
tip_df2 %>%
  select(label, Accession_ID, location, location_country, host, Collection_date) %>%
  head()


##############################
## 5. 按 Location 着色的圆形树
##############################

p_loc <- ggtree(tree, layout = "circular") %<+% tip_df2 +
  geom_tippoint(
    aes(color = location_country),
    size = 1.5,
    shape = 16,
    na.rm = TRUE
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 6. 按 Host 着色的圆形树
##############################

p_host <- ggtree(tree, layout = "circular") %<+% tip_df2 +
  geom_tippoint(
    aes(color = host),
    size = 1.5,
    shape = 16,
    na.rm = TRUE
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

##############################
## 7. 按 Date（Collection_date）着色的圆形树
##############################

## 若 Collection_date 成功转为 Date，则用连续色标：
p_date <- ggtree(tree, layout = "circular") %<+% tip_df2 +
  geom_tippoint(
    aes(color = Collection_date),
    size = 1.5,
    shape = 16,
    na.rm = TRUE
  ) +
  scale_color_viridis_c(
    option  = "plasma",
    na.value = "grey80",
    name    = "Collection date",
    breaks  = scales::pretty_breaks(5),                  # 生成 5 个合适的刻度
    labels  = function(x) format(as.Date(x), "%Y-%m-%d") # 把内部数值转成日期显示
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8)
  )

## 如果你的 date 实际上是年份（例如 2019、2020），可以改用因子离散色：
## tip_df2 <- tip_df2 %>% mutate(year = as.factor(Collection_date))
## p_date <- ggtree(tree, layout = "circular") %<+% tip_df2 +
##   geom_tippoint(aes(color = year), size = 1.5, shape = 16, na.rm = TRUE) +
##   theme(...)

##############################
## 8. 展示与保存
##############################

print(p_loc)
print(p_host)
print(p_date)

ggsave(
  filename = "/scr/u/dongw21/Chikungunya/chikv_III_Asian/chikv_III_Asian_tree_location_circular.pdf",
  plot     = p_loc,
  width    = 10,
  height   = 10,
  units    = "in"
)

ggsave(
  filename = "/scr/u/dongw21/Chikungunya/chikv_III_Asian/chikv_III_Asian_tree_host_circular.pdf",
  plot     = p_host,
  width    = 10,
  height   = 10,
  units    = "in"
)

ggsave(
  filename = "/scr/u/dongw21/Chikungunya/chikv_III_Asian/chikv_III_Asian_tree_date_circular.pdf",
  plot     = p_date,
  width    = 10,
  height   = 10,
  units    = "in"
)
