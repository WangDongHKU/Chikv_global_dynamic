# ==========================================
# Chikungunya I-WestAfrica 亚型空间传播可视化
# 参考: map_plot_2016-17.R 和 Fig4.R
# ==========================================
#
# 前提条件:
# 1. MCC树文件（West_African_MCC.tree）应包含位置信息（location1, location2）
#    这通常需要BEAST分析时使用phylogeography模型
# 2. 元数据文件（West_African.tsv）包含Accession_ID、Country、Date等信息
# 3. 地图数据文件（TM_WORLD_BORDERS_SIMPL-0.3.shp）用于绘制世界地图
# 4. 已安装所需的R包：seraphim, rgdal, diagram, viridis, raster, tidyverse等
#
# 输出:
# - chikv_I_WestAfrica_MCC.csv: MCC树提取的空间-时间信息
# - map_chikv_I_WestAfrica_YYYYMMDD.pdf: 空间传播地图
# ==========================================
.libPaths("/home/dongw21/R_libs/4.1.2")
.libPaths("/home/dongw21/R_libs/4.2.1")

# 设置工作目录
setwd("/scr/u/dongw21/Chikungunya/chikv_I_WestAfrica")

# 加载所需的R包
# 首先尝试加载seraphim（如果已安装）
if (!require("seraphim", quietly = TRUE)) {
  cat("警告: seraphim包未安装或无法加载\n")
  cat("提示: seraphim需要从GitHub安装: devtools::install_github('sdellicour/seraphim/unix_OS')\n")
}
library("seraphim")
# 加载其他必需的包
library(plyr)
library(rgdal)
library(diagram)
library(viridis)
library(raster)
library(tidyverse)
library(stringr)
library(dplyr)

# 尝试加载treeio作为备用（用于读取BEAST树文件）
if (!require(treeio, quietly = TRUE)) {
  cat("尝试安装treeio包...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("treeio", update = FALSE)
  library(treeio)
}

# ==========================================
# 1. 读取元数据（带经纬度坐标）
# ==========================================
cat("1. 读取元数据（带经纬度坐标）...\n")

# 直接读取带坐标的文件
coord_file <- "/scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/West_African_with_coordinates.tsv"
if (!file.exists(coord_file)) {
  stop(paste("错误: 未找到带坐标的文件:", coord_file, "\n请先运行 add_coordinates_to_tsv.R 生成该文件"))
}

cat("  读取带坐标文件:", coord_file, "\n")
meta_data <- read.delim(coord_file, sep = "\t", stringsAsFactors = FALSE)

# 检查是否包含经纬度列
if (!"Longitude" %in% colnames(meta_data) || !"Latitude" %in% colnames(meta_data)) {
  stop("错误: 元数据文件缺少Longitude或Latitude列")
}

cat("  ✓ 已读取带坐标的数据，包含", nrow(meta_data), "行\n")
cat("  有效坐标记录数:", sum(!is.na(meta_data$Longitude) & !is.na(meta_data$Latitude)), "\n")

# 获取数据中出现的所有国家（用于信息显示）
unique_countries <- unique(meta_data$Country)
cat("  找到", length(unique_countries), "个唯一国家\n")
cat("  国家列表:", paste(unique_countries, collapse = ", "), "\n")

# ==========================================
# 2. 提取MCC树的空间-时间信息
# ==========================================
cat("\n2. 提取MCC树的空间-时间信息...\n")

# 源文件mccExtractions.r
source("/scr/u/dongw21/africa-covid19-genomics-main/mccExtractions.r")

# 创建修复版本的mccExtractions函数（处理位置信息缺失/字段名不同的情况）
mccExtractions_fixed <- function(mcc_tre, mostRecentSamplingDatum) {
  # 小工具：尽量从annotation里解析经纬度
  parse_lon_lat <- function(annotations) {
    lon <- NA; lat <- NA
    # seraphim 标准字段
    if (!is.null(annotations$location2)) lon <- suppressWarnings(as.numeric(annotations$location2[1]))
    if (!is.null(annotations$location1)) lat <- suppressWarnings(as.numeric(annotations$location1[1]))
    # 部分BEAST/TreeAnnotator会把坐标放在 location 或 location.1/location.2
    if (is.na(lon) || is.na(lat)) {
      # 情况1: location 是长度>=2的数值向量
      if (!is.null(annotations$location) && length(annotations$location) >= 2) {
        lon <- suppressWarnings(as.numeric(annotations$location[1]))
        lat <- suppressWarnings(as.numeric(annotations$location[2]))
      }
      # 情况2: location 是以逗号分隔的字符 "lon,lat"
      if ((is.na(lon) || is.na(lat)) && !is.null(annotations$location) && is.character(annotations$location)) {
        if (grepl(",", annotations$location[1])) {
          parts <- suppressWarnings(as.numeric(strsplit(annotations$location[1], ",")[[1]]))
          if (length(parts) >= 2) {
            lon <- parts[1]; lat <- parts[2]
          }
        }
      }
      # 情况3: location.1 / location.2
      if (is.na(lon) && !is.null(annotations$location.1)) lon <- suppressWarnings(as.numeric(annotations$location.1[1]))
      if (is.na(lat) && !is.null(annotations$location.2)) lat <- suppressWarnings(as.numeric(annotations$location.2[1]))
    }
    c(lon, lat)
  }
  
  mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=11)
  colnames(mcc_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat","startNodeL","endNodeL","startYear","endYear")
  mcc_tab[,c("node1","node2")] = mcc_tre$edge
  mcc_tab[,c("length")] = mcc_tre$edge.length
  
  # 初始化位置信息为NA
  mcc_tab[,c("endLon","endLat")] = NA
  mcc_tab[,c("startLon","startLat")] = NA
  
  for (i in 1:length(mcc_tre$annotations)) {
    annotations = mcc_tre$annotations[[i]]
    lonlat <- parse_lon_lat(annotations)
    mcc_tab[i,"endLon"] = lonlat[1]
    mcc_tab[i,"endLat"] = lonlat[2]
    mcc_tab[i,"endNodeL"] = if(is.null(annotations$height)) NA else as.numeric(annotations$height)
  }
  
  for (i in 1:length(mcc_tre$annotations)) {
    index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
    if (length(index) > 0) {
      mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
      mcc_tab[i,"startNodeL"] = mcc_tab[index,"endNodeL"]
    } else {
      annotations = mcc_tre$root.annotation
      lonlat_root <- parse_lon_lat(annotations)
      mcc_tab[i,"startLon"] = lonlat_root[1]
      mcc_tab[i,"startLat"] = lonlat_root[2]
      mcc_tab[i,"startNodeL"] = if(is.null(annotations$height)) NA else as.numeric(annotations$height)
    }
    mcc_tab[i,"startYear"] = mostRecentSamplingDatum - mcc_tab[i,"startNodeL"]
    mcc_tab[i,"endYear"] = mostRecentSamplingDatum - mcc_tab[i,"endNodeL"]
  }
  mcc_tab = mcc_tab[order(mcc_tab[,"startYear"],decreasing=F),]
  mcc_tab1 = mcc_tab[1,]; mcc_tab2 = mcc_tab[2:dim(mcc_tab)[1],]
  mcc_tab2 = mcc_tab2[order(mcc_tab2[,"endYear"],decreasing=F),]
  mcc_tab = rbind(mcc_tab1, mcc_tab2)
  return(mcc_tab)
}

# 查找MCC树文件
mcc_file <- "West_African_MCC.tree"
if (!file.exists(mcc_file)) {
  stop(paste("错误: 未找到MCC树文件:", mcc_file))
}

cat("  读取MCC树文件:", mcc_file, "\n")

# 定义使用treeio读取MCC树的函数
read_mcc_tree_with_treeio <- function() {
  if (!require(treeio, quietly = TRUE)) {
    stop("无法读取MCC树文件。请安装treeio包: BiocManager::install('treeio')")
  }
  
  mcc_tre_beast <- treeio::read.beast(mcc_file)
  mcc_tre_phylo <- as.phylo(mcc_tre_beast)
  
  # 创建seraphim期望的格式
  mcc_tre <- mcc_tre_phylo
  
  # 提取annotations信息（从treeio的data slot）
  if (!is.null(mcc_tre_beast@data) && nrow(mcc_tre_beast@data) > 0) {
    # 获取节点数据
    node_data <- mcc_tre_beast@data
    
    # 创建annotations列表（对应每条边）
    mcc_tre$annotations <- list()
    for (i in 1:nrow(mcc_tre$edge)) {
      node_num <- mcc_tre$edge[i, 2]
      # 在node_data中查找对应的节点
      node_row <- which(node_data$node == node_num)
      if (length(node_row) > 0 && node_row[1] <= nrow(node_data)) {
        mcc_tre$annotations[[i]] <- list(
          location1 = if("location1" %in% colnames(node_data)) {
            as.numeric(node_data$location1[node_row[1]])
          } else if("location" %in% colnames(node_data)) {
            # 如果location是单个值，可能需要解析
            loc_val <- node_data$location[node_row[1]]
            if (is.character(loc_val) && grepl(",", loc_val)) {
              as.numeric(strsplit(loc_val, ",")[[1]][2])
            } else {
              as.numeric(loc_val)
            }
          } else NA,
          location2 = if("location2" %in% colnames(node_data)) {
            as.numeric(node_data$location2[node_row[1]])
          } else if("location" %in% colnames(node_data)) {
            loc_val <- node_data$location[node_row[1]]
            if (is.character(loc_val) && grepl(",", loc_val)) {
              as.numeric(strsplit(loc_val, ",")[[1]][1])
            } else {
              as.numeric(loc_val)
            }
          } else NA,
          height = if("height" %in% colnames(node_data)) {
            as.numeric(node_data$height[node_row[1]])
          } else if("height_median" %in% colnames(node_data)) {
            as.numeric(node_data$height_median[node_row[1]])
          } else NA
        )
      } else {
        mcc_tre$annotations[[i]] <- list(location1 = NA, location2 = NA, height = NA)
      }
    }
    
    # 设置root.annotation（根节点）
    root_node <- mcc_tre$edge[1, 1]
    root_row <- which(node_data$node == root_node)
    if (length(root_row) > 0 && root_row[1] <= nrow(node_data)) {
      mcc_tre$root.annotation <- list(
        location1 = if("location1" %in% colnames(node_data)) {
          as.numeric(node_data$location1[root_row[1]])
        } else NA,
        location2 = if("location2" %in% colnames(node_data)) {
          as.numeric(node_data$location2[root_row[1]])
        } else NA,
        height = if("height" %in% colnames(node_data)) {
          as.numeric(node_data$height[root_row[1]])
        } else NA
      )
    } else {
      mcc_tre$root.annotation <- list(location1 = NA, location2 = NA, height = NA)
    }
    
    cat("  ✓ 使用treeio::read.beast()成功读取MCC树\n")
    cat("  节点数据列名:", paste(colnames(node_data), collapse = ", "), "\n")
  } else {
    cat("  警告: 无法从树文件中提取节点数据\n")
    mcc_tre$annotations <- list()
    mcc_tre$root.annotation <- list(location1 = NA, location2 = NA, height = NA)
  }
  
  return(mcc_tre)
}


mcc_tre <- readAnnotatedNexus(mcc_file)


# 确定最远采样日期（从元数据中获取）
dates <- as.Date(meta_data$Date, format = "%Y-%m-%d")
mostRecentSamplingDatum <- as.numeric(format(max(dates, na.rm = TRUE), "%Y")) + 
  as.numeric(format(max(dates, na.rm = TRUE), "%j")) / 365.25

cat("  最远采样日期:", mostRecentSamplingDatum, "\n")

# 提取MCC树信息
cat("  提取MCC树信息...\n")

# 先检查MCC树是否包含位置信息
has_location_in_tree <- FALSE
if (length(mcc_tre$annotations) > 0) {
  first_annotation <- mcc_tre$annotations[[1]]
  if (!is.null(first_annotation$location1) && !is.null(first_annotation$location2)) {
    has_location_in_tree <- TRUE
    cat("  MCC树包含位置信息（location1, location2）\n")
  }
}

if (!has_location_in_tree) {
  cat("  警告: MCC树中未找到位置信息（location1, location2）\n")
  cat("  这通常意味着BEAST分析时未使用phylogeography模型\n")
  cat("  将尝试从元数据中补充位置信息...\n")
}

mcc_tab <- mccExtractions_fixed(mcc_tre, mostRecentSamplingDatum)


# 检查mcc_tab中的位置信息
cat("  MCC表维度:", dim(mcc_tab), "\n")
cat("  MCC表列名:", paste(colnames(mcc_tab), collapse = ", "), "\n")

# 检查位置信息是否缺失
has_location_info <- !all(is.na(mcc_tab[,"startLon"])) && !all(is.na(mcc_tab[,"endLon"]))

if (!has_location_info) {
  cat("  警告: MCC树中缺少位置信息，尝试从元数据中补充...\n")
  
  # 获取树的tip标签
  tip_labels <- mcc_tre$tip.label
  cat("  找到", length(tip_labels), "个tip标签\n")
  
  # 从tip标签中提取Accession_ID
  tip_df <- data.frame(
    label = tip_labels,
    Accession_ID = str_extract(tip_labels, "EPI_ISL_[0-9]+"),
    stringsAsFactors = FALSE
  )
  
  # 匹配元数据（直接使用TSV文件中的Longitude和Latitude）
  tip_df <- tip_df %>%
    left_join(meta_data, by = "Accession_ID")
  
  cat("  ✓ 使用TSV文件中的经纬度坐标\n")
  
  # 创建tip标签到坐标的映射（用于后续补充位置信息）
  # 创建一个数据框，包含tip标签和对应的坐标
  tip_coords_df <- tip_df %>%
    select(label, Longitude, Latitude) %>%
    filter(!is.na(Longitude) & !is.na(Latitude))
  
  cat("  已为", nrow(tip_coords_df), "个tip分配了坐标\n")
  if (nrow(tip_coords_df) < length(tip_labels)) {
    cat("  警告: 有", length(tip_labels) - nrow(tip_coords_df), "个tip没有坐标信息\n")
  }
  
  # 显示坐标匹配情况
  cat("  坐标匹配情况:\n")
  coord_summary <- tip_df %>%
    group_by(Country) %>%
    summarise(
      n_tips = n(),
      n_with_coords = sum(!is.na(Longitude) & !is.na(Latitude)),
      .groups = "drop"
    )
  print(coord_summary)
  
  # ==========================================
  # 为MCC树添加位置信息
  # ==========================================
  cat("\n  为MCC树节点添加位置信息...\n")
  
  # 创建tip节点编号到坐标的映射
  # tip节点编号是 1 到 Ntip(mcc_tre)
  n_tips <- length(mcc_tre$tip.label)
  n_internal <- ifelse(is.null(mcc_tre$Nnode), 
                       nrow(mcc_tre$edge) - n_tips + 1, 
                       mcc_tre$Nnode)
  n_total_nodes <- n_tips + n_internal
  
  # 创建节点编号到坐标的映射（先处理tip节点）
  node_coords <- matrix(NA, nrow = n_total_nodes, ncol = 2)
  colnames(node_coords) <- c("Longitude", "Latitude")
  
  # 为tip节点分配坐标
  for (i in 1:n_tips) {
    tip_label <- mcc_tre$tip.label[i]
    tip_row <- which(tip_df$label == tip_label)
    if (length(tip_row) > 0 && !is.na(tip_df$Longitude[tip_row[1]]) && !is.na(tip_df$Latitude[tip_row[1]])) {
      node_coords[i, "Longitude"] <- tip_df$Longitude[tip_row[1]]
      node_coords[i, "Latitude"] <- tip_df$Latitude[tip_row[1]]
    }
  }
  
  cat("  已为", sum(!is.na(node_coords[1:n_tips, "Longitude"])), "个tip节点分配坐标\n")
  
  # 为内部节点推断位置（使用子节点的平均位置）
  # 内部节点编号从 n_tips + 1 开始
  cat("  为内部节点推断位置...\n")
  
  # 从树的底部向上遍历，为每个内部节点计算其子节点的平均位置
  # 使用后序遍历（post-order traversal）
  internal_nodes <- (n_tips + 1):n_total_nodes
  
  # 函数：获取节点的所有后代tip节点
  get_descendant_tips <- function(node, tree) {
    if (node <= n_tips) {
      return(node)  # 如果是tip节点，直接返回
    }
    # 找到该节点的所有子节点
    children <- tree$edge[tree$edge[, 1] == node, 2]
    tips <- c()
    for (child in children) {
      if (child <= n_tips) {
        tips <- c(tips, child)
      } else {
        tips <- c(tips, get_descendant_tips(child, tree))
      }
    }
    return(tips)
  }
  
  # 为每个内部节点计算位置（使用其所有后代tip节点的平均位置）
  for (node in internal_nodes) {
    desc_tips <- get_descendant_tips(node, mcc_tre)
    # 获取这些tip节点的坐标
    desc_coords <- node_coords[desc_tips, , drop = FALSE]
    valid_coords <- desc_coords[!is.na(desc_coords[, "Longitude"]) & !is.na(desc_coords[, "Latitude"]), , drop = FALSE]
    
    if (nrow(valid_coords) > 0) {
      # 使用有效坐标的平均值
      node_coords[node, "Longitude"] <- mean(valid_coords[, "Longitude"])
      node_coords[node, "Latitude"] <- mean(valid_coords[, "Latitude"])
    }
  }
  
  cat("  已为", sum(!is.na(node_coords[internal_nodes, "Longitude"])), "个内部节点推断位置\n")
  
  # 更新mcc_tre$annotations，添加location1和location2
  # annotations对应每条边（edge），每条边连接一个父节点和一个子节点
  # 我们需要为每条边的子节点（end node）添加位置信息
  cat("  更新MCC树的annotations...\n")
  
  if (is.null(mcc_tre$annotations) || length(mcc_tre$annotations) == 0) {
    # 如果annotations不存在，创建新的
    mcc_tre$annotations <- list()
    for (i in 1:nrow(mcc_tre$edge)) {
      mcc_tre$annotations[[i]] <- list()
    }
  }
  
  # 为每条边添加位置信息（使用子节点的位置）
  for (i in 1:nrow(mcc_tre$edge)) {
    end_node <- mcc_tre$edge[i, 2]
    if (!is.na(node_coords[end_node, "Longitude"]) && !is.na(node_coords[end_node, "Latitude"])) {
      mcc_tre$annotations[[i]]$location1 <- node_coords[end_node, "Latitude"]  # location1是纬度
      mcc_tre$annotations[[i]]$location2 <- node_coords[end_node, "Longitude"]  # location2是经度
    }
  }
  
  # 更新root.annotation（根节点的位置）
  root_node <- mcc_tre$edge[1, 1]
  if (!is.na(node_coords[root_node, "Longitude"]) && !is.na(node_coords[root_node, "Latitude"])) {
    if (is.null(mcc_tre$root.annotation)) {
      mcc_tre$root.annotation <- list()
    }
    mcc_tre$root.annotation$location1 <- node_coords[root_node, "Latitude"]
    mcc_tre$root.annotation$location2 <- node_coords[root_node, "Longitude"]
  }
  
  cat("  ✓ MCC树位置信息已更新\n")
  
  # 重新提取MCC树信息
  cat("  使用更新后的MCC树重新提取信息...\n")
  mcc_tab <- mccExtractions_fixed(mcc_tre, mostRecentSamplingDatum)
  
  # 检查更新后的位置信息
  has_location_info_updated <- !all(is.na(mcc_tab[,"startLon"])) && !all(is.na(mcc_tab[,"endLon"]))
  if (has_location_info_updated) {
    cat("  ✓ 成功！已为", sum(!is.na(mcc_tab[,"startLon"])), "条边添加了位置信息\n")
  } else {
    cat("  警告: 更新后仍缺少位置信息\n")
  }
}

# 保存MCC表
write.csv(mcc_tab, "chikv_I_WestAfrica_MCC.csv", row.names = FALSE, quote = FALSE)
cat("  MCC表已保存到: chikv_I_WestAfrica_MCC.csv\n")

# ==========================================
# 3. 估算HPD区域（如果可用）
# ==========================================
cat("\n3. 估算HPD区域...\n")

localTreesDirectory <- "extracted_trees"
nberOfExtractionFiles <- 100  # 可以根据实际情况调整
prob <- 0.95
precision <- 0.025
startDatum <- min(mcc_tab[,"startYear"])

# 检查是否有提取的树文件
if (dir.exists(localTreesDirectory)) {
  extracted_files <- list.files(localTreesDirectory, pattern = "\\.txt$", full.names = FALSE)
  if (length(extracted_files) > 0) {
    cat("  找到", length(extracted_files), "个提取的树文件\n")
    nberOfExtractionFiles <- min(nberOfExtractionFiles, length(extracted_files))
    
    tryCatch({
      polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, 
                                                   prob, startDatum, precision))
      cat("  HPD区域计算完成，共", length(polygons), "个多边形\n")
    }, error = function(e) {
      cat("  错误: spreadGraphic2执行失败:", e$message, "\n")
      cat("  将跳过polygons步骤，继续生成地图（不包含HPD区域）\n")
      polygons <<- list()
    })
  } else {
    cat("  警告: extracted_trees目录中未找到提取的树文件\n")
    cat("  将跳过HPD区域计算\n")
    polygons <- list()
  }
} else {
  cat("  警告: extracted_trees目录不存在\n")
  cat("  将跳过HPD区域计算\n")
  polygons <- list()
}

# ==========================================
# 4. 加载空间边界数据
# ==========================================
cat("\n4. 加载空间边界数据...\n")

shapefile_path <- "/scr/u/dongw21/africa-covid19-genomics-main/data/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp"

if (!file.exists(shapefile_path)) {
  stop(paste("错误: 未找到shapefile文件:", shapefile_path))
}

cat("  读取shapefile:", shapefile_path, "\n")
my_spdf <- readOGR(
  dsn = dirname(shapefile_path),
  layer = "TM_WORLD_BORDERS_SIMPL-0.3",
  verbose = FALSE
)

# 移除南极洲（通过名称和纬度过滤）
cat("  移除南极洲...\n")
# 方法1: 通过名称过滤
antarctica_names <- c("Antarctica", "ANTARCTICA", "antarctica")
my_spdf <- my_spdf[!toupper(my_spdf@data$NAME) %in% toupper(antarctica_names), ]

# 方法2: 通过纬度过滤（移除南纬60度以南的区域）
# 获取每个多边形的边界框
bbox_list <- lapply(1:length(my_spdf), function(i) {
  bbox(my_spdf[i, ])
})
# 只保留北边界在南纬60度以北的多边形
keep_indices <- sapply(bbox_list, function(bbox_i) {
  bbox_i[2, 2] > -60  # 北边界纬度 > -60
})
my_spdf <- my_spdf[keep_indices, ]

cat("  已移除南极洲，剩余", length(my_spdf), "个多边形\n")

borders <- my_spdf
template_raster <- my_spdf

cat("  空间边界数据加载完成\n")

# ==========================================
# 5. 定义时间窗口与颜色标尺
# ==========================================
cat("\n5. 定义颜色标尺...\n")

# 数据实际时间范围
data_min_year <- floor(min(mcc_tab[,"startYear"]))
data_max_year <- ceiling(max(mcc_tab[,"endYear"]))
cat("  数据实际时间范围:", data_min_year, "-", data_max_year, "\n")

# 固定展示的时间范围（根据West Africa数据：2015-2024）
minYear <- 2015.0
maxYear <- 2025.0
cat("  可视化时间范围固定为:", minYear, "-", maxYear, "\n")

# 按固定时间窗口过滤数据
# 保留所有与时间范围有重叠的边（只要startYear或endYear在范围内，或边跨越时间范围）
mcc_tab_this <- mcc_tab %>% 
  data.frame() %>%
  dplyr::filter(
    # 边的开始时间在范围内，或结束时间在范围内，或边跨越整个时间范围
    (startYear >= minYear & startYear <= maxYear) |  # 开始时间在范围内
    (endYear >= minYear & endYear <= maxYear) |      # 结束时间在范围内
    (startYear < minYear & endYear > maxYear)         # 边跨越整个时间范围
  )

if (nrow(mcc_tab_this) == 0) {
  stop(paste("错误: 在", minYear, "到", maxYear, "时间范围内没有数据。"))
}

cat("  过滤后的数据行数:", nrow(mcc_tab_this), "\n")
cat("  原始数据行数:", nrow(mcc_tab), "\n")
cat("  被过滤掉的行数:", nrow(mcc_tab) - nrow(mcc_tab_this), "\n")

# 计算颜色索引
# 将endYear限制在[minYear, maxYear]范围内，以确保颜色索引在有效范围内
endYears_clamped <- pmax(minYear, pmin(maxYear, mcc_tab_this[,"endYear"]))
endYears_indices <- (((endYears_clamped - minYear) / (maxYear - minYear)) * 100) + 1
n_number_colours_needed <- max(round(endYears_indices))
if (n_number_colours_needed < 1) {
  n_number_colours_needed <- 1
}

cat("  需要的颜色数:", n_number_colours_needed, "\n")

# 使用magma颜色标尺
colour_scale <- magma(n_number_colours_needed)
endYears_colours <- colour_scale[pmin(length(colour_scale), pmax(1, round(endYears_indices)))]

# 为polygons设置颜色
if (length(polygons) > 0) {
  polygons_colours <- rep(NA, length(polygons))
  for (i in 1:length(polygons)) {
    date <- as.numeric(names(polygons[[i]]))
    if (date > maxYear) next
    polygon_index <- round((((date - minYear) / (maxYear - minYear)) * 100) + 1)
    if (polygon_index > 0 && polygon_index <= length(colour_scale)) {
      polygons_colours[i] <- paste0(substr(colour_scale[polygon_index], 1, 
                                          nchar(colour_scale[polygon_index]) - 2), "15")
    }
  }
} else {
  polygons_colours <- character(0)
}

# ==========================================
# 6. 生成传播历史地图
# ==========================================
cat("\n6. 生成传播历史地图...\n")

output_file <- sprintf('map_chikv_I_WestAfrica_%s.pdf', format(Sys.time(), "%Y%m%d"))
pdf(output_file, width = 16, height = 15, bg = "white")

ptsize <- 0.6
pitjit <- 0.08
par(mar = c(0, 0, 0, 0), oma = c(1.2, 3.5, 1, 0), mgp = c(0, 0.4, 0), 
    lwd = 0.2, bty = "o")

# 绘制地图背景（国家为白色，海洋为浅灰，边界为深灰）
# 限制纬度范围，不显示南极洲（南纬60度以南）
plot(
  template_raster,
  col   = "white",    # 国家填充白色
  box   = FALSE,
  axes  = FALSE,
  colNA = "white",   # 海洋背景白色
  legend = FALSE,
  ylim  = c(-60, 85)  # 限制纬度范围，不显示南极洲
)
plot(
  borders,
  add   = TRUE,
  lwd   = 0.3,
  border = "gray20",
  ylim  = c(-60, 85)  # 限制纬度范围
)

# 绘制HPD区域（如果存在）
if (length(polygons) > 0) {
  for (i in length(polygons):1) {
    if (as.numeric(names(polygons[[i]])) > maxYear) next
    if (!is.na(polygons_colours[i])) {
      plot(polygons[[i]], axes = F, col = polygons_colours[i], add = T, border = NA)
    }
  }
}

# 绘制传播路径
# 统计信息
total_edges <- dim(mcc_tab_this)[1]
edges_with_valid_coords <- 0
edges_skipped_no_coords <- 0

for (i in 1:total_edges) {
  # 检查经纬度是否有效（这是唯一需要跳过的条件）
  if (is.na(mcc_tab_this[i, "startLon"]) || is.na(mcc_tab_this[i, "startLat"]) ||
      is.na(mcc_tab_this[i, "endLon"]) || is.na(mcc_tab_this[i, "endLat"])) {
    edges_skipped_no_coords <- edges_skipped_no_coords + 1
    next
  }
  
  edges_with_valid_coords <- edges_with_valid_coords + 1
  
  # 绘制灰色背景箭头
  curvedarrow(
    cbind(mcc_tab_this[i, "startLon"], mcc_tab_this[i, "startLat"]), 
    cbind(mcc_tab_this[i, "endLon"], mcc_tab_this[i, "endLat"]), 
    arr.length = 0, arr.width = 0, lwd = 2.5, lty = 1, lcol = "grey22", 
    arr.col = NA, arr.pos = FALSE, curve = 0.3, dr = NA, endhead = F
  )
  # 绘制彩色箭头
  curvedarrow(
    cbind(mcc_tab_this[i, "startLon"], mcc_tab_this[i, "startLat"]), 
    cbind(mcc_tab_this[i, "endLon"], mcc_tab_this[i, "endLat"]), 
    arr.length = 0, arr.width = 0, lwd = 2, lty = 1, lcol = endYears_colours[i], 
    arr.col = NA, arr.pos = FALSE, curve = 0.3, dr = NA, endhead = F
  )
}

# 绘制节点点
for (i in dim(mcc_tab_this)[1]:1) {
  # 检查经纬度是否有效（这是唯一需要跳过的条件）
  if (is.na(mcc_tab_this[i, "startLon"]) || is.na(mcc_tab_this[i, "startLat"]) ||
      is.na(mcc_tab_this[i, "endLon"]) || is.na(mcc_tab_this[i, "endLat"])) {
    next
  }
  
  xs <- mcc_tab_this[i, "startLon"]
  ys <- mcc_tab_this[i, "startLat"]
  xe <- jitter(mcc_tab_this[i, "endLon"], pitjit)
  ye <- jitter(mcc_tab_this[i, "endLat"], pitjit)
  
  if (i == 1) {
    points(xs, ys, pch = 16, col = colour_scale[1], cex = ptsize)
    points(xs, ys, pch = 1, col = "gray90", cex = ptsize)
  }
  points(xe, ye, pch = 16, col = endYears_colours[i], cex = ptsize)
  points(xe, ye, pch = 1, col = "gray90", cex = ptsize)
}

# 输出绘制统计信息
cat("  绘制统计:\n")
cat("    总边数:", total_edges, "\n")
cat("    成功绘制的边数:", edges_with_valid_coords, "\n")
cat("    因缺少坐标而跳过的边数:", edges_skipped_no_coords, "\n")

# 添加坐标轴
xrange <- c(xmin(template_raster), xmax(template_raster))
yrange <- c(ymin(template_raster), ymax(template_raster))
rect(xrange[1], yrange[1], xrange[2], yrange[2], xpd = T, lwd = 0.2)
axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), 
     pos = ymin(template_raster), mgp = c(0, 0.2, 0), cex.axis = 0.5, 
     lwd = 0, lwd.tick = 0.2, padj = -0.8, tck = -0.01, col.axis = "gray30")
axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), 
     pos = xmin(template_raster), mgp = c(0, 0.5, 0), cex.axis = 0.5, 
     lwd = 0, lwd.tick = 0.2, padj = 1, tck = -0.01, col.axis = "gray30")

# 添加时间图例（年份格式，专注于颜色显示）
rast <- raster(matrix(nrow = 1, ncol = 2))
rast[1] <- minYear
rast[2] <- maxYear

# 生成年份标签（每2年一个标签）
legend_years <- seq(floor(minYear), ceiling(maxYear), by = 2)

plot(
  rast, legend.only = TRUE, add = TRUE, col = colour_scale, 
  legend.width = 0.5, legend.shrink = 0.35, 
  smallplot = c(0.40, 0.80, 0.14, 0.17),
  legend.args = list(
    text = "", 
    cex = 0.7, 
    line = 0.3, 
    col = "gray30", 
    font = 1
  ),
  horizontal = TRUE,
  axis.args = list(
    at = legend_years,
    labels = legend_years,
    cex.axis = 0.6, 
    lwd = 0, 
    lwd.tick = 0.2, 
    tck = -0.5, 
    col.axis = "gray30", 
    line = 0, 
    mgp = c(0, -0.02, 0)
  )
)

dev.off()

cat("  地图已保存到:", output_file, "\n")
cat("\n完成！\n")

