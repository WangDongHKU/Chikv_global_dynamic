#!/usr/bin/env Rscript
#remotes::install_github("YuLab-SMU/ggtree")
# 加载必要的包
# 如果包未安装，请运行: install.packages(c("ape", "RColorBrewer", "dplyr"))
setwd("/scr/u/dongw21/Chikungunya/chikv_II_ECSA")
# 尝试加载ggtree，如果失败则尝试安装或使用基础方法
if (!require(ggtree, quietly = TRUE)) {
  cat("警告: ggtree包未安装，尝试安装...\n")
  cat("注意: 分支着色功能需要ggtree包才能实现\n")
  
  # 尝试安装ggtree
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  tryCatch({
    BiocManager::install("ggtree", ask = FALSE, update = FALSE)
    if (require(ggtree, quietly = TRUE)) {
      cat("✓ ggtree包安装成功\n")
      use_ggtree <- TRUE
    } else {
      cat("✗ ggtree包安装失败，将尝试使用phytools包实现分支着色\n")
      use_ggtree <- FALSE
    }
  }, error = function(e) {
    cat("✗ 无法安装ggtree包:", e$message, "\n")
    cat("将尝试使用phytools包实现分支着色（如果可用）\n")
    use_ggtree <- FALSE
  })
} else {
  cat("✓ ggtree包已加载\n")
  use_ggtree <- TRUE
}

library(ape)

# 尝试加载其他包，如果失败则使用基础方法
if (require(RColorBrewer, quietly = TRUE)) {
  use_rcolorbrewer <- TRUE
} else {
  use_rcolorbrewer <- FALSE
  cat("警告: RColorBrewer包未安装，将使用基础颜色\n")
}

if (require(dplyr, quietly = TRUE)) {
  use_dplyr <- TRUE
} else {
  use_dplyr <- FALSE
  cat("警告: dplyr包未安装，将使用基础R函数\n")
}

# 设置工作目录
setwd("/scr/u/dongw21/Chikungunya/chikv_II_ECSA")

# 读取树文件 (NEXUS格式)
cat("正在读取树文件...\n")
tree <- read.nexus("ECSA_MCC.tree")

# 读取元数据
cat("正在读取元数据...\n")
metadata <- read.delim("ECSA_with_coordinates.tsv", sep = "\t", stringsAsFactors = FALSE)

# 确保Accession_ID列名正确
colnames(metadata)[1] <- "Accession_ID"

# 提取tip标签
tip_labels <- tree$tip.label

# 创建tip数据框并匹配元数据
tip_data <- data.frame(
  label = tip_labels,
  stringsAsFactors = FALSE
)

# 匹配元数据
if (use_dplyr) {
  tip_data <- tip_data %>%
    left_join(metadata, by = c("label" = "Accession_ID")) %>%
    mutate(
      Continent = ifelse(is.na(Continent) | Continent == "", "Unknown", Continent),
      Country = ifelse(is.na(Country) | Country == "", "Unknown", Country)
    )
} else {
  # 使用基础R方法
  tip_data <- merge(tip_data, metadata, by.x = "label", by.y = "Accession_ID", all.x = TRUE)
  tip_data$Continent[is.na(tip_data$Continent) | tip_data$Continent == ""] <- "Unknown"
  tip_data$Country[is.na(tip_data$Country) | tip_data$Country == ""] <- "Unknown"
}

# 创建颜色映射 - 使用更美观的颜色方案
# 为Continent创建颜色
unique_continents <- sort(unique(tip_data$Continent))
if (use_rcolorbrewer && length(unique_continents) <= 8) {
  continent_colors <- RColorBrewer::brewer.pal(max(3, length(unique_continents)), "Set2")
  continent_colors <- continent_colors[1:length(unique_continents)]
} else {
  continent_colors <- rainbow(length(unique_continents))
}
names(continent_colors) <- unique_continents

# 为Country创建颜色
unique_countries <- sort(unique(tip_data$Country))
if (use_rcolorbrewer && length(unique_countries) <= 12) {
  country_colors <- RColorBrewer::brewer.pal(max(3, length(unique_countries)), "Set3")
  country_colors <- country_colors[1:length(unique_countries)]
} else {
  country_colors <- rainbow(length(unique_countries))
}
names(country_colors) <- unique_countries

# 绘制树
cat("正在绘制进化树...\n")

if (use_ggtree) {
  # 使用ggtree绘制（更美观）
  library(ggplot2)
  
  # 检查必要的包
  if (!require(tidytree, quietly = TRUE)) {
    cat("错误: tidytree包未安装，ggtree需要此包\n")
    cat("请安装: install.packages('tidytree')\n")
    use_ggtree <- FALSE
  } else {
    library(tidytree)
  }
  
  if (use_ggtree && !use_dplyr) {
    cat("警告: dplyr包未安装，但ggtree需要dplyr\n")
    cat("尝试加载dplyr...\n")
    if (require(dplyr, quietly = TRUE)) {
      use_dplyr <- TRUE
    } else {
      cat("错误: 无法加载dplyr，将使用基础方法\n")
      use_ggtree <- FALSE
    }
  }
  
  if (use_ggtree) {
    # 转换为tibble格式以便处理
    tree_tbl <- as_tibble(tree)
    
    # 提取tip的Country信息
    cat("正在处理节点状态...\n")
    tip_country <- tree_tbl %>%
      dplyr::filter(!is.na(label)) %>%
      dplyr::left_join(tip_data %>% dplyr::select(label, Country), by = "label")
    
    # 检查并修复分支长度（在推断节点状态之前）
    cat("正在检查分支长度...\n")
    if (any(tree$edge.length < 0, na.rm = TRUE)) {
      cat("发现负分支长度，正在修复（设置为0）...\n")
      tree$edge.length[tree$edge.length < 0] <- 0
    }
    if (any(is.na(tree$edge.length)) || any(is.null(tree$edge.length))) {
      cat("发现NA/NULL分支长度，正在修复（设置为0）...\n")
      tree$edge.length[is.na(tree$edge.length) | is.null(tree$edge.length)] <- 0
    }
    
    # 推断内部节点的Country状态
    cat("正在推断内部节点状态...\n")
    tryCatch({
      country_fac <- factor(tip_country$Country)
      names(country_fac) <- tip_country$label
      # 按树的tip顺序重排，并只保留树上存在的tip
      country_fac <- country_fac[tree$tip.label]
      
      ace_res_country <- ape::ace(country_fac, tree, type = "discrete")
      node_states_country <- apply(ace_res_country$lik.anc, 1, which.max)
      country_node <- data.frame(
        node = (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)),
        Country_node = levels(country_fac)[node_states_country],
        stringsAsFactors = FALSE
      )
      cat("成功使用ape::ace推断内部节点Country状态\n")
    }, error = function(e) {
      cat("ape::ace推断Country失败:", e$message, "\n")
      cat("使用备用方法：内部节点使用NA（将显示为灰色）...\n")
      internal_nodes <- (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree))
      country_node <<- data.frame(
        node = internal_nodes,
        Country_node = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    
    # 合并tip和内部节点的状态
    country_all <- rbind(
      tip_country %>% dplyr::transmute(node, Country_branch = as.character(Country)),
      country_node %>% dplyr::transmute(node, Country_branch = as.character(Country_node))
    )
    
    # 准备完整的数据框
    tip_data_full <- tree_tbl %>%
      dplyr::left_join(country_all, by = "node") %>%
      dplyr::left_join(tip_data, by = "label")
    
    # 创建树图：分支和tip点都用Country着色
    cat("正在创建树图...\n")
    p <- ggtree(tree, layout = "rectangular", size = 0.3) %<+% tip_data_full +
      # 分支用Country着色（像参考文件一样）
      aes(color = Country_branch) +
      geom_tree(size = 0.3) +
      # tip点也用Country着色（填充），不显示tip labels
      geom_tippoint(
        aes(fill = Country), 
        size = 2.5, 
        shape = 21, 
        stroke = 0.8,
        alpha = 0.8
      ) +
      scale_color_manual(
        name = "Country",
        values = country_colors,
        na.value = "gray50",
        guide = guide_legend(
          order = 1,
          title.position = "top",
          override.aes = list(size = 2, linetype = 1)
        )
      ) +
      scale_fill_manual(
        name = "Country",
        values = country_colors,
        na.value = "gray50",
        guide = "none"  # 不显示fill图例，因为和color图例重复
      ) +
      theme_tree2() +
      theme(
        legend.position = c(0.95, 0.5),  # 中间右边
        legend.justification = c(1, 0.5),
        legend.box = "vertical",
        legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(20, 20, 20, 20),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.8, "cm")
      ) +
      labs(title = "Phylogenetic Tree of CHIKV II-ECSA Genotype")
    
    # 保存图片
    cat("正在保存图片...\n")
    ggsave("phylogenetic_tree.png", p, width = 20, height = 14, dpi = 300, bg = "white")
  }
  
} else {
  # 使用基础方法（ape包）绘制，尝试实现分支着色
  cat("警告: ggtree不可用，使用基础方法绘制\n")
  cat("注意: 基础方法可能无法完整显示分支着色，建议安装ggtree包\n")
  
  # 准备颜色向量 - 只使用Country
  tip_country_colors <- country_colors[tip_data$Country]
  tip_country_colors[is.na(tip_country_colors)] <- "gray50"
  
  # 尝试推断内部节点状态以实现分支着色
  cat("正在推断内部节点状态以实现分支着色...\n")
  tryCatch({
    country_fac <- factor(tip_data$Country)
    names(country_fac) <- tip_data$label
    country_fac <- country_fac[tree$tip.label]
    
    # 检查并修复分支长度
    if (any(tree$edge.length < 0, na.rm = TRUE)) {
      tree$edge.length[tree$edge.length < 0] <- 0
    }
    if (any(is.na(tree$edge.length))) {
      tree$edge.length[is.na(tree$edge.length)] <- 0
    }
    
    ace_res <- ape::ace(country_fac, tree, type = "discrete")
    node_states <- apply(ace_res$lik.anc, 1, which.max)
    country_levels <- levels(country_fac)
    
    # 为每条边分配颜色（使用子节点的Country）
    edge_country <- character(nrow(tree$edge))
    for (i in 1:nrow(tree$edge)) {
      child_node <- tree$edge[i, 2]
      if (child_node <= ape::Ntip(tree)) {
        # 是tip节点，使用其Country
        tip_label <- tree$tip.label[child_node]
        edge_country[i] <- as.character(country_fac[tip_label])
      } else {
        # 是内部节点，使用推断的Country
        internal_node_idx <- child_node - ape::Ntip(tree)
        edge_country[i] <- country_levels[node_states[internal_node_idx]]
      }
    }
    
    # 转换为颜色
    edge_colors <- country_colors[edge_country]
    edge_colors[is.na(edge_colors)] <- "gray50"
    
    cat("✓ 成功推断节点状态，分支将按Country着色\n")
    use_edge_colors <- TRUE
  }, error = function(e) {
    cat("✗ 推断节点状态失败:", e$message, "\n")
    cat("分支将显示为默认颜色\n")
    edge_colors <<- rep("black", nrow(tree$edge))
    use_edge_colors <<- FALSE
  })
  
  # 创建PDF输出
  pdf("phylogenetic_tree.pdf", width = 20, height = 14)
  
  # 绘制树（如果成功推断节点状态，使用edge.color参数）
  if (use_edge_colors && exists("edge_colors")) {
    plot(tree, type = "phylogram", show.tip.label = FALSE, edge.width = 0.5, edge.color = edge_colors)
  } else {
    plot(tree, type = "phylogram", show.tip.label = FALSE, edge.width = 0.5)
  }
  
  # 添加tip点 - 使用Country颜色填充
  tiplabels(
    pch = 21,
    bg = tip_country_colors,
    col = tip_country_colors,
    cex = 0.8,
    lwd = 1
  )
  
  # 添加标题
  title("Phylogenetic Tree of CHIKV II-ECSA Genotype", cex.main = 1.5)
  
  # 添加图例 - 只显示Country
  legend(
    "right",
    legend = names(country_colors),
    fill = country_colors,
    title = "Country",
    cex = 0.8,
    bg = "white",
    ncol = 1
  )
  
  dev.off()
  
  # 也尝试保存为PNG - 使用cairo类型（不需要X11）
  tryCatch({
    png("phylogenetic_tree.png", width = 20, height = 14, units = "in", res = 300, type = "cairo")
    if (use_edge_colors && exists("edge_colors")) {
      plot(tree, type = "phylogram", show.tip.label = FALSE, edge.width = 0.5, edge.color = edge_colors)
    } else {
      plot(tree, type = "phylogram", show.tip.label = FALSE, edge.width = 0.5)
    }
    tiplabels(
      pch = 21,
      bg = tip_country_colors,
      col = tip_country_colors,
      cex = 0.8,
      lwd = 1
    )
    title("Phylogenetic Tree of CHIKV II-ECSA Genotype", cex.main = 1.5)
    legend(
      "right",
      legend = names(country_colors),
      fill = country_colors,
      title = "Country",
      cex = 0.8,
      bg = "white",
      ncol = 1
    )
    dev.off()
    cat("PNG文件已成功创建\n")
  }, error = function(e) {
    cat("警告: 无法创建PNG文件 (", e$message, ")，已保存PDF格式\n")
  })
}

cat("\n进化树已成功绘制并保存！\n")
cat("生成的文件：phylogenetic_tree.png\n")
cat("Legend包含：Country\n")
