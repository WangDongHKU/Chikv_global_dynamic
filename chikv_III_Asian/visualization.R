library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tidytree)

# 文件路径
file1 <- "/scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian_alnMCC1.tree"
tree_phylo <- read.tree(file1)
