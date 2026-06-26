# --- 加载必要包 ---
library(ggvenn)
library(gridExtra)
library(ggplot2)
library(eulerr)
# --- 设置工作目录 ---
# setwd("D:\\文件1-25和文件38，网络毒理+分子模拟\\09.OMIM和genecard数据库并集")

# --- 创建输出文件夹 ---
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# --- 获取所有txt文件（不包含compound.txt）---
files <- list.files(pattern="\\.txt$")
files <- files[files != "compound.txt"]

# --- 读取所有txt文件中的基因信息，整理到geneList ---
geneList <- list()
for (inputFile in files) {
  rt <- read.table(inputFile, header=FALSE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
  geneNames <- unlist(strsplit(as.vector(rt[,1]), " "))
  geneNames <- trimws(geneNames)         # 去掉前后空格
  geneNames <- geneNames[geneNames != ""]# 去掉空字符串
  uniqGene <- unique(geneNames)
  setName <- tools::file_path_sans_ext(basename(inputFile))  # 文件名去扩展名作为集合名
  geneList[[setName]] <- uniqGene
}

# --- 计算所有txt基因名的并集 ---
unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)

# --- 绘制Venn图并保存 ---
pdf(file = file.path(outputDir, "venn.pdf"), width = 6, height = 6)
venn_plot <- ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)], # 适应多个集合
  set_name_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)],
  set_name_size = 6,
  text_size = 4.5
)
# 添加并集个数
grid.arrange(venn_plot, textGrob(paste("Union Gene Count:", unionCount), gp=gpar(fontsize=14)), ncol=1, heights=c(5, 1))
dev.off()

# --- 输出并集基因到output_folder ---
write.table(unionGenes, file = file.path(outputDir, "disease gene.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

setSizes <- sapply(geneList, length)
pdf(file = file.path(outputDir, "set_barplot.pdf"), width=6, height=4)
ggplot(data.frame(Set=names(setSizes), Size=setSizes), aes(x=Set, y=Size, fill=Set)) +
  geom_bar(stat="identity") +
  theme_minimal() + ylab("Gene Count") + xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

# geneList 已经是一个list，每个集合名下是基因名向量
fit <- euler(geneList)

outputDir <- "output_folder"  # 确认你的输出目录
pdf(file = file.path(outputDir, "euler.pdf"), width = 6, height = 6)
plot(fit, fills = list(fill = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3")[seq_along(geneList)], alpha = 0.7),
     labels = list(font = 2), edges = list(lwd = 2), quantities = TRUE)
dev.off()

library(ggplot2)

setSize <- sapply(geneList, length)
commonGenes <- Reduce(intersect, geneList)
nCommon <- length(commonGenes)
df <- data.frame(
  group = names(setSize),
  value = setSize - nCommon
)
df <- rbind(df, data.frame(group = "Common", value = nCommon))

# 保存为PDF
pdf("output_folder/flower_plot.pdf", width=6, height=6) # 指定保存路径和尺寸

print(
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_bar(stat = "identity", width=1) +
    coord_polar(start = 0) +
    theme_minimal() +
    labs(title = "Petal Diagram (Flower Plot)", y = "Gene Count", x = "")
)

dev.off()


cat("处理完成！已输出Venn图和并集基因列表。\n")
