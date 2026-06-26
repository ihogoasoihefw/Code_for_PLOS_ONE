# 如果尚未安装，请取消注释以下安装包的代码
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("dplyr")
# install.packages("broom")
# install.packages("ggsci")
# install.packages("corrplot")
# install.packages("ggridges")

# ------------------ 加载所需包 ------------------
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)
library(ggsci)
library(corrplot)
library(ggridges)

# ------------------ 设置工作目录与数据路径 ------------------
# setwd("D:\\文件1-25和文件38，网络毒理+分子模拟\\24.免疫浸润评分可视化2.差异")
inputFile <- "CIBERSORT-Results.csv"

# ------------------ 数据读取与预处理 ------------------
# 读取数据（假设每行代表一个样本，列为各免疫细胞的比例）
rt <- read.table(inputFile, header = TRUE, sep = ",", 
                 check.names = FALSE, row.names = 1)

# 根据行名判断组别：后缀为 _con 的为 Control，后缀为 _tre 的为 Treat
rt$Group <- ifelse(grepl("_con$", rownames(rt), ignore.case = TRUE), "Control",
                   ifelse(grepl("_tre$", rownames(rt), ignore.case = TRUE), "Treat", NA))
# 过滤掉组别为NA的样本（避免后续错误）
rt <- rt %>% filter(!is.na(Group))
# 添加样本名称列
rt$Sample <- rownames(rt)

# 获取 Control 与 Treat 的样本名，并构造样本顺序：Control -> gap -> Treat
control_samples <- rownames(rt)[rt$Group == "Control"]
treat_samples   <- rownames(rt)[rt$Group == "Treat"]
all_samples_ordered <- c(control_samples, "gap", treat_samples)

# ------------------ 数据转换与过滤 ------------------
data_long <- melt(rt, id.vars = c("Sample", "Group"), 
                  variable.name = "Immune", 
                  value.name = "Fraction")
# 过滤掉在任一组别中无数据的免疫细胞类型（避免统计检验报错）
data_long <- data_long %>%
  group_by(Immune, Group) %>%
  filter(n() > 0) %>%  # 确保每组至少有一个样本
  ungroup() %>%
  group_by(Immune) %>%
  filter(n_distinct(Group) == 2) %>%  # 确保每个免疫细胞都有两组数据
  ungroup()
# 指定因子水平，确保 gap 保留
data_long$Sample <- factor(data_long$Sample, levels = all_samples_ordered)

# -------- 导出表格1：长格式数据 --------
write.csv(data_long, file = "barplot_data_long.csv", row.names = FALSE)
cat("Barplot 长格式数据已保存为：barplot_data_long.csv\n")

# ------------------ 绘制箱线图 ------------------
# 计算每组样本数
countControl <- length(unique(data_long$Sample[data_long$Group == "Control"]))
countTreat   <- length(unique(data_long$Sample[data_long$Group == "Treat"]))
myCuteColors <- c("Control" = "#FFC0CB", "Treat" = "#87CEFA")

boxplot_cute <- ggboxplot(
  data_long,
  x = "Immune", 
  y = "Fraction",
  fill = "Group",           
  palette = myCuteColors,  
  xlab = "",
  ylab = "Fraction",
  legend.title = "Group",
  notch = FALSE,
  width = 0.8
) +
  stat_compare_means(
    aes(group = Group),
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 1),  # 修正size为linewidth
    axis.ticks = element_line(color = "black", linewidth = 1)  # 修正size为linewidth
  ) +
  labs(
    title = "Immune Cell Comparison",
    subtitle = paste0(" (Control n=", countControl, ", Treat n=", countTreat, ")")
  )

ggsave("immune_diff-points_n.pdf", boxplot_cute, width = 8, height = 6)
cat("箱线图已保存为：immune_diff-points_n.pdf\n")

# -------- 导出表格2：统计汇总表 --------
summary_table <- data_long %>%
  group_by(Immune, Group) %>%
  summarise(
    MeanFraction = mean(Fraction, na.rm = TRUE),
    MedianFraction = median(Fraction, na.rm = TRUE),
    SD = sd(Fraction, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  )
write.csv(summary_table, file = "boxplot_summary_table.csv", row.names = FALSE)
cat("箱线图统计汇总表已保存为：boxplot_summary_table.csv\n")

# ------------------ 计算免疫细胞p值 ------------------
pvalue_table <- data_long %>%
  group_by(Immune) %>%
  do(tidy(wilcox.test(Fraction ~ Group, data = .))) %>%
  select(Immune, p.value)
write.csv(pvalue_table, file = "immune_pvalues.csv", row.names = FALSE)
cat("免疫细胞在两组间的 p 值已保存为：immune_pvalues.csv\n")

# ------------------ 单免疫细胞箱线图（SCI风格） ------------------
output_folder <- "immune_boxplots_SCI"
if(!dir.exists(output_folder)) {
  dir.create(output_folder)
  cat("已创建文件夹：", output_folder, "\n")
} else {
  cat("文件夹已存在：", output_folder, "\n")
}

immune_types <- unique(data_long$Immune)
for (immune_cell in immune_types) {
  data_subset <- subset(data_long, Immune == immune_cell)
  y_pos <- max(data_subset$Fraction, na.rm = TRUE) * 1.1
  if (y_pos == 0) y_pos <- 0.001
  
  p <- ggplot(data_subset, aes(x = Group, y = Fraction, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +
    geom_jitter(shape = 21, color = "black", alpha = 0.7, width = 0.15, size = 2) +
    scale_fill_npg() +
    stat_compare_means(aes(group = Group), label = "p.format", method = "wilcox.test",
                       label.x = 1.5, label.y = y_pos) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 14, hjust = 0.5, color = "gray30"),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = paste("Immune Cell:", immune_cell),
      subtitle = paste0("Control n=", sum(data_subset$Group == "Control"),
                        ", Treat n=", sum(data_subset$Group == "Treat")),
      x = NULL,
      y = "Fraction"
    )
  
  output_file <- file.path(output_folder, paste0(immune_cell, "_boxplot_SCI.pdf"))
  ggsave(output_file, p, width = 6, height = 8)
  cat("箱线图已保存为：", output_file, "\n")
}

# ------------------ 绘制相关性热图 ------------------
mat <- rt[, !(colnames(rt) %in% c("Group", "Sample"))]
pdf("immune_corrplot.pdf", width = 11, height = 11)
corrplot(cor(mat, method = "spearman", use = "pairwise.complete.obs"),
         method = "circle", type = "upper",
         tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.7)
dev.off()
cat("相关性热图已保存为：immune_corrplot.pdf\n")

# ------------------ 绘制两组叠加的山脊图 ------------------
p_ridges <- ggplot(data_long, aes(x = Fraction, y = Immune, fill = Group, color = Group)) +
  geom_density_ridges(alpha = 0.5, position = "identity", scale = 0.8, linewidth = 0.8) +  # 修正size为linewidth
  theme_minimal() +
  labs(
    title = "Ridgeline Plot of Immune Cells (Overlayed)",
    x = "Fraction",
    y = "Immune Cell"
  )
ggsave("immune_ridges_overlay.pdf", p_ridges, width = 10, height = 8)
cat("Overlayed Ridgeline Plot 已保存为：immune_ridges_overlay.pdf\n")


p_box_with_points <- ggplot(data_long, aes(x = Group, y = Fraction, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(color = Group), 
              position = position_jitter(width = 0.2), 
              size = 2, alpha = 0.7) +
  stat_compare_means(label = "p.format", method = "wilcox.test", 
                     label.y.npc = 0.8) +
  scale_fill_npg() +
  scale_color_npg() +
  theme_classic(base_size = 14) +
  facet_wrap(~ Immune, scales = "free") +
  labs(title = "Boxplot with Sample Distribution Points",
       x = "Group", y = "Fraction")

ggsave("boxplot_with_points.pdf", p_box_with_points, width = 13, height = 15)
cat("箱线图（含样本分布点）已保存为：boxplot_with_points.pdf\n")

# ------------------ 绘制箱线图（带样本点） ------------------
boxplot_with_points <- ggboxplot(
  data_long,
  x = "Immune", 
  y = "Fraction",
  fill = "Group",           
  palette = myCuteColors,  
  xlab = "",
  ylab = "Fraction",
  legend.title = "Group",
  notch = FALSE,
  width = 0.8
) +
  stat_compare_means(
    aes(group = Group),
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  ) +
  geom_jitter(shape = 21, color = "black", alpha = 0.7, width = 0.15, size = 2) +  # 添加样本点
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 1),  # 修正size为linewidth
    axis.ticks = element_line(color = "black", linewidth = 1)  # 修正size为linewidth
  ) +
  labs(
    title = "Immune Cell Comparison",
    subtitle = paste0(" (Control n=", countControl, ", Treat n=", countTreat, ")")
  )

ggsave("immune_diff-points_n_with_samples.pdf", boxplot_with_points, width = 8, height = 6)
cat("带样本点的箱线图已保存为：immune_diff-points_n_with_samples.pdf\n")