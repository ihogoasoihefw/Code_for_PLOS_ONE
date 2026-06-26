# 如果尚未安装，请取消注释以下安装包的代码
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("dplyr")
# install.packages("broom")

# ------------------ 加载所需包 ------------------
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)

# ------------------ 数据准备 ------------------
inputFile <- "CIBERSORT-Results.csv"
# setwd("D:\\文件1-25和文件38，网络毒理+分子模拟\\23.免疫浸润可视化1.热图")

# 读取数据（假设每行代表一个样本，列为各免疫细胞的比例）
rt <- read.table(inputFile, header = TRUE, sep = ",", 
                 check.names = FALSE, row.names = 1)

# 根据行名判断组别：后缀为 _con 的为 Control，后缀为 _tre 的为 Treat
rt$Group <- ifelse(grepl("_con$", rownames(rt), ignore.case = TRUE), "Control",
                   ifelse(grepl("_tre$", rownames(rt), ignore.case = TRUE), "Treat", NA))


# 获取 Control 与 Treat 的样本名，并构造样本顺序：Control -> gap -> Treat
control_samples <- rownames(rt)[rt$Group == "Control"]
treat_samples   <- rownames(rt)[rt$Group == "Treat"]
all_samples_ordered <- c(control_samples, "gap", treat_samples)

# ------------------ 数据转换 ------------------
rt$Sample <- rownames(rt)
data_long <- melt(rt, id.vars = c("Sample", "Group"), 
                  variable.name = "Immune", 
                  value.name = "Fraction")
# 指定因子水平，确保 gap 保留
data_long$Sample <- factor(data_long$Sample, levels = all_samples_ordered)

# -------- 导出表格1：绘制堆叠条形图使用的长格式数据 --------
write.csv(data_long, file = "barplot_data_long.csv", row.names = FALSE)
cat("Barplot 长格式数据已保存为：barplot_data_long.csv\n")

# ------------------ 使用 Set3 调色板 ------------------
immune_types <- unique(data_long$Immune)
nColors <- length(immune_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

# ------------------ 绘制堆叠条形图 ------------------
barplot_with_gap <- ggplot(data_long, aes(x = Sample, y = Fraction, fill = Immune)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myColors) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = "Relative Percent", 
    fill = "Immune\nCell Type",
    title = "Immune Cell Distribution",
    subtitle = "Control vs. Treat "
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)))

# 计算对照组和实验组的样本数（不包含 gap）
control_count <- length(control_samples)
treat_count   <- length(treat_samples)

barplot_with_gap_annot <- barplot_with_gap +
  # 添加 Control 组下方的线段和文字（将size改为linewidth）
  annotate("segment", x = 0.5, xend = control_count + 0.5, 
           y = -0.04, yend = -0.04, color = "#D65DB1", linewidth = 5) +
  annotate("text", x = (control_count)/2 + 0.5, y = -0.08, 
           label = "Control", color = "#D65DB1", size = 7, fontface = "bold") +
  # 添加 Treat 组下方的线段和文字（将size改为linewidth）
  annotate("segment", x = control_count + 1.5, 
           xend = control_count + treat_count + 1.5, 
           y = -0.04, yend = -0.04, color = "#0089BA", linewidth = 5) +
  annotate("text", x = control_count + (treat_count)/2 + 1.5, y = -0.08, 
           label = "Treat", color = "#0089BA", size = 7, fontface = "bold")

ggsave("barplot_two_lines_set3_gap.pdf", barplot_with_gap_annot, 
       width = 12, height = 7.5)
cat("堆叠条形图已保存为：barplot_two_lines_set3_gap.pdf\n")
