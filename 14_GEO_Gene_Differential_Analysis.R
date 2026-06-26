# 安装和加载必要的R包
# if (!requireNamespace("BiocManager", quietly = TRUE)) 
#     install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("ggpubr")

# 加载必要的R包
library(limma)   # 用于微阵列数据分析
library(ggpubr)  # 用于绘制美观的统计图形
library(pROC)    # 用于绘制和分析ROC曲线
library(ggplot2) # 用于图形的增强
library(gridExtra)
# 设置工作目录，替换为您本地的路径
# setwd("D:\\文件1-25和文件38，网络毒理+分子模拟\\21.诊断基因的差异分析和ROC曲线")

# 读取标准化后的表达矩阵文件
expFile = "geneexp.csv"  # 文件路径
rt = read.csv(expFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

# 提取样本的分组标签，并转换为0和1
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))  # 从列名提取分组标签
y = ifelse(y == "con", 0, 1)  # 对照组为0，实验组为1

# 读取基因列表文件
geneFile = "IntersectionGenes.csv"  # 文件路径
geneRT = read.csv(geneFile, header = FALSE, sep = ",", check.names = FALSE)

# 从基因列表中提取所需的基因
selectedGenes = as.vector(geneRT[, 1])  # 提取基因名
rt_filtered = rt[selectedGenes, , drop = FALSE]  # 仅保留基因列表中的基因

# 设置比较组（对照组和实验组）样本数量
conNum = sum(y == 0)
treatNum = sum(y == 1)
Type = c(rep("con", conNum), rep("treat", treatNum))  # 对照组与实验组标签

# 设置调色板
colors <- c("con" = "#FF6347", "treat" = "#4682B4")  # 红色和蓝色

# 创建一个数据框来存储基因的差异分析结果
result_df <- data.frame(
  Gene = character(),
  P_Value = numeric(),
  AUC = numeric(),
  AUC_Lower_CI = numeric(),
  AUC_Upper_CI = numeric(),
  SE = numeric(),
  Mean_Con = numeric(),
  Mean_Treat = numeric(),
  stringsAsFactors = FALSE
)

# 设置进度条来跟踪分析进度
cat("开始基因差异分析...\n")
total_genes = nrow(rt_filtered)
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)  # 初始化进度条

# 遍历基因列表中的每个基因，绘制箱线图并保存
for (i in rownames(rt_filtered)) {
  
  # 获取该基因的表达数据，并创建数据框
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  
  # 执行t检验
  t_test_result = t.test(as.numeric(rt_filtered[i, ]) ~ Type)
  
  # 提取t检验的p值和标准误（SE）
  p_value = t_test_result$p.value
  SE = t_test_result$stderr
  
  # 计算每组的均值
  mean_con = mean(rt1[rt1$Type == "con", "expression"])
  mean_treat = mean(rt1[rt1$Type == "treat", "expression"])
  
  # 计算并绘制ROC曲线
  roc1 = roc(y, as.numeric(rt_filtered[i, ]))
  ci1 = ci.auc(roc1, method = "bootstrap")  # 计算AUC的95%置信区间
  ciVec = as.numeric(ci1)
  
  # 将每个基因的统计结果添加到result_df中
  result_df = rbind(result_df, data.frame(
    Gene = i,
    P_Value = p_value,
    AUC = ciVec[2],
    AUC_Lower_CI = ciVec[1],
    AUC_Upper_CI = ciVec[3],
    SE = SE,
    Mean_Con = mean_con,
    Mean_Treat = mean_treat
  ))
  
  # 绘制箱线图并添加显著性比较
  boxplot = ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    geom_jitter(color = "black", size = 2, width = 0.2) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(i, "expression")) +
    stat_compare_means(method = "t.test")  # 显示统计比较
  
  # 保存箱线图到PDF文件
  pdf(file = paste0("boxplot.", i, ".pdf"), width = 3.4, height = 4.5)
  print(boxplot)
  dev.off()
  
  # 绘制并保存ROC曲线
  pdf(file = paste0("ROC.", i, ".pdf"), width = 5, height = 5)
  plot(roc1, print.auc = TRUE, col = "orange", legacy.axes = TRUE, main = i)  # 绘制ROC曲线
  
  # 添加AUC的置信区间阴影
  polygon(c(roc1$specificities, rev(roc1$specificities)), 
          c(roc1$sensitivities, rep(0, length(roc1$specificities))),
          col = rgb(1, 0.647, 0, 0.3), border = NA)  # 设置阴影的颜色和透明度
  
  # 在图中添加置信区间文本
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "orange")
  dev.off()
  
  # 绘制并保存密度图
  density_plot = ggplot(rt1, aes(x = expression, fill = Type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(title = paste(i, "Density Plot"), x = "Expression", y = "Density")
  
  # 保存密度图
  pdf(file = paste0("density_plot.", i, ".pdf"), width = 6, height = 6)
  print(density_plot)
  dev.off()
  
  # 更新进度条
  setTxtProgressBar(pb, which(rownames(rt_filtered) == i))  # 更新进度条
}

# 关闭进度条
close(pb)

# 将结果保存为CSV文件
write.csv(result_df, file = "gene_analysis_results.csv", row.names = FALSE)

# 输出结果表格
cat("差异分析结果已保存\n")



# 创建一个空的数据框来存储所有基因的ROC曲线
roc_plots <- list()
boxplot_plots <- list()

# 为每个基因生成箱线图和ROC曲线
for (i in rownames(rt_filtered)) {
  
  # 获取该基因的表达数据，并创建数据框
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  
  # 绘制箱线图
  boxplot = ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    geom_jitter(color = "black", size = 2, width = 0.2) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(i, "expression")) +
    stat_compare_means(method = "t.test")  # 显示统计比较
  
  boxplot_plots[[i]] <- boxplot  # 将箱线图添加到列表中
  
  # 绘制ROC曲线
  roc1 = roc(y, as.numeric(rt_filtered[i, ]))
  roc_plot = ggroc(roc1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    ggtitle(paste(i, "AUC = ", round(auc(roc1), 3))) +  # 显示AUC值
    theme(legend.position = "none")
  
  roc_plots[[i]] <- roc_plot  # 将ROC曲线添加到列表中
}

# 合并箱线图和ROC曲线
combined_plots <- list()

for (i in 1:length(boxplot_plots)) {
  combined_plots[[i]] <- grid.arrange(
    boxplot_plots[[i]], 
    roc_plots[[i]], 
    ncol = 2, 
    top = paste("Gene:", names(boxplot_plots)[i])  # 标题
  )
  
  # 保存每个基因的合并图
  ggsave(
    filename = file.path(output_folder, paste0("combined_", names(boxplot_plots)[i], ".pdf")),
    plot = combined_plots[[i]],
    width = 10,
    height = 6
  )
}

cat("所有基因的合并图已保存！\n")


# 设置比较组（对照组和实验组）样本数量
conNum = sum(y == 0)
treatNum = sum(y == 1)
Type = c(rep("con", conNum), rep("treat", treatNum))  # 对照组与实验组标签

# 设置调色板
colors <- c("con" = "#FF6347", "treat" = "#4682B4")  # 红色和蓝色

# 创建一个数据框来存储所有基因的箱线图数据
combined_data <- data.frame(expression = numeric(), Type = character(), Gene = character())

# 为每个基因生成箱线图数据
for (i in rownames(rt_filtered)) {
  # 获取该基因的表达数据，并创建数据框
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  rt1$Gene = i  # 为数据框添加基因列
  combined_data = rbind(combined_data, rt1)  # 将数据合并到一个数据框中
}

# 自定义p值的星号显示
pvalue_to_asterisk <- function(pvalue) {
  if (pvalue <= 0.001) {
    return("***")
  } else if (pvalue <= 0.01) {
    return("**")
  } else if (pvalue <= 0.05) {
    return("*")
  } else {
    return("ns")  # 非显著
  }
}

# 绘制所有基因的箱线图（合并在一起）
combined_boxplot = ggplot(combined_data, aes(x = Gene, y = expression, fill = Type)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
  scale_fill_manual(values = colors, labels = c("Control Group", "Treatment Group")) +  # 添加图例标签
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom",  # 图例放在底部
    axis.text.x = element_text(angle = 90, hjust = 1)  # 旋转基因名以便显示
  ) +
  labs(
    y = "Expression", 
    x = "Gene", 
    title = "Diagnostic Model Differential Analysis"  # 添加标题
  ) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",  # 使用显著性符号显示
    label.x = 1.5  # 将p值标签放在箱线图中央
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 为y轴增加一些间隙

# 保存合并后的箱线图
ggsave("combined_boxplot_all_genes_with_labels_and_title.pdf", plot = combined_boxplot, width = 9, height = 6)

cat("所有基因的合并箱线图已保存！\n")
## ========= 多基因合并到一张ROC曲线图 ========= ##
# 假定你用的基因均在selectedGenes中

# 用于存储ROC对象和AUC
roc_list <- list()
auc_list <- c()
leglab <- c()
my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
             "purple", "gray40", "black", "magenta", "gold", "brown")
while(length(my_cols) < length(selectedGenes)) {
  my_cols <- c(my_cols, rainbow(length(selectedGenes) - length(my_cols)))
}

pdf("All_Genes_Combined_ROC.pdf", width = 6, height = 6)
par(mar = c(5, 6, 4, 2)+0.1, cex = 1.3)

# 先画空图
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "All Genes Combined ROC", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45)
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

for(i in seq_along(selectedGenes)) {
  gene <- selectedGenes[i]
  gene_exp <- as.numeric(rt_filtered[gene, ])
  # 绘制ROC
  roc1 <- roc(y, gene_exp, direction="auto")
  auc1 <- as.numeric(auc(roc1))
  # 确保AUC反映区分而不是方向
  if (auc1 < 0.5) {
    auc1 <- 1-auc1
    roc1 <- roc(y, gene_exp, direction=ifelse(roc1$direction==">","<",">"))
  }
  lines(1 - roc1$specificities, roc1$sensitivities, col = my_cols[i], lwd = 2.5)
  leglab <- c(leglab, paste0(gene, "  AUC=", sprintf("%.3f", auc1)))
}
legend("bottomright", legend = leglab, col = my_cols[1:length(selectedGenes)], 
       lwd = 2.5, lty = 1, bty = "n", cex = 1)
dev.off()
cat("所有基因合并ROC图已保存 All_Genes_Combined_ROC.pdf\n")


print(result_df)  # 打印最终结果
