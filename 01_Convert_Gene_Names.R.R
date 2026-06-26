# 设置工作目录（建议用正斜杠）
# working_dir <- "E:/网络毒理学最新课程/10.【纯生信无转录组】Q1区6.1，网络毒理学+毒性分析+孟德尔随机化+分子对接+分子动力学+单细胞/05.CHEMBL数据库中检索化合物的潜在靶点"
# setwd(working_dir)

# 加载包
library(data.table)
library(stringr)

# 读取数据（核心修正：列名是Accessions，不是UniProt Accessions）
chembl_file <- "gene.tsv"
ref_file <- "UniProt gene conversion reference file.tsv"

chembl <- fread(chembl_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
ref <- fread(ref_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# 验证列名（确认Accessions列存在）
cat("=== 验证chembl数据列名 ===\n")
print(colnames(chembl))

# 提取First_UniProt（修正列名，处理空值和多ID情况）
get_first_uniprot <- function(x) {
  if (is.na(x) || x == "" || trimws(x) == "") {
    return(NA_character_)
  } else {
    ids <- str_split(trimws(x), "\\|")[[1]]
    ids <- ids[ids != ""]
    return(ids[1])
  }
}

chembl$First_UniProt <- sapply(chembl$Accessions, get_first_uniprot)

# 查看提取结果（修正列名引用：`ChEMBL ID`）
cat("\n=== 提取后的First_UniProt前10行 ===\n")
# 关键修正：用`ChEMBL ID`代替ChEMBL.ID
print(head(chembl[, .(`ChEMBL ID`, Accessions, First_UniProt)], 10))

cat("\n=== 有效First_UniProt数量 ===\n")
print(sum(!is.na(chembl$First_UniProt)))

# 去空格，准备合并
chembl$First_UniProt <- trimws(chembl$First_UniProt)
ref$Entry <- trimws(ref$Entry)

# 合并参考文件（按First_UniProt和Entry匹配）
result <- merge(
  chembl,
  ref[, .(Entry, Gene_Name=`Gene Names`)],
  by.x="First_UniProt",
  by.y="Entry",
  all.x=TRUE,
  sort=FALSE
)

# 查看匹配结果
cat("\n=== Gene_Name匹配情况 ===\n")
match_stats <- table(is.na(result$Gene_Name))
print(match_stats)
cat(sprintf("匹配成功率：%.2f%%\n", (match_stats[["FALSE"]]/nrow(result))*100))

# 导出结果
fwrite(result, "gene_with_genename_corrected.tsv", sep="\t", na="")
cat("\n结果已导出为：gene_with_genename_corrected.tsv\n")

# 查看匹配成功的示例数据（修正列名引用：`ChEMBL ID`）
cat("\n=== 匹配成功的前5行示例 ===\n")
# 关键修正：用`ChEMBL ID`代替ChEMBL.ID
print(head(result[!is.na(result$Gene_Name), .(`ChEMBL ID`, First_UniProt, Name, Gene_Name)], 5))