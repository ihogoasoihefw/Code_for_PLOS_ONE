# setwd("C:\\Users\\wode3\\Desktop\\20250507\\32.Merging ligand and receptor data using R")

# 文件路径
gro1 <- "protein_receptor_processed.gro"  # 受体
gro2 <- "ligands.gro"      # 配体
out_gro <- "complex.gro"

# 读取文件内容
lines1 <- readLines(gro1)
lines2 <- readLines(gro2)

# 解析原子数
n1 <- as.integer(lines1[2])
n2 <- as.integer(lines2[2])

# 原子坐标行（第3行到倒数第2行）
atoms1 <- lines1[3:(2 + n1)]
atoms2 <- lines2[3:(2 + n2)]

# 盒子尺寸（最后一行，通常用受体的盒子尺寸）
box_line <- lines1[length(lines1)]

# 合并
all_atoms <- c(atoms1, atoms2)
total_atoms <- length(all_atoms)

# 新gro文件内容
new_gro <- c(
  "Complex: receptor + ligand",   # 新标题
  as.character(total_atoms),      # 新原子数
  all_atoms,                     # 所有原子坐标
  box_line                       # 盒子尺寸
)

# 写入新文件
writeLines(new_gro, out_gro)
cat("已合并gro文件，输出为：", out_gro, "\n")



# Step 1: 读取两个文件内容
lig_content <- readLines("ligands.itp")        # 读取lig.itp文件
topol_content <- readLines("protein_receptor.top")    # 读取topol.top文件

# Step 2: 合并内容
merged_content <- c( topol_content,lig_content)

# Step 3: 写入合并后的文件，文件名以.top为后缀
writeLines(merged_content, "merged_top.top")

