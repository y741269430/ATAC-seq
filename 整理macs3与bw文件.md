# 整理macs3与bw文件
将macs3的分析结果，以及bw文件，整理到一个目标目录，进行下游分析    
## 1. 批量创建文件夹，并复制narrowPeak文件以及bw文件到对应的目录    
```bash
vim atac1_mkdirNmv.sh
```
```bash
#!/bin/bash

# 设置变量
group_file="Group_Info.txt"
source_dir_macs3="macs3"
source_dir_bw="bw"
target_base_peaks="Output/04.Peak_calling"
target_base_qc="Output/03.QC"

# 创建所有需要的目标基础目录
mkdir -p "$target_base_peaks"
mkdir -p "$target_base_qc"
mkdir -p Output/05.Peak_anno Output/06.Motif Output/07.GO_KEGG

# 读取文件并处理MACS3 peak文件
echo "开始处理MACS3 peak文件..."
while IFS=$'\t' read -r sample_id group_id; do
    # 跳过注释行和空行
    [[ "$sample_id" =~ ^# ]] && continue
    [[ -z "$sample_id" ]] && continue
    
    echo "处理样本: $sample_id (组: $group_id) - MACS3 peak文件"
    
    # 创建样本文件夹
    sample_dir="$target_base_peaks/$sample_id"
    mkdir -p "$sample_dir"
    sample_dir2="Output/05.Peak_anno/$sample_id"
    mkdir -p "$sample_dir2"

    # 源文件路径
    source_file="$source_dir_macs3/${sample_id}_peaks.narrowPeak"
    
    # 检查源文件是否存在
    if [[ -f "$source_file" ]]; then
        # 复制文件
        cp "$source_file" "$sample_dir/"
        echo " ✅ 已复制: $(basename "$source_file") -> $sample_dir/"
    else
        echo " ❌ 未找到: $source_file"
    fi
    echo "  ---"
done < "$group_file"

# 读取文件并处理bigWig文件
echo "开始处理bigWig文件..."
while IFS=$'\t' read -r sample_id group_id; do
    # 跳过注释行和空行
    [[ "$sample_id" =~ ^# ]] && continue
    [[ -z "$sample_id" ]] && continue
    
    echo "处理样本: $sample_id (组: $group_id) - bigWig文件"
    
    # 创建样本文件夹
    sample_dir="$target_base_qc/$sample_id"
    mkdir -p "$sample_dir"

    # 源文件路径
    source_file="$source_dir_bw/${sample_id}.bw"
    
    # 检查源文件是否存在
    if [[ -f "$source_file" ]]; then
        # 复制文件
        cp "$source_file" "$sample_dir/"
        echo " ✅ 已复制: $(basename "$source_file") -> $sample_dir/"
    else
        echo " ❌ 未找到: $source_file"
    fi
    echo "  ---"
done < "$group_file"

echo "所有处理完成！"
```
```bash
nohup bash atac1_mkdirNmv.sh >log1.txt 2>&1 &
```

## 2. 批量修改narrowPeak文件，去除线粒体行以及更改其他XY染色体名称，去除blacklist        
```bash
vim atac2_fixed_chr.sh
```
```bash
#!/bin/bash

# 定义黑名单文件路径
BLACKLIST="/home/jjyang/downloads/genome/mm39_GRCm39/mm10_blacklist.bed"

for file in Output/04.Peak_calling/*/*_peaks.narrowPeak; do

    # 提取文件所在目录和样本名
    file_dir=$(dirname "$file")
    file_basename=$(basename "$file")  # 如：Sample1_peaks.narrowPeak
    sample_name=${file_basename%%_peaks.narrowPeak}  # 如：Sample1

    echo " 正在处理样本: $sample_name, 文件: $file"

    # Step 1: 使用 awk 清理染色体信息并生成临时文件
    tmp_file="${file%.*}_peaks.narrowPeak_peaks.tmp.narrowPeak" 
    awk -F'\t' -v OFS='\t' '
        {
            chrom = $1
            if (chrom == "MT") next  
            if (chrom == "X") $1 = "chrX"
            else if (chrom == "Y") $1 = "chrY"
            print
        }' "$file" > "$tmp_file"

    # Step 2: 使用 bedtools intersect -v 去除黑名单区域
    output_file="${file_dir}/${sample_name}_peaks.noBlacklist.narrowPeak"

    echo " 正在使用 bedtools 过滤黑名单..."
    srun -p all -c 10 singularity exec \
        --bind /public1,/public23,/mnt,/prog1,/public1/Softwares \
        /prog1/Container/ACGT101_ATAC_seq/sif/ACGT101_ATAC_seq_1.0.sif \
        bedtools intersect -v -a "$tmp_file" -b "$BLACKLIST" > "$output_file"

    # 删除临时文件
    rm -f "$tmp_file"
    echo " 已删除临时文件: $tmp_file"
    echo "  ---"
done

echo "所有样本处理完成！"
```
```bash
nohup bash atac2_fixed_chr.sh >log2.txt 2>&1 &
```
