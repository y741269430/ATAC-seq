# 绘制deeptools的热图  
## 目录 ####
- 0.gtf转bed
- 1.bam转bw
- 2.computeMatrix 计算
- 3.绘图  

## 0.gtf转bed
```bash
cat gencode.vM27.annotation.gtf | awk -F '[\t *;]' '/^chr/{if($3=="transcript"){print $1,$4-1,$5,$19,$22,$7,$3}}' OFS="\t" > mm_gene_vM27.bed
head mm_gene_vM27.bed
```
它是这样的  
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/mm_gene_vM27bed_image.png" width="600" />

- 假如我有特定的基因，那么先构建一个数据框
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/test_csv_bed_image.png" width="600" />
然后生成该数据框的bed

```bash
awk -F',' 'NR==1 {for(i=1;i<=NF;i++) if($i ~ /SYMBOL/) col=i} NR>1 {print $col}' test.csv | grep -Ff - mm_gene_vM27.bed > test.bed
head test.bed
```
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/test_csv_bed_image2.png" width="600" />

然后就可以用上述bed文件去绘制热图   


## 1.bam转bw
前面提到了利用deeptools将bam文件转为bw文件 ([计算重复样本的 peak 之间的重叠的坐标位置.md](https://github.com/y741269430/ATAC-seq/blob/main/%E8%AE%A1%E7%AE%97%E9%87%8D%E5%A4%8D%E6%A0%B7%E6%9C%AC%E7%9A%84%20peak%20%E4%B9%8B%E9%97%B4%E7%9A%84%E9%87%8D%E5%8F%A0%E7%9A%84%E5%9D%90%E6%A0%87%E4%BD%8D%E7%BD%AE.md#5%E4%BD%BF%E7%94%A8deeptools-%E5%B0%86bam%E8%BD%AC%E4%B8%BAbw%E7%94%A8%E4%BA%8Eigv%E5%8F%AF%E8%A7%86%E5%8C%96))

举例，我们这里将bam文件进行合并，然后再转为bw  

```bash
# 合并前先确保bam文件是position排序而不是name排序！
nohup samtools sort -@ 10 -O bam CTRL-sorted-name.bam -o CTRL_1.last.bam &
nohup samtools sort -@ 10 -O bam Treatment-sorted-name.bam -o Treatment_1.last.bam &

# 相同样本合并bam
nohup samtools merge -@ 10 mergebam/CTRL.bam CTRL_1.last.bam CTRL_2.last.bam &
nohup samtools merge -@ 10 mergebam/Treatment.bam Treatment_1.last.bam Treatment_2.last.bam &

# 构建index
nohup samtools index -@ 10 mergebam/CTRL.bam &
nohup samtools index -@ 10 mergebam/Treatment.bam &

# deeptools的bam转bw
nohup bamCoverage --bam mergebam/CTRL.bam -o mergebam/CTRL.bw --binSize 10 --normalizeUsing RPKM & 
nohup bamCoverage --bam mergebam/Treatment.bam -o mergebam/Treatment.bw --binSize 10 --normalizeUsing RPKM &

# 顺便把peak也call了，与单独样本call的peak，一同放在igv查看(这里的peak文件也可以接着做homer的motif预测)
nohup macs3 callpeak -f BAMPE -t mergebam/CTRL.bam -g mm -n mergebam/CTRL -B -q 0.1 &  
nohup macs3 callpeak -f BAMPE -t mergebam/Treatment.bam -g mm -n mergebam/Treatment -B -q 0.1 &  

```

- [homer传送门](https://github.com/y741269430/homer)

## 2.computeMatrix 计算   
参考：  
- [computeMatrix](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)  
- [deepTools 使用指南](https://www.jianshu.com/p/2b386dd437d3)
- [使用deeptools生成ChIP-seq信号热图与谱图](https://www.jianshu.com/p/1f38674b2475)

<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/computeMatrix_modes.png" width="600" />

### single input files (reference-point mode)    
```bash
nohup computeMatrix reference-point \
-S CTRL.bw \
-R /home/jjyang/downloads/genome/mm39_GRCm39/mm_gene_vM27.bed \
--referencePoint TSS \
-a 3000 -b 3000 \
-p 10 \
--skipZeros \
-o matrix_TSS_CTRL.gz &

nohup computeMatrix reference-point \
-S Treatment.bw \
-R /home/jjyang/downloads/genome/mm39_GRCm39/mm_gene_vM27.bed \
--referencePoint TSS \
-a 3000 -b 3000 \
-p 10 \
--skipZeros \
-o matrix_TSS_Treatment.gz &
```

### multiple input files (scale-regions mode)   
```bash
nohup computeMatrix scale-regions -S \
CTRL.bw \
Treatment.bw  \
-R /home/jjyang/downloads/genome/mm39_GRCm39/mm_gene_vM27.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
-p 10 \
--skipZeros \
-o matrix.test.gz &
```
## 3.绘图  
参考： 
- [plotProfile](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html)   
```bash
plotProfile -m matrix.test.gz \
-out profile_split.png \
--numPlotsPerRow 1 \
--plotTitle "Test data profile"

plotProfile -m matrix.test.gz \
-out profile_group.png \
--plotType=fill \
--perGroup \
--plotTitle "Test data profile"
```
参考：  
- [plotHeatmap](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html)   
```bash  
plotHeatmap -m matrix.test.gz -out ExampleHeatmap.pdf

plotHeatmap -m matrix_TSS_CTRL.gz -out Heatmap_CTRL.png
plotHeatmap -m matrix_TSS_Treatment.gz -out Heatmap_Treatment.png
```  




















