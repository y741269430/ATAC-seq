# 绘制deeptools的热图  
## 目录 ####
- 0.挑选基因集，构建gtf，用作绘制TSS富集热图展示     
- 1.bam转bw
- 2.computeMatrix 计算
- 3.绘图  

## 0.挑选基因集，构建gtf，用作绘制TSS富集热图展示     
```bash
vim geneList2gtf.sh
```
```bash
#!/bin/bash
#genefile 需要有geneID列！

genefile=Treatment_vs_CTRL.gene_p005.txt
gtf=/Homo_sapiens.GRCh38.112.chr.gtf

awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="geneId") col=i} NR>1 {print $col}' "$genefile" > gene_list.txt
grep -Fwf gene_list.txt "$gtf" > filtered.gtf
```
```bash
bash geneList2gtf.sh
```

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
nohup bamCoverage --bam mergebam/CTRL.bam -o mergebam/CTRL.bw --binSize 20 --normalizeUsing RPKM & 
nohup bamCoverage --bam mergebam/Treatment.bam -o mergebam/Treatment.bw --binSize 20 --normalizeUsing RPKM &

nohup bamCoverage --bam mergebam/CTRL.bam -o mergebam/CTRL.bw --binSize 10 --effectiveGenomeSize 2723414844 --normalizeUsing RPGC --outFileFormat bigwig -p 10 &
nohup bamCoverage --bam mergebam/Treatment.bam -o mergebam/Treatment.bw --binSize 10 --effectiveGenomeSize 2723414844 --normalizeUsing RPGC --outFileFormat bigwig -p 10 &

# 顺便把peak也call了，与单独样本call的peak，一同放在igv查看(这里的peak文件也可以接着做homer的motif预测)
nohup macs3 callpeak -f BAMPE -t mergebam/CTRL.bam -g mm -n mergebam/CTRL -B -q 0.1 &  
nohup macs3 callpeak -f BAMPE -t mergebam/Treatment.bam -g mm -n mergebam/Treatment -B -q 0.1 &  

```
- [homer传送门](https://github.com/y741269430/homer)    

    
#### bamCoverage里的 --binSize 20是什么意思？     
​​--binSize 20​​ 表示：​​将基因组分成 20 bp 的窗口（bin），并计算每个 20 bp 区间内的平均测序深度（reads 覆盖度）​​。    
​​bamCoverage的作用​​ 是把 ​​BAM 文件（原始测序比对结果）转换成 BigWig 文件（信号强度文件）​​，而 --binSize决定了 ​​BigWig 文件里信号的分辨率​​。    
​​20 bp 的 bin​​ 意味着：​​每 20 个碱基计算一次覆盖度​​，最终 BigWig 文件里的信号是 ​​20 bp 一个点​​。    
如果改成 --binSize 10，就是 ​​10 bp 一个点​​（信号更精细，但文件更大）；    
如果改成 --binSize 50或 100，就是 ​​50 bp 或 100 bp 一个点​​（信号更平滑，但可能丢失细节）。    

#### bamCoverage里面为什么需要做归一化，做哪种归一化？    
归一化是为了让 BigWig 文件里的信号更稳定​​，避免因为测序深度差异导致信号忽高忽低。不同归一化的区别如下，目前默认根据基因组大小，选择`--normalizeTo1x`进行归一化：       

| **归一化方式**               | **参数**                          | **作用**                                                                 | **适用场景**                                                                 | **示例命令**                                                                 |
|------------------------------|-----------------------------------|--------------------------------------------------------------------------|------------------------------------------------------------------------------|------------------------------------------------------------------------------|
| **1x 基因组覆盖度标准化**     | `--normalizeTo1x <基因组大小>`    | 将信号标准化到“30亿条reads均匀覆盖基因组”的水平（如人类基因组约3Gbp，参数设为`3000000000`） | 比较不同测序深度的样本（如CTRL vs HR），消除测序量差异影响                     | `bamCoverage --bam input.bam -o output.bw --binSize 10 --normalizeTo1x 3000000000` |
| **无归一化（原始信号）**      | （不添加任何归一化参数）           | 直接输出原始reads覆盖度，信号值受测序深度直接影响                           | 仅用于单样本分析或测序深度相同的样本比较                                       | `bamCoverage --bam input.bam -o output.bw --binSize 10`                        |
| **RPKM/CPM归一化**            | `--normalizeUsing RPKM` 或 `CPM`  | 按“每千碱基每百万reads”（RPKM）或“每百万reads”（CPM）标准化                | 适用于需要与其他RNA-seq等转录组数据比较的场景（但ChIP-seq较少用）              | `bamCoverage --bam input.bam -o output.bw --binSize 10 --normalizeUsing RPKM`  |
| **BPM（百分位数中位数标准化）**| `--normalizeUsing BPM`            | 按百分位数调整信号分布，减少极端高/低信号的影响                              | 信号分布不均匀时（如某些区域reads堆积严重）                                    | `bamCoverage --bam input.bam -o output.bw --binSize 10 --normalizeUsing BPM`   |
关键说明：    
1. **`--normalizeTo1x <基因组大小>`** 是ChIP-seq/ATAC-seq分析的**最常用归一化方法**（尤其用于多样本比较），例如人类基因组参数设为 `3000000000`（3Gbp）。    
2. 若**不指定归一化参数**，输出的BigWig文件信号值直接反映原始测序深度，不同样本间不可直接比较。    
3. 其他归一化方式（如RPKM/BPM/SES）适用于特定场景，需根据实验设计选择。    

## 2.computeMatrix 计算    
参考：  
- [computeMatrix](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)  
- [deepTools 使用指南](https://www.jianshu.com/p/2b386dd437d3)
- [使用deeptools生成ChIP-seq信号热图与谱图](https://www.jianshu.com/p/1f38674b2475)

<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/computeMatrix_modes.png" width="600" />

### single input files (reference-point mode)    
```bash
nohup computeMatrix reference-point \
--referencePoint TSS \
-S CTRL.bw \
-R /home/jjyang/downloads/genome/mm39_GRCm39/mm_gene_vM27.bed \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
-p 10 \
--skipZeros \
 --missingDataAsZero \
--binSize 20  \
-o matrix_TSS_CTRL.gz &

nohup computeMatrix reference-point \
--referencePoint TSS \
-S Treatment.bw \
-R /home/jjyang/downloads/genome/mm39_GRCm39/mm_gene_vM27.bed \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
-p 10 \
--skipZeros \
 --missingDataAsZero \
--binSize 20  \
-o matrix_TSS_Treatment.gz &

nohup computeMatrix reference-point \
--referencePoint TSS \
-S CTRL.bw Treatment.bw \
-R Homo_sapiens.GRCh38.112.chr.gtf \
--samplesLabel CTRL Treatment \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
-p 10 \
--skipZeros \
--missingDataAsZero \
--binSize 20 \
-o matrix_bin20.gz &

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
plotHeatmap -m matrix_TSS_CTRL.gz --colorList 'white,white,red' -out TSS_CTRL.pdf
plotHeatmap -m matrix_TSS_CTRL.gz -out TSS_CTRL.pdf

plotHeatmap -m matrix_bin20.gz --colorList 'white,white,red' -out TSS_TreVSCtrl.pdf
plotHeatmap -m matrix_bin20.gz -out TSS_TreVSCtrl.pdf
plotHeatmap -m matrix_bin20.gz --zMin 0 --zMax 3 -out TSS_TreVSCtrl.pdf
```  




















