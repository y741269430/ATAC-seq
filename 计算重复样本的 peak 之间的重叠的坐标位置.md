# 计算重复样本的 peak 之间的重叠的坐标位置   

## 目录 #### 
- 1.IDR 计算peaks之间的overlaping  
- 2.bedtools 计算peaks之间的overlaping，输出bed文件  
- 3.intervene 计算peaks之间的overlaping  
- 4.deeptools 计算peaks之间的overlaping和correlation  
- 5.使用deeptools 将bam转为bw（用于igv可视化）
- 6.deeptools 计算bam PE FragmentSize 统计片段长度
- 7.利用wigCorrelate计算样本相关性    

---
## 1.IDR 计算peaks之间的overlaping  
参考  
- [07_handling-replicates-idr.md](https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lessons/07_handling-replicates-idr.md)
- [08_handling-replicates](https://rkhetani.github.io/Intro-to-ChIPseq/08_handling-replicates)
- https://www.jianshu.com/p/d8a7056b4294

>To run IDR the recommendation is to run MACS2 less stringently to allow a larger set of peaks to be identified for each replicate. In addition the narrowPeak files have to be sorted by the `-log10(p-value)` column.   

对peak的log10pvalue进行排序
```bash
cat filenames | while read i; do sort -k8,8nr macs3/${i}_peaks.narrowPeak > macs3/sorted_${i}_peaks.narrowPeak; done &
```
创建环境
```bash
conda create -n idr
conda activate idr
conda install -c bioconda idr

idr --samples ../macs3/a1_rep1_peaks.narrowPeak ../macs3/b1_rep2_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./a1_b1_overlaps \
--plot \
--log-output-file a1_b1_overlaps.log
```

The output file format mimics the input file type, with some additional fields. Note that the **first 10 columns are a standard narrowPeak file**, pertaining to the merged peak across the two replicates. 

**Column 5 contains the scaled IDR value, `min(int(log2(-125IDR), 1000)`** For example, peaks with an IDR of 0 have a score of 1000, peaks with an IDR of 0.05 have a score of int(-125log2(0.05)) = 540, and IDR of 1.0 has a score of 0.

**Columns 11 and 12 correspond to the local and global IDR value, respectively.** 
* The **global IDR** is the value used to calculate the scaled IDR number in column 5, it _is analogous to a multiple hypothesis correction on a p-value to compute an FDR_. 
* The **local IDR** is akin to the posterior probability of a peak belonging to the irreproducible noise component. You can read [this paper](http://projecteuclid.org/euclid.aoas/1318514284
) for more details. 

**Columns 13 through 16 correspond to Replicate 1 peak data** and **Columns 17 through 20 correspond to Replicate 2 peak data.**

More detail on the output can be [found in the user manual](https://github.com/nboley/idr#output-file-format). Also, if you have any unanswered questions check out posts in the [Google groups forum](https://groups.google.com/forum/#!forum/idr-discuss).  

## 2.bedtools 计算peaks之间的overlaping，输出bed文件   
具体参考   
- [07_handling_peaks_bedtools.md](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/07_handling_peaks_bedtools.md)
- https://www.jianshu.com/p/f8bbd51b5199  

- `-wo`: Write the original A (file 1) and B (file 2) entries plus the number of base pairs of overlap between the two features.  
- `-f`: Minimum overlap required as a fraction of A. The value ranges from 0 to 1. We will use 0.3, requiring the overlap region being at least 30% of A.  
- `-r`: Require that the fraction overlap be reciprocal for A and B. Together with the `-f` flag above, we require the overlap region being at least 30% of B as well.  

取narrowPeak的前三列（计算peak的重叠区域仅仅是计算其位置的重复区间，跟里面peak的高低，形状之类的毫不相干，所以只取坐标位置即可）  
```bash
cat filenames | while read i; do cut -f1-3 macs3/${i}_peaks.narrowPeak > macs3/top3_${i}_peaks.narrowPeak; done &
```
```bash
bedtools intersect \
-wo -f 0.3 -r \
-a ../macs3/a1_peaks.narrowPeak \
-b ../macs3/b1_peaks.narrowPeak \
-wo > a1b1_overlaps.bed
```
> ## Other approaches for assessing peak reproducibility
> Historically, the ENCODE standard was using the overlaps that we described above but with a set of given criteria. This was developed based on experience with accumulated ENCODE ChIP-seq data, albeit with a much smaller sample size back then. In the paper [Landt et al, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/) describe the approach as:
> 
> _"...either 80% of the top 40% of the peaks identified from one replicate using an acceptable scoring method should overlap the list of peaks from the other replicate, OR peak lists scored using all available reads from each replicate should share more than 75% of regions in common."_ 
> 
> Since then, the field has moved towards more statistically motivated approaches like the [Irreproducibility Discovery Rate (IDR)](https://sites.google.com/site/anshulkundaje/projects/idr). The IDR framework was developed by Qunhua Li and Peter Bickel's group. It compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility. You can read more about IDR and how it works in this [linked lesson](handling-replicates-idr.md).
> 
> IDR analysis is extensively used by the ENCODE and modENCODE projects and is part of their ChIP-seq guidelines and standards. However, more recently there has been dicussion about the two approaches converging on similar results and so it remains to be seen what the gold standard will be.


## 3.intervene 计算peaks之间的overlaping  
```bash
conda create -n intervene
conda activate intervene
conda install -c bioconda intervene

nohup intervene venn  -i ../macs3/narrow/*.narrowPeak --save-overlaps &
nohup intervene upset -i ../m1/*_rmBL.narrowPeak --output ./ &
```

## 4.deeptools 计算peaks之间的overlaping和correlation  
```bash
conda create -n deeptools
conda activate deeptools
conda install -c bioconda deeptools

nohup multiBamSummary bins --bamfiles bam/*last.bam --minMappingQuality 30 --labels CTRL_1 CTRL_2 CTRL_3 -out readCounts.npz --outRawCounts readCounts.tab && 
  
nohup plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --log1p --removeOutliers -p scatterplot -o scatterplot_SpM.pdf --outFileCorMatrix Spearman.tab &
```

## 5.使用deeptools 将bam转为bw（用于igv可视化）
```bash
vim g3_bam2bw.sh

#!/bin/bash
## bam to bw ##

cat filenames | while read i; 
do
nohup bamCoverage --bam ./bam/${i}.last.bam -o ./bw/${i}.bw --binSize 10 --normalizeUsing RPKM & 
# nohup bamCoverage --bam ./bam/${i}.last.bam -o ./bw2/${i}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads & 
done
```

## 可视化  
```bash
nohup multiBigwigSummary bins -b *.bw -o test.npz && plotCorrelation -in test.npz --corMethod spearman --skipZeros --log1p --removeOutliers -p scatterplot -o scatterplot_SpM.pdf --outFileCorMatrix Spearman.tab &
```
    
## 6.deeptools 计算bam PE FragmentSize 统计片段长度  
```bash
nohup bamPEFragmentSize -hist fragmentSize_CON.png \
-T "Fragment size of CON" \
--maxFragmentLength 1000 \
-b CTRL_1.last.bam CTRL_2.last.bam \
--samplesLabel CTRL_1 CTRL_2 &
```

## 7.利用wigCorrelate计算样本相关性   
下载对应的程序，或者用conda 环境
```bash  
mamba install ucsc-wigcorrelate

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
chmod +x wigCorrelate
```
执行   
```bash
~/downloads/wigCorrelate peak1.bw peak2.bw
```








