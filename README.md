# ATAC-seq workflow

本文主要参考：  
- [ATAC-seq-data-analysis-from-FASTQ-to-peaks](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/)
- [Intro-to-ChIPseq-flipped](https://github.com/hbctraining/Intro-to-ChIPseq-flipped)
- [Analytical Approaches for ATAC-seq Data Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8191135/)
- [From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)
- https://nf-co.re/chipseq/2.0.0/
- [ATAC-seq_QC_analysis](https://github.com/Zhang-lab/ATAC-seq_QC_analysis)

## 目录 ####
- 0.创建conda环境用于ATACseq分析（也可以用mamba）
- 0.利用bowtie2构建小鼠基因组（mm39）索引（构建一次以后都不用做了）  
- 1.开始——激活conda环境
- 2.利用Trim adaptors去除接头
- 3.比对到mm39 
- 4.生成raw bam (optional) 
- 5.sam to bam 同时去除 ChrM   
- 6.macs2（现已弃用）  
- 7.使用macs3进行call peak   
- narrowPeak和bed文件格式   
- 8.macs3 peak文件转 bw（用于igv可视化） 
- Remove blacklist  
- fastqc质控  
- louvain 聚类  
 
---
## 0.创建conda环境用于ATACseq分析（也可以用mamba）
我这里创建了多个环境，防止软件之间的冲突   

```bash
# 1
conda create -n atac
conda activate atac
mamba install samtools
mamba install bowtie2
mamba install sambamba
mamba install bedtools
mamba install macs3
mamba install trim-galore
mamba install bioconda::ucsc-bedclip
mamba install bioconda::ucsc-bedgraphtobigwig

# 2
conda create -n mqc python=3.6
conda activate mqc
pip install multiqc

# 3
conda create -n idr
conda activate idr
conda install -c bioconda idr

# 4
```
## 0.利用bowtie2构建小鼠基因组（mm39）索引（构建一次以后都不用做了）  

```bash
conda activate atac  
cd /home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/  
nohup bowtie2-build GRCm38.primary_assembly.genome.fa /home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39 &
```

## 1.开始——激活conda环境

创建文件夹   
```bash
conda activate atac  
mkdir -p raw clean trim bam macs2 macs2/narrow macs3 macs3/narrow
```
    
生成一个filenames的文件，用来记录输出的文件名称（样本名称），例如：  
```bash
ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 > filenames
```

## 2.利用Trim adaptors去除接头
```bash
vim pre_trim.sh

#!/bin/bash
## trim_galore ##

cat filenames | while read i; 
do
# paired end
nohup trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 1 --paired ./raw/${i}*_1.fq.gz ./raw/${i}*_2.fq.gz -o ./trim &
  
# single end
# nohup trim_galore -q 25 --phred33 --length 20 -e 0. 1 --stringency 1 ./raw/${i}*_1.fq.gz -o ./trim &
done
```

## 3.比对到mm39 
```bash
vim atac1_bw2.sh

#!/bin/bash
## Alignment to mm39 ##

mm39="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

cat filenames | while read i; 
do
nohup bowtie2 -p 4 --very-sensitive -X 2000 -k 10 \
-x ${mm39} \
-1 trim/${i}*_val_1.fq.gz \
-2 trim/${i}*_val_2.fq.gz \
-S ./bam/${i}.sam 2> ./bam/${i}_map.txt & 
done
```

## 4.生成raw bam (optional) 
```bash
vim atac2_sam2bamop.sh

#!/bin/bash
## sam to bam (samtools) ##
## sorted by position (samtools) ##

cat filenames | while read i; 
do
nohup samtools view -@ 4 -h ./bam/${i}.sam | samtools sort -@ 4 -O bam -o ./bam/${i}-sorted-pos.bam &&
samtools index -@ 4 ./bam/${i}-sorted-pos.bam & 
done

# mtReads=$(samtools idxstats ./bam/${i}-sorted-pos.bam | grep 'chrM' | cut -f 3)
# totalReads=$(samtools idxstats ./bam/${i}-sorted-pos.bam | awk '{SUM += $3} END {print SUM}')
# echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
# samtools flagstat -@ 10 ./bam/${i}-sorted-pos.bam > ./bam/${i}-sam.stat &
```

## 5.sam to bam 同时去除 ChrM   
```bash
vim atac2_sam2lastbam.sh

#!/bin/bash
## sam to bam (samtools) ##
## remove ChrM & sorted by position (samtools) ##
## remove duplication (sambamba) ##

cat filenames | while read i; 
do
nohup samtools view -@ 4 -h ./bam/${i}.sam | grep -v chrM | samtools sort -@ 4 -O bam -o ./bam/${i}-rmChrM-sorted-pos.bam && 
sambamba markdup -r -t 40 --overflow-list-size 600000 ./bam/${i}-rmChrM-sorted-pos.bam ./bam/${i}.last.bam & 
done

# samtools flagstat -@ 10 ./bam/${i}-rmChrM-sorted-pos.bam > ./bam/${i}-rmChrM-sorted-pos.stat &
# samtools flagstat -@ 10 ./bam/${i}.last.bam > ./bam/${i}.last.stat &
```

## 6.macs2（现已弃用）  
```bash
vim atac3_macs2.sh

#!/bin/bash
## peak calling (macs2) ##

cat filenames | while read i; 
do
nohup macs2 callpeak -t ./bam/${i}.last.bam -g mm --nomodel --shift -75 --extsize 150  -n ./macs2/${i} -q 0.1 --keep-dup all &  
done
```
    
## 7.使用macs3进行call peak   
```bash
vim atac4_macs3.sh

#!/bin/bash
## peak calling (macs3) ##

cat filenames | while read i; 
do
nohup macs3 callpeak -f BAMPE -t ./bam/${i}.last.bam -g mm -n ./macs3/${i} -B -q 0.1 &  
done
```

## narrowPeak和bed文件格式   

参考   
- [narrowPeak.html](https://macs3-project.github.io/MACS/docs/narrowPeak.html) 
- [01_Introduction_to_peak_files.md](https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/01_Introduction_to_peak_files.md)    

- A narrowPeak (.narrowPeak) file is used by the ENCODE project to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. The narrowPeak file is a BED6+4 format, which means the first 6 columns of a standard BED file with 4 additional fields:
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/narrowPeak.png" width="800" />

- BED files require at least 3 fields indicating the genomic location of the feature, including the chromosome and the start and end coordinates. However, there are 9 additional fields that are optional, as shown in the image below.   
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/bed_file.png" width="800" />

- Each row in the narrowPeak file represents a called peak. Below is an the example of a narrowPeak file, displaying the coordinate and statistical information for a handful of called peaks.  

---

## 8.macs3 peak文件转 bw（用于igv可视化） 
参考：  
- [Build-Signal-Track](https://github.com/macs3-project/MACS/wiki/Build-Signal-Track)
- [bedGraph to bigWig](https://gist.github.com/taoliu/2469050)  
- bedGraphToBigWig: error while loading shared libraries: libssl.so.1.0.0: cannot open shared object file: No such file or directory https://github.com/macs3-project/MACS/issues/505  
- chromInfo.txt  https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/chromInfo.txt.gz

以下脚本不用下载，直接跳到执行那一步，因为前期conda环境已经安装过这些脚本了   
>首先下载两个脚本：bedGraphToBigWig 和 bedClip（已过期）  
>```bash
>wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
>wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip  
>```
>    
>赋予可执行权限（已过期）    
>```bash
>chmod +x bedGraphToBigWig
>chmod +x bedClip
>```
>    
>添加到环境变量中（打开bashrc，添加到最后一行，保存，source）（已过期）    
>```bash
>vim ~/.bashrc
>export PATH="$PATH:/home/jjyang/downloads/bedClip"
>export PATH="$PATH:/home/jjyang/downloads/bedGraphToBigWig"
>source ~/.bashrc
>```
    
执行以下脚本  

- `-m`: FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.  
- `-m`: FE` 表示计算倍增富集。其他选项可以是 `logLR` 计算对数似然比，`subtract` 用于从处理样本中减去噪声。  
- `-p`: sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.  
- `-p`: 设置伪计数。这个数值将被加到“每百万读段的堆积”值上。在生成倍增富集轨迹时，您不需要它，因为对照组的 lambda 值将总是大于 0。但在计算对数似然比时，为了避免出现 `log(0)` 的情况，我们会添加伪计数。因为我将精度设置为五位小数，所以这里我使用 `0.00001`。  

Run MACS2 bdgcmp to generate fold-enrichment and logLR track (Then you will have this bedGraph file for fold-enrichment(FE) and logLR)   
```bash
vim atac5_bdgcmp.sh
#!/bin/bash
## peak calling (macs3) ##

cat filenames | while read i; 
do
nohup macs3 bdgcmp -t ./macs3/${i}_treat_pileup.bdg -c ./macs3/${i}_control_lambda.bdg -o ./macs3/${i}_FE.bdg -m FE &&
macs3 bdgcmp -t  ./macs3/${i}_treat_pileup.bdg -c ./macs3/${i}_control_lambda.bdg -o ./macs3/${i}_logLR.bdg -m logLR -p 0.00001 &
done
```

Fix the bedGraph and convert them to bigWig files. (And you will have these bigwig files)    
```bash
vim atac6.sh
#!/bin/bash

# check commands: slopBed, bedGraphToBigWig and bedClip

# which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
# which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
# which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
       echo "Need 2 parameters! <bedgraph> <chrom info>"
     exit
     fi
     
     F=$1
     G=$2
     
     bedtools slop -i ${F} -g ${G} -b 0 | /home/jjyang/downloads/bedClip stdin ${G} ${F}.clip
     
     LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip
     
     /home/jjyang/downloads/bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}
     
     rm -f ${F}.clip ${F}.sort.clip
```

```bash 
vim atac7.sh
#!/bin/bash
cat filenames | while read i; 
do
nohup bash atac6.sh ./macs3/${i}_FE.bdg /home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/mm10.chrom.sizes &
done  
```

最后运行以下脚本即可   
```bash 
bash atac7.sh
```

## Remove blacklist  
参考  
- [07_handling_peaks_bedtools.md](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/07_handling_peaks_bedtools.md)

> **How were the 'blacklists compiled?** These blacklists were empirically derived from large compendia of data using a combination of automated heuristics and manual curation. Blacklists were generated for various species and genome versions including human, mouse, worm and fly. The lists can be [downloaded here](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/). For human, they used 80 open chromatin tracks (DNase and FAIRE datasets) and 12 ChIP-seq input/control tracks spanning ~60 cell lines in total. These blacklists are applicable to functional genomic data based on short-read sequencing (20-100bp reads). These are not directly applicable to RNA-seq or any other transcriptome data types.
> 
> More information about the blacklist region is described in this [paper](https://www.nature.com/articles/s41598-019-45839-z). This is a more recent resource and the authors compiled blacklists that can be [downloaded here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). _This is the source for the bed file used in this workshop._ 

The black lists were downloaded from https://www.encodeproject.org/annotations/ENCSR636HFF/   

```bash 
vim atac4_rmblackls.sh

#!/bin/bash
## remove blacklist (bedtools) ##

Blacklist="/home/jjyang/downloads/genome/mm39_GRCm39/ENCFF547MET.bed"

cat filenames | while read i; 
do
nohup bedtools intersect \
-v \
-a ./macs3/${i}_peaks.narrowPeak \
-b ${Blacklist} | awk '{if($0~"chr") print}' \
> ./macs3/narrow/${i}_rmBL.narrowPeak &
done
```

## fastqc质控  
```bash
nohup fastqc -q -t 30 raw/*.fq.gz -o fqc/ &
nohup fastqc -q -t 30 trim/*.fq.gz -o trim_fqc/ &

nohup multiqc fqc/*.zip -o mqc/ &
nohup multiqc trim_fqc/*.zip -o trim_mqc/ &
```
## louvain 聚类  
```bash
conda activate atac
python ./louvain.py net_drg_select.csv net_drg_select.gexf 1 &
```
