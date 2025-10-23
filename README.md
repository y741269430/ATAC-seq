# ATAC-seq workflow

本文主要参考：  
- [ATAC-seq-data-analysis-from-FASTQ-to-peaks](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/)
- [Intro-to-ChIPseq-flipped](https://github.com/hbctraining/Intro-to-ChIPseq-flipped)
- [Analytical Approaches for ATAC-seq Data Analysis](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cphg.101)
- [From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)
- https://nf-co.re/chipseq/2.0.0/
- [ATAC-seq_QC_analysis](https://github.com/Zhang-lab/ATAC-seq_QC_analysis)
- [Harvardinformatics ATAC-seq](https://github.com/harvardinformatics/ATAC-seq)

基本原理：    
- 1.当发生DNA复制、基因转录时，DNA的致密高级结构变为松散状态，这部分打开的染色质被称之为开放染色质（open chromatin）。开放染色质可以和一些调控蛋白如转录因子和辅因子结合，染色质的这一特性就叫做染色质可接近性（chromatin accessibility）。染色质可接近性反映了染色质转录活跃程度。        
- 2.真核生物染色质基本结构单位是核小体，由 147bp 的DNA缠绕在组蛋白八聚体上。核小体与核小体之间由 20-90bp 的DNA linkers所连接。因此，Tn5酶转座过程，可发生在核小体之间 20-90bp 的任何地方，从连接区中产生 < 90bp 的短片段，或核小体被夹在Tn5切割的DNA中间产生更大的片段。     
- 3.文库核酸片段包含原始DNA插入片段和来自adapter两端的 135bp 测序的序列。这就形成了从 200bp - 1000bp 左右的文库片段。主要的片段堆积在160bp-200bp之间的形成高耸的peak。
- 4.如果加入的Tn5酶不足以从细胞中获得染色质，或者是Tn5与DNA的比例太低，使转座反应不充分，无法获得理想范围的DNA片段，此时片段大小估计为800bp左右，该现象被称为under-tagmentation。
- 5.解决办法：裂解均匀；延长在冰上的裂解时间；增加转座酶反应时间；在裂解前加入额外的PBS冲洗步骤，防止其他物质抑制裂解。        

生信分析流程：   
- 1.当我们数据下机之后，得到的fastq文件。使用`FastQC`软件对raw data进行质量评估。后续clean data同样需要评估。   
- 2.使用`Trimmomatic`软件对原始数据进行质控（这一步主要是去除3’端的接头污染、去除低质量序列（保留MPAQ >= 30））。得到clean data。
- 3.将clean data使用`bowtie2`软件与基因组进行比对，得到的sam文件使用`samtools`转换成bam。
- 4.得到的bam文件，获取其唯一比对以及去重复reads的结果bam文件。
- 5.使用`Deeptools`绘制TSS, Peak center 或GeneBody富集热图（依组学而定），展示数据在这些区域及前后3kb上的富集情况。
- 6.使用`MACS2`或`MACS3`进行peak calling。
- 7.使用`IDR`软件进行样品间高可信度的峰筛选.
- 8.将bam文件转换成bigwig文件，使用`IGV`进行可视化。
- 9.使用r包`ChIPseeker`对peak进行注释。
- 10.使用`homer`或`MEME`进行motif预测。
- 11.使用`MAnorm`（无生物学重复）或`DiffBind`（有生物学重复）进行差异peak分析.
- 12.假如测序深度高的情况下，可以考虑做足迹分析，即评估一个peak上是否有转录因子保护的足迹，如果有，这个peak会出现凹陷状。找出这类型的peak，然后再做motif预测，就可以找到具体是哪些转录因子在起调控作用了。比单一只做motif预测更稳健。（即：motif分析提供的是潜在的，统计相关的结合可能性；而足迹分析提供了更直接的，物理性的结合证据。）

结果指标：
- 1.文库的Q20（90%以上），Q30（80%以上）
- 2.peak注释的基因组占比（promoter占比大于20%，特殊的样本可大于10%）
- 3.FRiP值（大于0.2）
- 4.比对率（模式生物需大于85%）
- 5.片段大小分布在200，400，600，出现明显峰型
- 6.deeptools热图展示（peak主要富集在TSS区域）
- 7.call peak数无固定范围

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
mamba install trimmomatic
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
```bash
conda activate atac
```
创建文件夹  
```bash
mkdir -p trim bam macs3 logs bed
```
    
生成一个filenames的文件，用来记录输出的文件名称（样本名称），例如：  
```bash
ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 > filenames
```

## ~~2.利用Trim adaptors去除接头~~ 
```bash
vim pre_trim.sh
```
```bash
#!/bin/bash
## trim_galore ##

cat filenames | while read i; 
do
# paired end
nohup trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 1 -j 40 --paired ./raw/${i}*_1.fq.gz ./raw/${i}*_2.fq.gz -o ./trim &
  
# single end
# nohup trim_galore -q 25 --phred33 --length 20 -e 0. 1 --stringency 1 -j 40 ./raw/${i}*_1.fq.gz -o ./trim &
done
```

## 2.1 利用trimmomatic去除接头(Illumina)  
```bash
vim pre_trim.sh
```
```bash
#!/bin/bash
## trimmomatic ##

cat filenames | while read i; 
do
nohup trimmomatic PE -phred33 -threads 4 \
./RawData/${i}*_R1.fastq.gz \
./RawData/${i}*_R2.fastq.gz \
trim/${i}_forward_paired.fq.gz \
trim/${i}_forward_unpaired.fq.gz \
trim/${i}_reverse_paired.fq.gz \
trim/${i}_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &

done
```

## 3.比对到mm39 
```bash
vim atac1_bw2.sh
```
```bash
#!/bin/bash

mm39="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

cat filenames | while read i; 
do
nohup bowtie2 -p 4 --very-sensitive -X 2000 -k 10 \
-x ${mm39} \
-1 trim/${i}_forward_paired.fq.gz \
-2 trim/${i}_reverse_paired.fq.gz \
-S ./bam/${i}.sam 2> ./logs/${i}_map.txt & 
done
```
`-X` 设定最长的插入片段长度. Default: 500     
`--very-sensitive`  提高识别的敏感度但是会降低比对速度     
`-k` 默认设置下, bowtie2搜索出了一个read不同的比对结果, 并报告其中最好的比对结果(如果好几个最好的比对结果得分一致, 则随机挑选出其中一个). 而在该模式下, bowtie2最多搜索出一个read k个比对结果, 并将这些结果按得分降序报告出来。    

方法二：
```bash
#!/bin/bash

mm39="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

cat filenames | while read i; 
do
nohup bowtie2 -p 4 --very-sensitive -X 1000 \
-x ${mm39} \
-1 trim/${i}_forward_paired.fq.gz \
-2 trim/${i}_reverse_paired.fq.gz \
-S ./bam/${i}.sam 2> ./logs/${i}_map.txt & 
done
```

## hisat2/bowtie2 整合输出比对率         
输出比对率文件
```bash
conda activate rnaseq
python MappingRateOutput.py logs/ alignment_res
```

## ~~4.生成raw bam (optional)~~ 
```bash
vim atac2_sam2bamop.sh
```
```bash
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
# samtools flagstat -@ 10 ./bam/${i}-sorted-pos.bam > ./logs/${i}-sam.stat &
```

## 5.sam to bam 同时去除 ChrM   
```bash
vim atac2_sam2lastbam.sh
```
```bash
#!/bin/bash
## sam to bam (samtools) ##
## remove ChrM & sorted by position (samtools) ##
## remove duplication (sambamba) ##

cat filenames | while read i; 
do
nohup samtools view -@ 4 -h ./bam/${i}.sam | grep -v chrM | samtools sort -@ 4 -O bam -o ./bam/${i}-rmChrM-sorted-pos.bam && 
sambamba markdup -r -t 40 --overflow-list-size 600000 ./bam/${i}-rmChrM-sorted-pos.bam ./bam/${i}.last.bam & 
done

# samtools flagstat -@ 10 ./bam/${i}-rmChrM-sorted-pos.bam > ./logs/${i}-rmChrM-sorted-pos.stat &
# samtools flagstat -@ 10 ./bam/${i}.last.bam > ./logs/${i}.last.stat &
```

## ~~6.macs2（现已弃用）~~  
```bash
vim atac3_macs2.sh
```
```bash
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
```
```bash
#!/bin/bash
## peak calling (macs3) ##

cat filenames | while read i; 
do
nohup macs3 callpeak -f BAMPE -t ./bam/${i}.last.bam -g mm -n ./macs3/${i} -B -q 0.1 &  
done
```
对于ATAC-seq来说，-f BAMPE 参数仅分析“正确配对”的比对，即下图中的A。    
假如将bam转成bed，再去做call peak，这里通常是数据质量比较差的时候所作出的选择。理由是它不仅保留了双端“正确配对”的比对，也保留了单端的比对，即下图中的A和B。    
(Note that the BEDTools program bamtobed cannot be used here, since its output is in a nonstandard BED format that MACS2 cannot analyze.)    
[Harvardinformatics ATAC-seq](https://github.com/harvardinformatics/ATAC-seq/blob/master/README.md#alignments-to-analyze)    
<img src="https://github.com/y741269430/ATAC-seq/blob/main/img/alignments.png" width="500" />    


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

## 8.Remove blacklist  
参考  
- [07_handling_peaks_bedtools.md](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/07_handling_peaks_bedtools.md)

> **How were the 'blacklists compiled?** These blacklists were empirically derived from large compendia of data using a combination of automated heuristics and manual curation. Blacklists were generated for various species and genome versions including human, mouse, worm and fly. The lists can be [downloaded here](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/). For human, they used 80 open chromatin tracks (DNase and FAIRE datasets) and 12 ChIP-seq input/control tracks spanning ~60 cell lines in total. These blacklists are applicable to functional genomic data based on short-read sequencing (20-100bp reads). These are not directly applicable to RNA-seq or any other transcriptome data types.
> 
> More information about the blacklist region is described in this [paper](https://www.nature.com/articles/s41598-019-45839-z). This is a more recent resource and the authors compiled blacklists that can be [downloaded here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). _This is the source for the bed file used in this workshop._ 

The black lists were downloaded from https://www.encodeproject.org/annotations/ENCSR636HFF/   

```bash 
vim blackls_rm.sh

#!/bin/bash
## remove blacklist (bedtools) ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

Blacklist="/home/jjyang/downloads/genome/mm39_GRCm39/ENCFF547MET.bed"

cat filenames | while read i; 
do
input_file="$input_dir/${i}_peaks.narrowPeak"
output_file="$output_dir/${i}_rmBL.narrowPeak"

nohup bedtools intersect \
-v \
-a "$input_file" \
-b ${Blacklist} | awk '{if($0~"chr") print}' \
> "$output_file" &
done
```

## 9. bam 转bw（用于igv可视化）+ 绘制TSS 富集热图        
[绘制deeptools的热图](https://github.com/y741269430/ATAC-seq/blob/main/%E7%BB%98%E5%88%B6deeptools%E7%9A%84%E7%83%AD%E5%9B%BE.md)

## 10. homer 的motif预测      
[homer](https://github.com/y741269430/homer)    

## 11. macs3 peak文件转 bw（用于igv可视化）    
参考：  
- [Build-Signal-Track](https://github.com/macs3-project/MACS/wiki/Build-Signal-Track)
- [bedGraph to bigWig](https://gist.github.com/taoliu/2469050)  
- bedGraphToBigWig: error while loading shared libraries: libssl.so.1.0.0: cannot open shared object file: No such file or directory https://github.com/macs3-project/MACS/issues/505  
- chromInfo.txt  https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/chromInfo.txt.gz

以下脚本不用下载，直接跳到执行那一步，因为前期conda环境已经安装过这些脚本了   
>~~首先下载两个脚本：bedGraphToBigWig 和 bedClip~~  
>```bash
>wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
>wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip  
>```
>    
>~~赋予可执行权限（已过期）~~    
>```bash
>chmod +x bedGraphToBigWig
>chmod +x bedClip
>```
>    
>~~添加到环境变量中（打开bashrc，添加到最后一行，保存，source）~~    
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

如果嫌麻烦也可以参考：[计算重复样本的 peak 之间的重叠的坐标位置.md](https://github.com/y741269430/ATAC-seq/blob/main/%E8%AE%A1%E7%AE%97%E9%87%8D%E5%A4%8D%E6%A0%B7%E6%9C%AC%E7%9A%84%20peak%20%E4%B9%8B%E9%97%B4%E7%9A%84%E9%87%8D%E5%8F%A0%E7%9A%84%E5%9D%90%E6%A0%87%E4%BD%8D%E7%BD%AE.md#5%E4%BD%BF%E7%94%A8deeptools-%E5%B0%86bam%E8%BD%AC%E4%B8%BAbw%E7%94%A8%E4%BA%8Eigv%E5%8F%AF%E8%A7%86%E5%8C%96)

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
