# macs3 peak文件转 bw（用于igv可视化）    
参考：  
- [Build-Signal-Track](https://github.com/macs3-project/MACS/wiki/Build-Signal-Track)
- [bedGraph to bigWig](https://gist.github.com/taoliu/2469050)  
- bedGraphToBigWig: error while loading shared libraries: libssl.so.1.0.0: cannot open shared object file: No such file or directory https://github.com/macs3-project/MACS/issues/505  
- chromInfo.txt  https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/chromInfo.txt.gz
---

## 1.Run MACS2 main program        
```bash
macs2 callpeak -t ${Sample_name}_IP.bam -c ${Sample_name}_Input.bam  -B --nomodel --extsize 147 --SPMR -g ce -n ${Sample_name}
```
* `--nomodel` 和 `--extsize 147` 参数指示 MACS2 使用 147 bp 作为片段大小来堆积测序读段。        
* `-g ce` 参数让 MACS2 将线虫（C. elegans）基因组视为背景。对于果蝇请设为 `dm`，对于人类请设为 `hs`。        
* `-B --SPMR` 参数要求 MACS2 以 bedGraph 格式生成“每百万读段的片段堆积”信号文件。        
运行后，您将在同一目录下得到以下两个 bedGraph 文件：
```     
${Sample_name}_treat_pileup.bdg
${Sample_name}_control_lambda.bdg
```
## 2.Run MACS2 bdgcmp to generate fold-enrichment and logLR track        
```bash
vim atac5_bdgcmp.sh
```
```bash
#!/bin/bash
## peak calling (macs3) ##

cat filenames | while read i; 
do
nohup macs3 bdgcmp -t ./macs3/${i}_treat_pileup.bdg -c ./macs3/${i}_control_lambda.bdg -o ./macs3/${i}_FE.bdg -m FE &&
macs3 bdgcmp -t  ./macs3/${i}_treat_pileup.bdg -c ./macs3/${i}_control_lambda.bdg -o ./macs3/${i}_logLR.bdg -m logLR -p 0.00001 &
done
```
```bash
bash atac5_bdgcmp.sh
```
* `-m FE` 表示计算富集倍数。其他可选参数包括用于计算对数似然比的 `logLR` 和用于从处理样本中扣除背景噪声的 subtract。            
* `-p` 设置伪计数。这个数值将被添加到“每百万读段堆积”值中。生成富集倍数轨迹时不需要此参数，因为对照的 λ 值总是大于 0。但为了避免计算对数似然比时出现 log(0) 的情况，我们需要添加一个伪计数。由于我将精度设置为小数点后 5 位，因此这里使用了 `0.00001`。        

## 3.Fix the bedGraph and convert them to bigWig files. (And you will have these bigwig files)    
使用以下脚本将基因组索引生成chromInfo.txt      
```bash
awk '{print $1 "\t" $2 "\t" $1}' Homo_sapiens.GRCh38.dna_sm.chromosome.chr.fa.fai > ./chromInfo.txt
```
写入以下脚本    
```bash
vim bedGraphToBigWig.sh
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
创建批量运行的命令    
```bash 
vim atac_macs3tobw.sh
```
```bash 
#!/bin/bash
cat filenames | while read i; 
do
nohup bash bedGraphToBigWig.sh ./macs3/${i}_FE.bdg /home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/mm10.chrom.sizes &
done  
```
最后运行以下脚本即可       
```bash 
bash atac_macs3tobw.sh
```

这里面是使用bam文件转bigwig：[计算重复样本的 peak 之间的重叠的坐标位置.md](https://github.com/y741269430/ATAC-seq/blob/main/%E8%AE%A1%E7%AE%97%E9%87%8D%E5%A4%8D%E6%A0%B7%E6%9C%AC%E7%9A%84%20peak%20%E4%B9%8B%E9%97%B4%E7%9A%84%E9%87%8D%E5%8F%A0%E7%9A%84%E5%9D%90%E6%A0%87%E4%BD%8D%E7%BD%AE.md#5%E4%BD%BF%E7%94%A8deeptools-%E5%B0%86bam%E8%BD%AC%E4%B8%BAbw%E7%94%A8%E4%BA%8Eigv%E5%8F%AF%E8%A7%86%E5%8C%96)
