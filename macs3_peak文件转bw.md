# macs3 peak文件转 bw（用于igv可视化）    
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
