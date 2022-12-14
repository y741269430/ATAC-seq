# ATAC-seq

## 0. Build source used for ATAC-seq  

    conda create -n atac
    conda activate atac

    conda install -c bioconda trim-galore 
    conda install -c bioconda bowtie2 
    conda install -c bioconda macs2
    conda install -c bioconda samtools
    conda install -c bioconda sambamba
    conda install -c bioconda bedtools
    conda install -c bioconda picard
    
    conda create -n macs3 python=3.8
    conda activate macs3
    conda install -c maximinio macs3
    
    conda create -n mqc python=3.6
    conda activate mqc
    pip install multiqc
    
## 0. Build the bowtie2 reference genome index (mm39)

    cd /home/yangjiajun/downloads/genome/mm39_GRCm39/ucsc_fa/
    
    nohup bowtie2-build GRCm38.primary_assembly.genome.fa \
    /home/yangjiajun/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39 & 

## 1. Activate the source and create the folder  
    
    conda activate atac  
    
    mkdir -p raw clean trim bam macs2 macs2/narrow macs3 macs3/narrow  
    
## 2. Write the filenames  

    ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 > filenames

## 3. Trim adaptors

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

## 4. Alignment to mm39  

    vim atac1_bw2.sh

    #!/bin/bash
    ## Alignment to mm39 ##

    mm39="/home/yangjiajun/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

    cat filenames | while read i; 
    do
    nohup bowtie2 -p 4 --very-sensitive -X 2000 -k 10 \
        -x ${mm39} \
        -1 trim/${i}*_val_1.fq.gz \
        -2 trim/${i}*_val_2.fq.gz \
        -S ./bam/${i}.sam 2> ./bam/${i}_map.txt & 
    done

## Generate raw bam (optional) 

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


## 5. sam to bam and remove ChrM    

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

## 6. macs2 

    vim atac3_macs2.sh

    #!/bin/bash
    ## peak calling (macs2) ##

    cat filenames | while read i; 
    do
    nohup macs2 callpeak -t ./bam/${i}.last.bam -g mm --nomodel --shift -75 --extsize 150  -n ./macs2/${i} -q 0.1 --keep-dup all &  
    done
    
## 7. macs3 

    vim atac4_macs3.sh

    #!/bin/bash
    ## peak calling (macs3) ##

    cat filenames | while read i; 
    do
    nohup macs3 callpeak -f BAMPE -t ./bam/${i}.last.bam -g mm -n ./macs3/${i} -B -q 0.1 &  
    done
    
## Remove blacklist  

The black lists was downloaded from https://www.encodeproject.org/annotations/ENCSR636HFF/  

    vim atac4_rmblackls.sh

    #!/bin/bash
    ## remove blacklist (bedtools) ##

    cat filenames | while read i; 
    do
    nohup bedtools intersect -v -a ./macs2/${i}_peaks.narrowPeak -b /home/yangjiajun/downloads/mm10.blacklist_ENCFF547MET.bed | awk '{if($0~"chr") print}' > ./macs2/narrow/${i}_rmBL.narrowPeak & 
    done

## fastqc  

    nohup fastqc -q -t 30 raw/*.fq.gz -o fqc/ &
    nohup fastqc -q -t 30 trim/*.fq.gz -o trim_fqc/ &

    nohup multiqc fqc/*.zip -o mqc/ &
    nohup multiqc trim_fqc/*.zip -o trim_mqc/ &

## IDR计算overlap  

    conda create -n idr
    conda activate idr
    conda install -c bioconda idr

    nohup sort -k8,8nr p1.narrowPeak > b1 && sort -k8,8nr p2.narrowPeak > b2 && idr --samples b1 b2 --input-file-type narrowPeak --rank p.value --output-file b12 --plot --log-output-file b12.log &

## intervene 计算peaks之间的overlaping  

    conda create -n intervene
    conda activate intervene
    conda install -c bioconda intervene

    nohup intervene venn  -i ../macs2/narrow/*.narrowPeak --save-overlaps &
    nohup intervene upset -i ../m1/*_rmBL.narrowPeak --output ./ &

## deeptools 计算peaks之间的overlaping和correlation  

    conda create -n deeptools
    conda activate deeptools
    conda install -c bioconda deeptools
    
    nohup multiBamSummary bins --bamfiles bam/*last.bam --minMappingQuality 30 --labels BL6-TG15-ATAC-CT BL6-TG16-ATAC-CT BL6-TG-ATAC-C2 BL6-TG-ATAC-C4 BL6-TG-ATAC-C5 BL6-TG-ATAC-C6 -out readCounts.npz --outRawCounts readCounts.tab && 

    nohup plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --log1p --removeOutliers -p scatterplot -o scatterplot_SpM.pdf --outFileCorMatrix Spearman.tab &

    vim g3_bam2bw.sh

    #!/bin/bash
    ## bam to bw ##

    cat filenames | while read i; 
    do
    nohup bamCoverage --bam ./bam/${i}.last.bam -o ./bw/${i}.bw --binSize 10 --normalizeUsing RPKM & 
    # nohup bamCoverage --bam ./bam/${i}.last.bam -o ./bw2/${i}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads & 
    done

### 可视化  

    nohup multiBigwigSummary bins -b *.bw -o test.npz && plotCorrelation -in test.npz --corMethod spearman --skipZeros --log1p --removeOutliers -p scatterplot -o scatterplot_SpM.pdf --outFileCorMatrix Spearman.tab &

## louvain 聚类  

    conda activate atac

    python ./louvain.py net_drg_select.csv net_drg_select.gexf 1 &
