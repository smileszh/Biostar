


# hisat2
官方教程：https://daehwankimlab.github.io/hisat2/

## 构建索引


``` bash
REF=~/database/Human/Homo_sapiens.GRCh38.dna.toplevel.fa

make -f ~/src/run/hisat2.mk REF=${REF} index
```

## 比对

单端数据

``` bash
DATA=data
REF=genome/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
HISAT2=results/hisat2
CLEANDATA=${DATA}/reads
Run=SRR1343245

make -f src/run/hisat2.mk \
        REF=${REF} \
        R1=${CLEANDATA}/${Run}_1.trimmed.fastq \
        BAM=${HISAT2}/{sample}.bam \
        run
```


``` bash
DATA=data
design=${DATA}/design.csv
REF=genome/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
CLEANDATA=${DATA}/reads
HISAT2=results/hisat2
                
cat ${design} | parallel --header : --colsep , \
        make -f src/run/hisat2.mk \
        REF=${REF} \
        R1=${CLEANDATA}/{Run}_1.trimmed.fastq \
        BAM=${HISAT2}/{sample}.bam \
        run
```

双端数据


``` bash
DATA=data
REF=genome/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
HISAT2=results/hisat2
CLEANDATA=${DATA}/reads
Run=SRR1343245

make -f src/run/hisat2.mk \
        REF=${REF} \
        R1=${CLEANDATA}/${Run}_1.trimmed.fastq \
        R2=${CLEANDATA}/${Run}_2.trimmed.fastq
        BAM=${HISAT2}/{sample}.bam \
        run
```



``` bash
DATA=data
design=${DATA}/design.csv
REF=genome/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
CLEANDATA=${DATA}/reads
HISAT2=results/hisat2
                
cat ${design} | parallel --header : --colsep , \
        make -f src/run/hisat2.mk \
        REF=${REF} \
        R1=${CLEANDATA}/{Run}_1.trimmed.fastq \
        R2=${CLEANDATA}/{Run}_2.trimmed.fastq \
        BAM=${HISAT2}/{sample}.bam \
        run
```

[sra.mk](src/run/hisat2.mk) 
