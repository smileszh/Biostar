# (PART) fastq 质量控制 {-}



# fastp
官方教程：https://github.com/OpenGene/fastp

使用示范:

单端数据

``` bash
SRR=SRR12351448
DATA=data
READS=${DATA}/reads

make -f src/run/fastp.mk P1=${READS}/${SRR}_1.fastq run
```


``` bash
DATA=data
READS=${DATA}/reads
design=${DATA}/design.csv

cat ${design} | parallel --header : --colsep , \
        make -f src/run/fastp.mk \
        P1=${READS}/{Run}_1.fastq \
        run
```


双端数据

``` bash
SRR=SRR12351448
DATA=data
READS=${DATA}/reads

make -f src/run/fastp.mk P1=${READS}/${SRR}_1.fastq P2=${READS}/${SRR}_2.fastq run
```


``` bash
DATA=data
READS=${DATA}/reads
design=${DATA}/design.csv

cat ${design} | parallel --header : --colsep , \
        make -f src/run/fastp.mk \
        P1=${READS}/{Run}_1.fastq \
        P2=${READS}/{Run}_2.fastq \
        run
```

[sra.mk](src/run/fastp.mk) 

