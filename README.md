# volcano_LTR
The pipeline to characterize the LTR-RTs family, classify and predict the burst families

## Installation

`perl >= 5.0`​

Categorical dependence: `samtools`,`RepeatMasker`​,`blast`​,`fasttree`​

Quantitative dependence:`telescope`,`r-optparse`,`r-readr`,`r-dplyr`

`R >= 4.0`​

We recommend installing with `conda`

待補充

***Give conda yml***

manual installation

```shell
# Categorical dependence
conda install -c bioconda samtools 
conda install -c bioconda blast
conda install -c bioconda cd-hit
conda install -c bioconda fasttree
conda install -c bioconda mafft
# Quantitative dependence
conda install -c bioconda telescope
conda install -c conda-forge r-base r-optparse r-readr r-dplyr
```


```shell
git clone https://github.com/Suosihe/volcano_LTR.git
chmod a+x prepare.sh volcano.sh tel.sh
export PATH=/path/to/volcano
```

## Usage

### Preparation - EDTA-dependent operation

When you have finished `EDTA`, you can use our script to get the result file of `LTR_retriver` in `EDTA` to be used as the input file for subsequent classification and quantification.

`-d` option is the name of the folder you created.

`-E` option is the location to run `EDTA`, there are usually folders such as `EDTA.raw`, `EDTA.final` under the location.

`-f` option is the genome file used to run `EDTA`, which is usually in fasta format with a 'mod' suffix.


```shell
prepare.sh 
Usage: prepare.sh [-d dir_prefix] [-E /path/to/EDTA e.g. ~/A_thaliana/EDTA]
       [-f /path/to/fasta(for EDTA) e.g. ~/A_thaliana/EDTA/A_thaliana.fa.mod]
```

### Preparation - LTR_retriever-dependent operation

If you have run the `LTR retriever` manually, the `pass.list` file required for `volcano.sh` can be found in the output directory. 

The genome size file can be obtained by fai indexing. For example

```shell
samtools faidx xxx.fasta
tail -n 1 xxx.fasta.fai|cut -f 3 > len
```

### Classification

問題

***Question: Families names***

Using the volcano.sh script, four input files are essential.

`-l` option is the result of `LTR_retriver` (in `EDTA`) and ends with `pass.list`.

`-f` option is the reference genome, which usually ends in mod if you are using `EDTA` as an input file.

`-n` option is the prefix you defined.

`-s` option is the genome size file. If you are using `prepare.sh`, the file name is 'len'.

The optional parameters：

`-c` sequence identity threshold for cd-hit default 0.8

`-L` alignment coverage for the longer sequence in cd-hit, default 0.8

`-T` threads for cd-hit default 0, use all CPUs

`-M` memory limit (in MB) for cd-hit default 0 for no limit

`-p` parallel number for RepeatMasker default 80

`-v` 'div' in RepeatMasker Masks only those repeats with (number) percent diverged from consensus, default 40

```
volcano.sh
Usage: volcano.sh
                   Required parameters:
               [-l reference.fa.pass.list] [-f reference.fa]
               [-n prefix] [-s genome_size_bp_file]
                   Optional parameters:
               [-c sequence identity threshold for cd-hit] default 0.8
               [-L alignment coverage for the longer sequence in cd-hit], default 0.8
               [-T threads for cd-hit] default 0, use all CPUs
               [-M memory limit (in MB) for cd-hit] default 0 for no limit
               [-p parallel number for RepeatMasker] default 80
               [-v 'div' in RepeatMasker] Masks only those repeats with (number) percent diverged from consensus, default 40
```

### More precise categorisation

Import `copia.RT.tree` and `gypsy.RT.tree` in itol, then import the colour annotation file, sort and manually export by category (colour).

***待補充***

### Quantification

問題：

***1. 能否成功得到st.sam
2. 定量結果加上family名？***

Find `pass.list.gff` in `EDTA_raw/LTR` directory or `pass.list.gff` in `LTR retriever` directory.

The `sorted.sam` file is obtained by comparing the RNA-seq file with the reference genome.

The file can be obtained with the following command:


```shell
hisat2-build -p 10 xxx.fasta.mod xxxmod
hisat2 -x xxxmod -1 xxxRNAseq1_R1.fq.gz -2 xxxRNAseq1_R2.fq.gz -S xxx1.sam -p 8
hisat2 -x xxxmod -1 xxxRNAseq2_R1.fq.gz -2 xxxRNAseq2_R2.fq.gz -S xxx2.sam -p 8
hisat2 -x xxxmod -1 xxxRNAseq3_R1.fq.gz -2 xxxRNAseq3_R2.fq.gz -S xxx2.sam -p 8
samtools view -@ 4 -b xxx1.sam |samtools sort -@ 4 -O BAM -o xxx1.st.bam
samtools view -@ 4 -b xxx2.sam |samtools sort -@ 4 -O BAM -o xxx2.st.bam
samtools view -@ 4 -b xxx3.sam |samtools sort -@ 4 -O BAM -o xxx3.st.bam
samtools merge -O SAM merged.st.sam xxx1.st.bam xxx2.st.bam
```

`-p` option is a prefix for the results file.

`-n` option is the genome size file.

```
tel.sh
Usage: tel.sh     Required parameters:
               [-g pass.list.gff3]  [-s sorted.sam]  [-p prefix]  [-n number]

* [pass.list.gff3] must be defined in the script. It is in your LTR_retriver path or [EDTA.raw/LTR/]. 

* [sorted.sam] must be defined in the script. It is a sorted sam file usually contains alignment results. 

* [prefix] must be defined in the script. It is the prefix you defined for your output. 

* [number, genome_size_bp] must be defined in the script. It is a file which only contains genome size (bp) number. It can be got from index file. 

* -h, --help: print this help and exit
```

The obtained `prefix_RPKM.tsv` is the quantitative file.
