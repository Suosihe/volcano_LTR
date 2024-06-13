# volcano_LTR
The pipeline to characterize the LTR-RTs family, classify and predict the burst families.

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

You need the `LTR_retriver` results to run the `volcano` pipeline. If you run `EDTA`, you may found the input files in `*mod.EDTA.raw/LTR/`. To get more details for input-file-preparation, you can see [Wiki](https://github.com/Suosihe/volcano_LTR/wiki).

### Classification

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

Using the volcano.sh script, four input files are essential.

`-l` option is the result of `LTR_retriver` (in `EDTA`) and ends with `pass.list`.

`-f` option is the reference genome, which usually ends in mod if you are using `EDTA` as an input file.

`-n` option is the prefix you defined.

`-s` option is the genome size file. It only contains the length number of the genome.

The optional parameters：

`-c` sequence identity threshold for cd-hit default 0.8

`-L` alignment coverage for the longer sequence in cd-hit, default 0.8

`-T` threads for cd-hit default 0, use all CPUs

`-M` memory limit (in MB) for cd-hit default 0 for no limit

`-p` parallel number for RepeatMasker default 80

`-v` 'div' in RepeatMasker Masks only those repeats with (number) percent diverged from consensus, default 40

If you want to prepare the input files more conveniently, you may visit [Wiki](https://github.com/Suosihe/volcano_LTR/wiki).


### Quantification

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

問題：

***定量結果加上family名？***

Find `*pass.list.gff3` in `EDTA_raw/LTR` directory or `*pass.list.gff3` in `LTR retriever` directory.

The `*sorted.sam` file is obtained by comparing the RNA-seq file with the reference genome. If you are unsure about the file, you can refer to the [Wiki](https://github.com/Suosihe/volcano_LTR/wiki). 

`-p` option is a prefix for the results file.

`-n` option is the genome size file.

## Output

### Classification

### Quantification

The obtained `prefix_RPKM.tsv` is the quantitative file.
