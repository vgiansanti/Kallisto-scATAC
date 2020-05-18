# Kallisto-scATAC

This folder contains the commands and scripts to reproduce the analysis of scATAC-seq data with Kallisto-bustools that are described in the manuscript *Fast analysis of scATAC-seq data using a predefined set of genomic regions* (in preparation).

## Data

The PBMCs data used in the analysis are the one provided by the [10x Genomics public datasets](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k).

The files needed for the analysis are:

* FASTQs files
* Peak by cell matrix (filtered), in HDF5 format
* Peaks (BED)
* Clustering analysis (to retrieve the whitelist file)

In addition, a fasta file for hg19 genome is also needed, it can be retrieved from [UCSC genome browser](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz). DNase I Hypersensitive Sites can be retrieved from the appropriate [ENCODE file](http://big.databio.org/papers/RED/supplement/dhs112_v3.bed.gz). As an alternative, DNase I Hypersensitive Sites can also be retreieved from UCSC genome browser: use the Table Browser tool to retrieve, in BED format, the table `wgEncodeRegDnaseClusteredV3` (under Mammal -> Human -> hg19 -> Regulation -> DNase Clusters). 

We provide, in this repository, the bed files and the name maps used in the paper, both for DHS approach and the DHS500 approach, under the folder [reference](https://github.com/vgiansanti/Kallisto-scATAC/tree/master/reference).

The [K562_data](https://github.com/vgiansanti/Kallisto-scATAC/tree/master/K562_data) folder contains intead the count matrices in h5ad format of the K562 cell line analyzed. There is also the whitelist used to convert cellular barcodes to SRR accession numbers. 

## Building a kallisto index for ATAC analysis

In order to build an index, one needs to extract DNA sequence for a peak list. This can be done using [bedtools](https://github.com/arq5x/bedtools2) or any other suitable tool. Once bed files and genome fasta files are in place, just issue the following command, assuming that `peaks` is the prefix for bed files:

```bash
$ bedtools getfasta -fi hg19.fa -fo ${peaks}.fa -bed ${peaks}.bed
$ kallisto index -i ${peaks}.idx ${peaks}.fa
```

In order to produce the `DHS500` index, we just merged the DHS data as follows:

```bash
$ bedtools sort -i ${DHS_sites}.bed | bedtools merge -d 500 > DHS500.bed
```

and indexed it as described above. 
Kallisto analysis requires the list of entries in the index and a file mapping entries to summarized entries. The former can be extracted from fasta headers like this:

```bash
$ grep '>' ${peaks}.fa | tr -d '>' > ${peaks}.names.txt
```

or from bed files like this:

```bash
$ awk '{print $1":"$2"-"$3}' ${peaks}.bed > ${peaks}.names.txt
```

In order to create a map file we generally collate twice the `${peaks}.names.txt`:

```bash
$ awk '{print $1"\t"$1}' ${peaks}.names.txt > ${peaks}.map.txt
```

You can find the list of all maps used in the paper in this repository.

## Running kallisto for scATAC-seq 

We described in the paper several options to run kallisto with scATAC-seq data. If `n` is the number of nucleotides used for simulate the UMI, the command to pseudoalign in the *forward* configuration would be:

```bash
$ kallisto bus -t 8 -i ${peaks}.idx -o ${peaks}_${n}_fwd -x 1,0,16:2,0,${n}:0,0,0 \
pbmc_10k_R1.fastq.gz \
pbmc_10k_R2.fastq.gz \
pbmc_10k_R3.fastq.gz
```

and the command for the *reverse* configuration (with better results) is:

```bash
$ kallisto bus -t 8 -i ${peaks}.idx -o ${peaks}_${n}_rev -x 1,0,16:2,0,${n}:0,0,0 \
pbmc_10k_R3.fastq.gz \
pbmc_10k_R2.fastq.gz \
pbmc_10k_R1.fastq.gz
```

Once kallisto has finished, the count matrix can be generated like this:

```bash
$ cd ${kallisto_output}
$ bustools correct -w whitelist.txt -p output.bus | \
bustools sort -T tmp -t 4 -p - | \
bustools count -o counts/${prefix} -g ${peaks}.map.txt \
-e matrix.ec -t ${peaks}.names.txt --genecounts  -
```

where `$prefix` is the prefix for the resulting files (`${prefix}.genes.txt`, `${prefix}.barcodes.txt` and `${prefix}.mtx`). Remember to adjust the paths of `whitelist.txt`, `${peaks}.map.txt` and `${peaks}.names.txt` according to your system.

Their output files can be now analyzed with your preferred tool. Here you can find two example analysis:

* Kallisto_scATAC.ipynb: this notebook contains the analysis performed on PBMCs data with Scanpy (python commands)

* 01_kallisto_scATAC_clustering.Rmd: this contains the analysis performed on PBMCs data with Seurat (R commands)

The other files available are:

* Analysis_of_K562_cells.ipynb: the notebook with the analysis performed on the K562 cell line. The count matrices and the whitelist are the ones in the [K562_data](https://github.com/vgiansanti/Kallisto-scATAC/tree/master/K562_data) folder.

* Figure1C_MI_evaluation.ipynb: the observation of MI variation on PBMCs data can be reproduced by running this notebook.

* Figure2_Circular_staked_barplot.R: is the script to make the bar plot of the paper. 


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a [Creative Commons Attribution 4.0 International
License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
