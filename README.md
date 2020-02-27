# Kallisto-scATAC

This folder contains all the command and scripts to reprocude the analysis of scATAC-seq data with Kallisto-bustools that are described in the manuscript *Fast analysis of scATAC-seq data using a predefined set of genomic regions*.

## Data

The data used in the analysis are the one provided by the 10x Genomics public datasets (https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k)

The files are:

* Kallisto_scATAC.ipynb: this notebook contains the analysis performed with Scanpy (python commands)

* 01_kallisto_scATAC_clustering.Rmd: this contains the analysis performed with Seurat (R commands)

## Kallisto-bustools command for scATAC-seq

Here are reported the commands to build a count matrix for scATAC-seq data with Kallisto.

The **_kallisto index_** can be obtained from any collection of annotated peaks in fasta format.

**_kallisto bus_** for ATAC requires to modify the fastq files. 

To completely understand how these steps are performed, please see the Mehods section of the paper.

##### Generate a Kallisto index

```
kallisto index –i Dnase.idx –k 31 Dnase.fa --make-unique  
```

##### Run Kallisto
```
kallisto bus –t 8 –i Dnase.idx –o bus_output/ -x 0,0,16:0,16,26:1,0,0 –t 4 
atac_pbmc_10k_R2.fastq atac_pbmc_10k_R1.fastq 
```

#### Build the count matrix
```
bustools correct –w ../whitelist_correct.txt –p output.bus | bustools sort –T tmp/ -t 4 –p -| 
bustools count –o genecount/genes –g ../transcript_to_genes.txt –e matrix.ec 
–t transcripts.txt --genecounts -
```
