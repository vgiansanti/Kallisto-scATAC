# Kallisto-scATAC

This folder contains all the command and scripts to reprocude the analysis of scATAC-seq data with Kallisto-bustools that are described in the manuscript *Fast analysis of scATAC-seq data using a predefined set of genomic regions*.

## Data

The data used in the analysis are the one provided by the 10x Genomics public datasets (https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k)

The files are:

* Kallisto_scATAC.ipynb: this notebook contains the analysis performed with Scanpy (python commands)

* 01_kallisto_scATAC_clustering.Rmd: this contains the analysis performed with Seurat (R commands)

## Kallisto-bustools command for scATAC-seq

##### Generate a Kallisto index

```
kallisto index –i Dnase.idx –k 31 Dnase.fa --make-unique  
```

##### Run Kallisto
```
kallisto bus –t 8 –i Dnase.idx –o bus_output/ -x 0,0,16:0,16,26:1,0,0 –t 4 
atac_v1_pbmc_10k_S1_L001_R2_001_mod.fastq atac_v1_pbmc_10k_S1_L001_R1_001.fastq 
```

#### Build the count matrix
```
bustools correct –w ../whitelist_correct.txt –p output.bus | bustools sort –T tmp/ -t 4 –p -| 
bustools count –o genecount/genes –g ../transcript_to_genes.txt –e matrix.ec 
–t transcripts.txt --genecounts -
```
