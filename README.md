# sc-TE-RNA
*name of article*

by *name of authors*




## Description of scripts used in the pipeline
# Part I: Align reads to referance genome
|**bash Script**|[alinerscripts.sh](scripts/alinerscripts.sh) |
| ------------- |--------------|
| **Input**|  [human/hg38 genome sequence](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)  |  
| |  [annotation of referance genome on GTF format](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) |
| | cDNAfragmentSequence.fastq and CellBarcodeUMIsequence.fastq for each sample [cDNAfragment example ](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_2.fastq.gz)and [ CellBarcodeUMI example](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_1.fastq.gz)|
| |v3 cell barcode whitelist file [3M-february-2018.txt](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz)|
| |[repeatmasker.bed](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1333082877_7QhGx7WKpxCENEJGnGjP7lvsrSxl&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED)|
|**Output**| count matrixes on h5ad format |  
| **Dependencies**| STAR 2.7.10a,scTE(needs python >=3.6),bedtools|
|**Summary**| Genome indexes were generated and Patients' paired-end FASTQ files were mapped to them.transposable element were quantified using scTE package algorithm. |

<br/>
<br/>
<br/>

# Part II: Preprocessing, analysis, and exploration of  scRNA-seq data
|**R Script**|[Seuratscripts.R](scripts/Seuratscripts.R) |
| ------------- |--------------|
| **Input**|  count matrixes on h5ad format[ examples of count matrix](https://github.com/mitra-frn/sc-TE-RNA/blob/main/Inputs/count_matrix_example.zip)  |  
|**Output**|  Defferential expression tables for each cluster in compares to other clusters|  
| **Dependencies**| zellkonverter,SingleCellExperiment,Seurat|
|**Summary**| count matrixes were preprocessed(ig.,filtering,normalizaton).The bach effects between samples were removed and all-samples of each condition(ie,healthy,mild and severe conditions) were integrated.After scaling data, PCA was performed and The KNN graph was conducted based on the PCA-reduced data, and unsupervised clustering was performed.Differential expression analysis was conducted based on the 'MAST' method.|

<br/>
<br/>
<br/>

# Part III: Scoring pathways and correlation test
|**R Script**|[scoringscripts.R](scripts/scoringscripts.R) |
| ------------- |--------------|
| **Input**|  seurat-object of samples[ example of seurat metadata]|
| | gene-sets of pathways[ example](genset) |  
|**Output**|  correlation tables and relating p.adjvalue|  
| **Dependencies**| Seurat,psych|
|**Summary**| Average expression levels of  gene-sets pathways in neutrophil clusters and upregulated genes in these clusters(logfc>0.7 and p.adj<0.05)  were calculated and the correlation between these pathways and upregulated genes in neutrophil clusters was evaluated.|

