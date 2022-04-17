# sc-TE-RNA
*name of article*
by *name of authors*




## Description of the R scripts used in the pipeline
# Part I: Align reads to referance genome
### 1. generate human indexs
|**bash Script**|[aligning.sh](https://github.com/am) |
| ------------- |--------------|
| **Input**|  [human/hg38 genome sequence](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)  |  
| |  [annotation of referance genome on GTF format](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) |
| | cDNAfragmentSequence.fastq and CellBarcodeUMIsequence.fastq for each sample [cDNAfragment example](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_2.fastq.gz)[CellBarcodeUMI example](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_1.fastq.gz)|
|**Output**|  |  
| **Dependencies**| STAR 2.7.10a,|
|**Summary**|.|
