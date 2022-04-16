# sc-TE-RNA
*name of article*
by *name of authors*




## Description of the R scripts used in the pipeline
# Part I: Align reads to referance genome
### 1. generate human indexs
|**R Script**|[aligning.R](https://github.com/am) |
| ------------- |--------------|
| **Input**|  [human/hg38 genome sequence](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)  |  
| |  [annotation on GTF format](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) |
|**Output**| A |  
| **Dependencies**| data.table, FactoMineR, factoextra|
|**Summary**|Preparation of a binary data matrix for molecular initiating events (MIEs) from compound-gene interactions in CTD. The interactions subtypes related to metabolism were grouped as metabolism and the interaction types related to transport are grouped as transport. Performing multiple correspondence analysis on the resulting matrix uisng FactoMineR and factoextra. Selection of reaction,binding,activity,expression,metabolic processing as the more distant types of the interaction based on the plot of MCA. For the compounds with more than 50 gene interactions the less informative gene interactions will be removed.|
