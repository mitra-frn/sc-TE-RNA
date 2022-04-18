# sc-TE-RNA
*name of article*
by *name of authors*




## Description of the R scripts used in the pipeline
# Part I: Align reads to referance genome
|**bash Script**|[aligning.sh](https://github.com/am) |
| ------------- |--------------|
| **Input**|  [human/hg38 genome sequence](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)  |  
| |  [annotation of referance genome on GTF format](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) |
| | cDNAfragmentSequence.fastq and CellBarcodeUMIsequence.fastq for each sample [cDNAfragment example ](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_2.fastq.gz)and [ CellBarcodeUMI example](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/054/SRR12570154/SRR12570154_1.fastq.gz)|
| v3 cell barcode whitelist file [3M-february-2018.txt](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz) ,[repeatmasker.bed](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1333082877_7QhGx7WKpxCENEJGnGjP7lvsrSxl&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED)|
|**Output**| count matrixes on h5ad format |  
| **Dependencies**| STAR 2.7.10a,scTE(needs python >=3.6),bedtools|
|**Summary**| Genome indexes were generated and Patients' paired-end FASTQ files were mapped to them.transposable element were quantified using scTE package algorithm. |
