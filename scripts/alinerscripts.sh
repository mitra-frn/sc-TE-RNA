#build genome-index for STARaligner
#input:hg38.fa(human/hg38 genome sequence),hg38.knownGene.gtf(annotation of referance genome on GTF format)
#output:genome indexes are built and stored in GENOMEDIR
 
STAR-2.7.9a/source/STAR --runThreadN 14 --runMode genomeGenerate --genomeDir GENOMEDIR 
--genomeFastaFiles hg38.fa --sjdbGTFfile hg38.knownGene.gtf --sjdbOverhang 100
---------------------------------------------
#align paire-end fastq to genome index
#input:f(names of samples without 1/2.fastq.gz),3M-february-2018.txt(v3 cell barcode whitelist file),GENOME2(directory of genome-index)
#output:sorted alignments in binary BAM format



for f *_2.fastq.gz;do STAR-2.7.9a/source/STAR --soloType Droplet --soloCBwhitelist 3M-february-2018.txt --soloCBlen 16
 --soloUMIstart 17 --soloUMIlen 12 --genomeDir GENOME2 --readFilesIn $f ${f/_2.fastq.gz/_1.fastq.gz} --readFilesCommand zcat
 --runThreadN 15 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${f/_2.fastq.gz/''} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 
--outMultimapperOrder Random --outSAMmultNmax 1 --outSAMattributes NH HI nM AS CR UR CB UB;done
--------------------------------------------
#build genome-inedx for count transposable elements
#input:hg38.fa(human/hg38 genome sequence),hg38.knownGene.gtf(annotation of referance genome on GTF format),repeatmasker.bed(UCSC genome browser Repeatmasker track)
#output:hg38.exclusive.idx(indexes in exclusive/nointron mood)


scTE_build -g hg38.fa -te repeatmasker.bed -gene hg38.knownGene.gtf -m exclusive/nointron -o hg38.exclusive.idx/hg38.nointron.idx

----------------------------------------------
#Analysis of 10x style scRNA-seq data generated by STARSolo
#input:f(sorted alignments in binary BAM format),hg38.exclusive.idx(genome indexes in exclusive mood)
#output:(count-matrix on h5ad format)


for f in *.bam;do scTE -i $f -o ${f/.bam/''} -x hg38.exclusive.idx --hdf5 True -CB CR -UMI UR  --expect-cells 200000;done
