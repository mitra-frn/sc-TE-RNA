for f in *_2.fastq.gz;do newdata/star/STAR-2.7.8a/source/STAR --soloType Droplet 
--soloCBwhitelist whitelistv3version.txt --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 
--genomeDir newdata/star/calu3/genom2/ --readFilesIn $f ${f/_2.fastq.gz/_1.fastq.gz} 
--readFilesCommand zcat --runThreadN 15 
--outSAMtyp BAM SortedByCoordinate --outFileNamePrefix ${f/_2.fastq.gz/''} 
--winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 
--outMultimapperOrder Random --outSAMmultNmax 1 --outSAMattributes NH HI nM AS CR UR CB UB;done