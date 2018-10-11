mkdir sortbyname -p
xargs -a bamlist.txt -n1 -P 2 -I FILE ${SAMTOOLS_PATH}/samtools sort -n combined/FILE.bam sortbyname/FILE.bam
