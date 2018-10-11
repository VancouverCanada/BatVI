mkdir oneside -p
mkdir unmapped -p

xargs -a bamlist.txt -n1 -P 2 -I FILE ${SAMTOOLS_PATH}/samtools view -F4 -f 8 -bo oneside/FILE.bam sam/FILE.bam  
xargs -a bamlist.txt -n1 -P 2 -I FILE ${SAMTOOLS_PATH}/samtools view -f 4 -bo unmapped/FILE.bam sam/FILE.bam
