mkdir combined -p
while read FILE
do
	${SAMTOOLS_PATH}/samtools merge -f combined/$FILE.bam oneside/$FILE.bam unmapped/$FILE.bam 
done<bamlist.txt
