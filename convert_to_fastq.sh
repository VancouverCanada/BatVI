mkdir $1/fastq -p
xargs -a bamlist.txt -n1 -P 2 -I FILE java -jar ${PICARD_PATH}/SamToFastq.jar INPUT=$1/FILE.sam FASTQ=$1/fastq/FILE_1.fq SECOND_END_FASTQ=$1/fastq/FILE_2.fq VALIDATION_STRINGENCY=SILENT 
