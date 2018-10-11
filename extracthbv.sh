F=$1
F=${F##*/}
echo processing $F..
${SAMTOOLS_PATH}/samtools view  $1.bam| grep HBV-B2 > $F.sam
java -jar ${PICARD_PATH}/SamToFastq.jar INPUT=$F.sam F=fastq/${F}_1.fq F2=fastq/${F}_2.fq VALIDATION_STRINGENCY=SILENT
rm $F.sam
