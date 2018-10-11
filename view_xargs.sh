echo Processing $1..
${SAMTOOLS_PATH}/samtools view sortbyname/$1.bam.bam > sortbyname/$1.sam
