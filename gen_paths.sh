P_tmp=`which blastn 2>/dev/null`;P_tmp=${P_tmp%/blastn}
if [ "$P_tmp" != "" ]
then
	echo BLAST_PATH=$P_tmp
fi

P_tmp=`which bwa 2>/dev/null`;P_tmp=${P_tmp%/bwa}
if [ "$P_tmp" != "" ]
then
	echo BWA_PATH=$P_tmp
fi
P_tmp=`which samtools 2>/dev/null`;P_tmp=${P_tmp%/samtools}
if [ "$P_tmp" != "" ]
then
	echo SAMTOOLS_PATH=$P_tmp
fi

P_tmp=`which bedtools 2>/dev/null`;P_tmp=${P_tmp%/bedtools}
if [ "$P_tmp" != "" ]
then
	echo BEDTOOLS_PATH=$P_tmp
fi

