#N=0;while read F;do N=$((N+1));sed -i "s/@/@$N:170:/" unbug/fastq/${F}_1.fq;done<bamlist170.txt
#N=0;while read F;do N=$((N+1));sed -i "s/@/@$N:170:/" unbug/fastq/${F}_2.fq;done<bamlist170.txt
#N=0;while read F;do N=$((N+1));sed -i "s/@/@$N:800:/" unbug/fastq/${F}_1.fq;done<bamlist800.txt
#N=0;while read F;do N=$((N+1));sed -i "s/@/@$N:800:/" unbug/fastq/${F}_2.fq;done<bamlist800.txt


#N=0;while read F;do N=$((N+1));sed -i "s/^/$N:170:/" unbug/${F}.sam;done<bamlist170.txt
#N=0;while read F;do N=$((N+1));sed -i "s/^/$N:800:/" unbug/${F}.sam;done<bamlist800.txt
