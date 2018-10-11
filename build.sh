if [ "$1" == "clean" ]
then
	cd BatMis-3.00
	make clean
	cd ../batindel
	make clean
	cd ../bin
	rm cluster_bp
	rm *.o
	cd ../msapipeline/msa
	rm *.class
else
	echo ==========================================================
	echo                COMPILING BatMis
	echo ==========================================================
	cd BatMis-3.00
	./configure
	make
	make copy
	cd -

	echo ==========================================================
	echo                COMPILING BATINDEL-lite
	echo ==========================================================
	cd batindel
	./configure
	make
	cd -

	echo ==========================================================
	echo                COMPILING Binaries 
	echo ==========================================================
	cd bin
	./build.sh
	cd -

	echo ==========================================================
	echo                COMPILING MSA
	echo ==========================================================
	cd msapipeline/msa
	javac *.java
	cd -

  ./manualcompile.sh

	command -v blastn >/dev/null 2>&1 || { echo >&2 "Please check if BLAST is installed"; }
	command -v bwa >/dev/null 2>&1 || { echo >&2 "Please check if BWA is installed"; }
fi
