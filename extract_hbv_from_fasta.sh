set -e
set -u
set -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ ! -e $1 ]
then
	echo "Cannot find directory $1..";exit 1;
fi
cd $1
mkdir log -p
if [ $? -ne 0 ];then echo "Cannot create directory $? ..";exit 1;fi
if [ -e "log/step1" ]
then 
	echo "Step 1 seems complete. Moving to next step.."
	cd $OLDPWD 
else
	if [ ! -e "filelist.txt" ]
	then
		echo "Please provide a list of fasta files (gzipped or plain) to process in the file filelist.txt"
		exit 100
	fi

	mkdir out -p
	if [ $? -ne 0 ]
	then 
		echo "Error Creating out directory.."
		exit 100
	fi

	while read FILE
	do
		Arr=(`echo $FILE | tr ";" " "`)
		File1=${Arr[0]}
		File2=${Arr[1]}
		Insertsize=${Arr[2]}
		if [ "$Insertsize" == "" ];then Insertsize=1000;fi

		if [ "$File2" == "" ]
		then
			time $DIR/batindel/src/penguin -g ${INDEX} -q $File1 -o out/${File1} --threads $Threads
		else
			time $DIR/batindel/src/penguin -g ${INDEX} -q $File1 -q $File2 -o out/${File1##*/} --threads $Threads -s $Insertsize
		fi
#--- Error Handle ---
    if [ $? -ne 0 ]
    then
      exit 1
    fi
#--------------------
	done<filelist.txt

	touch log/step1
	cd  $OLDPWD
fi
