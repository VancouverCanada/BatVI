#!/bin/bash

Threads=1
bflag=no
cargument=none

# options may be followed by one colon to indicate they have a required argument
if ! options=`getopt -u -o alt: -l along,filterdup,log:,threads: -- "$@"`
then
	exit 1
fi

eval set -- "$options"

set -- $options

Log_File="batvirun.log"
Dup_Filter="";
while [ $# -gt 0 ]
do
	case $1 in
		-a|--along) aflag="yes" ;;
		-l|--log) 
			Log_File=$2 ; shift;;
		-t|--threads) 
			Threads=$2 ; shift;;
		-f|--filterdup) 
			Dup_Filter=1 ; shift;;
		(--) shift; break;;
		(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
		(*) break;;
	esac
	shift
done
#if [ "$2" == "" ];then echo "Please enter the number of threads..";exit 1;fi
