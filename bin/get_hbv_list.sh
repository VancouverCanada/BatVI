set -e
set -u
set -o pipefail

RUN=$1
mkdir ${RUN}/blast ||true 
cat head.hbv.txt tail.hbv.txt > ${RUN}/blast/unmapped.hbv.txt || return 1
awk '/^@/' ${RUN}/blast/unmapped.hbv.txt >${RUN}/tmp/hbv.list.txt || return 1
sed -i 's/\/2//' ${RUN}/tmp/hbv.list.txt || return 1
sed -i 's/\/1//' ${RUN}/tmp/hbv.list.txt || return 1
sort ${RUN}/tmp/hbv.list.txt -u |sed 's/^@//' > ${RUN}/tmp/hbv.final.list || return 1
