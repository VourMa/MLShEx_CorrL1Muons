#!/bin/bash
condorFile=ResolutionsCondor_$1_$2_$3_$4_$5.sub
cp ResolutionsTemplate.sub $condorFile

if [ ! -f "CondorLog.txt" ]; then
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" > CondorLog.txt
else
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" >> CondorLog.txt
fi

chmod 755 $condorFile
echo "$condorFile created!"

sed -i "s|DATASET|$1|" $condorFile; shift
sed -i "s|ERA|$1|" $condorFile; shift
sed -i "s|PERFORMONWHICHGUYS|$1|" $condorFile; shift

sed -i "s|QUEUE|$1|" $condorFile; shift
sed -i "s|NJOBS|$1|" $condorFile; shift

if [ ! -d "Resolutions_out" ]; then
  mkdir "Resolutions_out"
fi
if [ ! -d "Resolutions_err" ]; then
  mkdir "Resolutions_err"
fi

condor_submit $condorFile

rm $condorFile

condor_q
