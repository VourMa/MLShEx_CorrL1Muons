#!/bin/bash
condorFile=ScoutedDataCondor_$1_$2_$3.sub
cp ScoutedDataTemplate.sub $condorFile

if [ ! -f "CondorLog.txt" ]; then
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" > CondorLog.txt
else
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" >> CondorLog.txt
fi

chmod 755 $condorFile
echo "$condorFile created!"

sed -i "s|EXTRATEXT|$1|" $condorFile; shift

sed -i "s|QUEUE|$1|" $condorFile; shift
sed -i "s|NJOBS|$1|" $condorFile; shift

if [ ! -d "ScoutedData_out" ]; then
  mkdir "ScoutedData_out"
fi
if [ ! -d "ScoutedData_err" ]; then
  mkdir "ScoutedData_err"
fi

condor_submit $condorFile

rm $condorFile

condor_q
