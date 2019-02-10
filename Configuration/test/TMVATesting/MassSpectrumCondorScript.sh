#!/bin/bash
condorFile=MassSpectrumCondor_$1_$2_$3_$4_$5.sub
cp MassSpectrumTemplate.sub $condorFile

if [ ! -f "CondorLog.txt" ]; then
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" > CondorLog.txt
else
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" >> CondorLog.txt
fi

chmod 755 $condorFile
echo "$condorFile created!"

sed -i "s|DATASET|$1|" $condorFile; shift
sed -i "s|ERA|$1|" $condorFile; shift
sed -i "s|TFSPLIT|$1|" $condorFile; shift

sed -i "s|QUEUE|$1|" $condorFile; shift
sed -i "s|NJOBS|$1|" $condorFile; shift

if [ ! -d "MassSpectrum_out" ]; then
  mkdir "MassSpectrum_out"
fi
if [ ! -d "MassSpectrum_err" ]; then
  mkdir "MassSpectrum_err"
fi

condor_submit $condorFile

rm $condorFile

condor_q
