#!/bin/bash
condorFile=TMVARegressionCondor_$1_$2_$3_$4_$5.sub
cp TMVARegressionTemplate.sub $condorFile

if [ ! -f "CondorLog.txt" ]; then
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" > CondorLog.txt
else
	echo "$condorFile sent at $(date +'%T') on $(date +'%d-%m-%y')" >> CondorLog.txt
fi

chmod 755 $condorFile
echo "$condorFile created!"

sed -i "s|TF|$1|" $condorFile; shift
sed -i "s|ERA|$1|" $condorFile; shift
sed -i "s|GUYS|$1|" $condorFile; shift

sed -i "s|QUEUE|$1|" $condorFile; shift
sed -i "s|NJOBS|$1|" $condorFile; shift

if [ ! -d "TMVARegression_out" ]; then
  mkdir "TMVARegression_out"
fi
if [ ! -d "TMVARegression_err" ]; then
  mkdir "TMVARegression_err"
fi

condor_submit $condorFile

rm $condorFile

condor_q
