#!/bin/bash

echo "Phenotype: $1"
echo "Start: $2, end: $3"
echo "NumCovar: $4"
jobFile=/home/shirleyhui/Work/PNC/jobfiles/plink_assoc_${1}_${2}-${3}-covar-${4}.txt
cat /dev/null > $jobFile;
for i in $(seq $2 $3); do
   echo "Rscript /home/shirleyhui/Work/PNC/runGWAS.R $1 $i $4" >> $jobFile 
done 

jobScript=/home/shirleyhui/Work/PNC/run_gwas_bin_${1}_${2}-${3}-covar-${4}.sh
cat /dev/null > $jobScript;
echo "#!/bin/bash" >> $jobScript
echo "mkdir $1-covar$4" >> $jobScript
echo "cd $1-covar$4" >> $jobScript
echo "parallel -j 10 --joblog run_gwas_bin_${1}_${2}-${3}.log < $jobFile" >> $jobScript
chmod u+x $jobScript
