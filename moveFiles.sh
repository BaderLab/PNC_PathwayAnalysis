#!/bin/bash
phenotype=$1
contBin=$2
covars=$3
numperm=$4

echo "Phenotype: ${phenotype}"
echo "Type: ${contBin}"
echo "Covars: ${covars}"
echo "Num perm: ${numperm}"

if [ $covars = "0" ]
then
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_nocovar/perm
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_nocovar/obs
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_nocovar/genes
else
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_covars${3}/perm
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_covars${3}/obs
   mkdir -p  gsea/${contBin}/gsea_${phenotype}_covars${3}/genes
fi

for (( i=0; i<=$numperm; i++ ))
do
    if [ $i -eq 0 ] 
    then
       if [ $covars = "0" ]
       then
          mv Gsea-${phenotype}_${contBin}-nocovar_${i}/*gsea.txt gsea/${contBin}/gsea_${phenotype}_nocovar/obs
          mv Gsea-${phenotype}_${contBin}-nocovar_${i}/*genes.txt gsea/${contBin}/gsea_${phenotype}_nocovar/genes
       else
          mv Gsea-${phenotype}_${contBin}-covars${covars}_${i}/*gsea.txt gsea/${contBin}/gsea_${phenotype}_covars${3}/obs
          mv Gsea-${phenotype}_${contBin}-covars${covars}_${i}/*genes.txt gsea/${contBin}/gsea_${phenotype}_covars${3}/genes
       fi
    else
       if [ $covars = "0" ]
       then
          mv Gsea-${phenotype}_${contBin}-nocovar_${i}/*gsea.txt gsea/${contBin}/gsea_${phenotype}_nocovar/perm
          mv Gsea-${phenotype}_${contBin}-nocovar_${i}/*genes.txt gsea/${contBin}/gsea_${phenotype}_nocovar/genes
       else
          mv Gsea-${phenotype}_${contBin}-covars${covars}_${i}/*gsea.txt gsea/${contBin}/gsea_${phenotype}_covars${3}/perm
          mv Gsea-${phenotype}_${contBin}-covars${covars}_${i}/*genes.txt gsea/${contBin}/gsea_${phenotype}_covars${3}/genes
       fi
    fi
done

if [ $covars = "0" ]
then
   v=$(ls gsea/${contBin}/gsea_${phenotype}_nocovar/obs/*.txt | wc -l)
   echo "Num files in obs: " ${v}
   v=$(ls gsea/${contBin}/gsea_${phenotype}_nocovar/perm/*.txt | wc -l)
   echo "Num files in perm: " ${v}    
   v=$(ls gsea/${contBin}/gsea_${phenotype}_nocovar/genes/*.txt | wc -l)
   echo "Num files in genes: " ${v}    
else
   v=$(ls gsea/${contBin}/gsea_${phenotype}_covars${3}/obs/*.txt | wc -l)
   echo "Num files in obs: " ${v}
   v=$(ls gsea/${contBin}/gsea_${phenotype}_covars${3}/perm/*.txt | wc -l)
   echo "Num files in perm: " ${v}
   v=$(ls gsea/${contBin}/gsea_${phenotype}_covars${3}/genes/*.txt | wc -l)
   echo "Num files in genes: " ${v}
fi
