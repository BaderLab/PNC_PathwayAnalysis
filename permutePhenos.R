# permutePhenos.R
# @description 
# Permutes phenotype status in the phenoFile. Subjects with status -9 are not included in the permutation
# Each permutation is written out to a file with the same format as phenoFile
# Note: pemutations are written to files to ensure reproducibility 
# TODO: Update paramters to be command line supplied, but currently set the following variables:
# numPerm: Total number of permutations
# startPermAt: Which number to start the permutations as (if the first n permutations has already been computed, set this to be n+1)
# rootDir: Set this to be the root directory
# Author: Shirley Hui

numPerm =1000
startPermAt = 101
rootDir <- "/home/shirleyhui/Work/PNC"
phenoFile <- paste(rootDir,"/refdata/pheno/case_control_status_180706.txt",sep="")
pheno <- read.delim(phenoFile,sep="\t",header=TRUE)
set.seed(105)

for (num in startPermAt:numPerm)
{
   print(sprintf("Making permutation %i",num))
   phenoPerm <- pheno[,1:2]
   for (i in 3:ncol(pheno))
   {
      #print(sprintf("Shuffling column %i",i))
      phenoCol <- pheno[,i]
      ix = which((phenoCol == 1 | phenoCol==2),arr.ind=TRUE)
      six = sample(ix)
      phenoCol[ix] = pheno[six,i]

      phenoPerm <- cbind(phenoPerm,phenoCol)
   }
   names(phenoPerm) <- names(pheno)
   write.table(phenoPerm,paste(rootDir,"/refdata/pheno/perm",numPerm,"-bin/pheno_bin_",num,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}
