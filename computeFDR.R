# computeFDR.R
# @description
# Compute the FDRs for a given phenotype, assumes gsea pipeline has been run previously (i.e. gsea + nes)
# TODO: Update to remove type as it is always binary in this case
# Set rootdir to point to root directory of user.
# @param phenotype The phenotype to compute FDRs for
# @param type Type of regression that was run (note. always set to bin for binary)
# @param covars How many covariates were used for regression (format: 0, for no covariates, 1-n for n covariates)
# @param numPerm Number of permutations gsea was run for
# Author: Shirley Hui 
args <- commandArgs(TRUE)

phenotype = args[1]
type = args[2]
covars = args[3]
numPerm = args[4]

cat(paste("* Processing phenotype: ",phenotype,",",type,",",covars,",",numPerm,sep=""))
rootdir = "/home/shirleyhui/Work/PNC/"

nesMat <- c()
for (i in 0:numPerm)
 {
      nesFile <- paste(rootdir,"nes/",phenotype,"/nes-",phenotype,"-",type,"-covars",covars,"-",i,".txt",sep="")
      r <- read.delim(nesFile,sep="\t",header=F)
      nesMat <- cbind(nesMat,r[,2])
 }
pathwayNames <- r[,1]
totalObs = length(which(nesMat[,1]>=0))
nesObs = nesMat[,1]
nesAll = nesMat[,2:dim(nesMat)[2]]
totalAll = length(which(nesAll >=0))
numGtObs = c()
fdrs = c()
pix = 1
for (nesStar in nesObs)
{   
    ix = which(nesObs>=0)
    nesObs = nesObs[ix]
    ix = which(nesObs >= nesStar)
    numGtObs = length(ix)
    numGtAll = 0
    for (i in 1:dim(nesAll)[1])
    {
       ix = which(nesAll[i,]>=0)
       ix = which(nesAll[i,ix] >= nesStar)
       numGtAll = numGtAll + length(ix)      
    } 
    percentGtObs = numGtObs/totalObs
    percentGtAll = numGtAll/totalAll
    fdr = percentGtAll/percentGtObs
    fdrs = c(fdrs,fdr)
    print(sprintf("%s\t%1.2f",as.character(pathwayNames[pix]),fdr))
    pix=pix+1
}
outFile = paste(rootdir,"fdr/fdr-",phenotype,"-",type,"-covars",covars,"-",numPerm,".txt",sep="")
results = cbind(as.character(pathwayNames),fdrs)
rix = order(results[,2])      
write.table(results[rix,],outFile,quote=F,row.names=F,col.names=F,sep="\t")
