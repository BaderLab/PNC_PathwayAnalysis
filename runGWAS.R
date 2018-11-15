# runGWAS.R
# @description
# Runs gwas for pheno of interest using logistic regression and 0-n PC covariates supplied by the user 
# @param pheno Phenotype to run gwas for
# @param permNum How many permutations to run gwas for
# @param numCovar Number of covariates to use with logistic regression (format: 0, no covariates or 1-n)
# Authors: Shirley Hui, Shraddha Pai

args = commandArgs(trailingOnly=TRUE)
pheno <- args[1]
permNum <- as.integer(args[2])
numCovar <- args[3]

cat(paste("* Running pheno:",pheno,", num perm:",permNum,",",numCovar,"\n"))

rootDir <- "/home/shirleyhui/Work/PNC/"
genoFile <- sprintf("%s/geno/PNC_imputed_merged.CLEAN_CEUTSI_5sd_FINAL",rootDir)
phenoFile <- sprintf("%s/refdata/pheno/case_control_status_180706.txt",rootDir)
if (permNum > 0) phenoFile <- sprintf("%s/refdata/pheno/perm1000-bin/pheno_bin_%i.txt",rootDir,permNum)
pcaFile <- sprintf("%s/geno/PNC_imputed_merged.CLEAN_CEUTSI_5sd_FINAL.pca.eigenvec",rootDir)

plink <- "/home/shirleyhui/plink_1.90/plink"

runGWAS <- function(usePheno,permNum =0,withCovariates=FALSE,covarNum) {
        cat(paste("* Running with covariates flag: ",withCovariates,sep=""))
        outDir <- sprintf("%s/%s-covar%s/",rootDir,usePheno,covarNum)
        baseOut <- sprintf("%s/%s",outDir,basename(genoFile))

	pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
	pheno <- pheno[,c("IID",sprintf("%s_ageSexCorr",usePheno))]

        # Read fam and pc files
        fam <- read.table(sprintf("%s.fam",genoFile),sep=" ",h=F,as.is=T)
        if (withCovariates) {
           pcdat <- read.delim(pcaFile,sep=" ",h=F,as.is=T)
        }
        pheno <- subset(pheno, IID %in% fam[,2])
        not_in_pheno <- setdiff(fam[,2], pheno$IID)
        not_in_pheno <- cbind(not_in_pheno,NA)
        colnames(not_in_pheno) <- names(pheno)
        pheno <- rbind(pheno, not_in_pheno)

        # Get ids that match between fam and pheno
        midx <- match(fam[,2],pheno$IID)
        if (all.equal(pheno$IID[midx],fam[,2])!=TRUE) {
           cat("pheno/fam didn't match."); browser()
        }
        pheno <- pheno[midx,]

        if (withCovariates) {
           # Make sure pcdat matches fam
           if (all.equal(pcdat[,2],fam[,2])!=TRUE) {
              cat("fam/pcdat don't match"); browser()
           }
           midx <- match(fam[,2],pcdat[,2])
           if (all.equal(pcdat[midx,2],fam[,2])!=TRUE) {
              cat("can't match fam/pcdat"); browser()
           }
           pcdat <- pcdat[midx,]
           covDat <- cbind(fam[,1:2],pcdat[,3:12])
        }

        # Write subject ids that have non NA phenotypes
        keepsamps <- pheno$IID[which(!is.na(pheno[,2]))]
        keepSampsFile <- sprintf("%s/keepsamps.perm%i.txt",outDir,permNum)
        write.table(cbind(0,keepsamps),file=keepSampsFile,sep=" ",col=F,row=F,quote=F)

	# writing binary status
	fam[,6] <- -9
	ctrl <- pheno$IID[which(pheno[,2]==1)]
	case <- pheno$IID[which(pheno[,2]==2)]
	fam[which(fam[,2]%in% case),6] <- 2
	fam[which(fam[,2] %in% ctrl),6] <- 1
	print(table(fam[,6]))

        if (withCovariates) {
           covFile <- sprintf("%s/%s.perm%i.cov",outDir,basename(genoFile),permNum)
           colnames(covDat) <- c("FID","IID",sprintf("PC%i",1:10))
           write.table(covDat,file=covFile,sep="\t",col=F,row=F,quote=F)
        }

        famFile <- sprintf("%s.perm%i.fam",baseOut,permNum)
        write.table(fam,file=famFile,sep="\t",col=F,row=F,quote=F)

        # gwas call
        if (withCovariates) {
           cmd <- sprintf("%s --bfile %s --fam %s --logistic --covar %s --covar-number %s --allow-no-sex --adjust --keep %s --out %s/%s_binAssoc.perm%i",plink,genoFile,famFile,covFile,covarNum,keepSampsFile,outDir,usePheno,permNum)
        } else {
           cmd <- sprintf("%s --bfile %s --fam %s --logistic --allow-no-sex --adjust --keep %s --out %s/%s_binAssoc.perm%i",plink,genoFile,famFile,keepSampsFile,outDir,usePheno,permNum)
        }
        print(cmd)
        system(cmd,intern=TRUE,wait=TRUE)
        browser()
        cat("* Cleaning up files...")
        unlink(sprintf("%s/%s_binAssoc.perm%i.assoc.logistic",outDir,usePheno,permNum))
        unlink(sprintf("%s/%s_binAssoc.perm%i.log",outDir,usePheno,permNum))
        unlink(sprintf("%s/%s_binAssoc.perm%i.nosex",outDir,usePheno,permNum))
        if (withCovariates) {
	   unlink(sprintf("%s/PNC_imputed_merged.CLEAN_CEUTSI_5sd_FINAL.perm%i.cov",outDir,permNum))
        }
	unlink(sprintf("%s/PNC_imputed_merged.CLEAN_CEUTSI_5sd_FINAL.perm%i.fam",outDir,permNum))
	unlink(sprintf("%s/keepsamps.perm%i.txt",outDir,permNum))
        unlink(sprintf("run_gwas_bin_%s_*.log",usePheno))

        adj_pvalue_outfile <- sprintf("%s/%s_binAssoc.perm%i.assoc.logistic.adjusted",outDir,usePheno,permNum)
        r = read.table(adj_pvalue_outfile,sep="")
        if (withCovariates) {
           write.table(r[2:dim(r)[1],2:3],sprintf("%s/%s_binAssoc.perm%i.assoc.logistic.adjusted.covar%s.txt",outDir,usePheno,permNum,covarNum),sep="\t",quote=F,row.names=F,col.names=F)
        } else {
           write.table(r[2:dim(r)[1],2:3],sprintf("%s/%s_binAssoc.perm%i.assoc.logistic.adjusted.nocovar.txt",outDir,usePheno,permNum),sep="\t",quote=F,row.names=F,col.names=F)
        }
        unlink(adj_pvalue_outfile)
}
if (numCovar == 0) {
   withCovariates = FALSE
} else {
   withCovariates = TRUE
}
runGWAS(pheno,permNum,withCovariates,numCovar)

