* Instructions on how to run the Pathway Analysis pipeline for PNC data
* Nov 9, 2018
* By Shirley Hui

A) PERMUTE PHENOTYPES

1) Permute phenotypes

File: ROOT_DIR/permutePhenos.R
- This file permutes phenotype status in the phenoFile (see code). 
- Subjects with status -9 are not included in the permutation

Run: 
> Rscript perumtePhenos.R 

 * Note: Set variables inside script before running

Output: 
 * Each permutation is written out to a file with the same format as phenoFile
 * Note: pemutations are written to files to ensure reproducibility 

B) PERFORM GWAS (genome-wide associations)

1) Make gwas batch and submit scripts

File: ROOT_DIR/make_plink_assoc_bin_jobfiles.sh
- This script make a txt file that contains the call to run gwas (logistic regression + PC covariates) to compute snp-phenotype associations using plink.

Run:
> ./make_plink_assoc_bin_jobfiles.sh <phenotype> 0 <num_perms> <covars>

Example: ./make_plink_assoc_bin_jobfiles.sh volt_svt 0 100 1-5

Output:
 * a file called ROOT_DIR/jobfiles/plink_assoc_<phenotype>_0-<num_perms>-covar-<covars>.txt that contains the call to run one permutation of gwas (one call per line)
 * a file called ROOT_DIR/run_gwas_bin_<phenotype>_0-<num_perms>-covar-<covars>.sh that calls the txt file above via parallel

2) Compute gene-phenotype associations 
- Run ROOT_DIR/run_gwas_bin_<phenotype>_0-<num_perms>-covar-<covars>.sh script from previous
- The script runs GWAS (runGWAS.R) 0 to <num_perms> times for pheno of interest using logistic regression and 0-n PC covariates supplied by the user 

Output:
 * pvalues for each SNP outputted to a file for each permuation (0 to <num_perms>)

C) PERFORM PATHWAY GSEA

The following assumes the pipeline is running on a high 

1) Make gsea batch and submit scripts

File: ROOT_DIR/code/src/MakeJobs.java
- This script will make a bunch of .sh files that contain the 8 parallel jobs that will run when it is called.
- If num perm = 100, 13 .sh files will be made.  If num perm= 1000, 125 .sh files will be made.  Note code doesn't really support num perms not equal to 100 or 1000.
- If you make any code changes, to compile script, run the ./compile.sh at the prompt

Run: 
> module load java
> java -cp ../out MakeJobs <phenotype> <num_perms> bin <covars>

Example: java -cp ../out MakeJobs volt_svt 100 bin 1-5

Output:
  * a bunch of job files located in ROOT_DIR/job_scripts/gsea/run-gsea-<phenotype>-<perm>-bin-covars<covars>.sh
  * a master submit script located in ROOT_DIR/submit-gsea-<phenotype>-bin-covars<covars>.sh

2) Run gsea 
- Run submit script from previous
- You will get a bunch of warnings (which is ok) about how the submitted job time is the minimum of 15 mins.
- This will create several directories (one for each permutation).  In each directory there will be two output files (*gsea.txt - gsea output, *genes.txt - detailed output)

3) Move gsea output files into one directory

File: ROOT_DIR/moveFiles.sh
- This will move all output files into ROOT_DIR/gsea/bin/gsea-<phenotype>_covars<covars> where:
  * gsea.txt files get moved to obs directory,for perm0 or perm directory, for all other permutations
  * genes.txt files get moved to genes directory for all perms

Run: 
> moveFiles.sh <phenotype> bin <covars> <numperms>

Example: moveFiles.sh volt_svt bin 1-5 100

- Delete empty directories and run txt and log files:
> rm -r Gsea-<phenotype>_bin-covars<covars>_*
> rm run_gsea_<phenotype>-<batch_num>-bin-covars<covars>*

4) Normalize ES scores

File: ROOT_DIR/mygseaCalcNES.py
- This script will normalize the enrichment scores from step 2 above.

Run: 
> module load python/2.7.14-anaconda5.1.0
> python mygseaCalcNES.py <phenotype> bin <covars> <numperm>   
  
Example: python mygseaCalcNES.py volt_svt bin 1-5 100

Output:
- This will create a bunch of nes*.txt in the ROOT_DIR/nes/<phenotype> directory

5) Compute FDR score

File: ROOT_DIR/compareFDR.R
- This script will compute FDR scores using NES from step 4)

Run: 
> module load r/3.4.3-anaconda5.1.0
> Rscript computeFDR.R <phenotype> bin <covars> <numperm>
  
Example: Rscript computeFDR.R volt_svt bin 1-5 8

Output: - This will create a file fdr-<phenotype>-bin-covars<covars>-<numperm>.txt in the ROOT_DIR/fdr directory
