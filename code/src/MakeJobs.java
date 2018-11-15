import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: shirleyhui
 * Date: Jul 22, 2016
 * Time: 1:18:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class MakeJobs{
    private static String SRC_DIR = "/home/shirleyhui/Work/PNC/";
    private static String ROOT_DIR = "/home/shirleyhui/Work/PNC/";
    public MakeJobs()
    {
    }
    
    public static void main(String[] args)
    {
        MakeJobs j = new MakeJobs();

        /* Make GSEA job and submit files */
        String phenotype = args[0];
        int numPerm = Integer.parseInt(args[1]);
        String type = args[2];
        String covars = args[3];
        try
        {
           String dirName = "";
           if (!covars.equals("0"))
              dirName = ROOT_DIR + "pvalues/perm"+numPerm + "/bin/"+covars;
           else
              dirName = ROOT_DIR + "pvalues/perm"+numPerm + "bin";
           
           System.out.println("Looking in: " + dirName);
           File dir = new File(dirName);
           File[] files = dir.listFiles();
           for (int jj = 0 ; jj < files.length;jj++)
           {
              String filename = files[jj].getName();
              if (filename.startsWith(phenotype))
              {
                 int startIx = 1;
                 //int ix = filename.indexOf("_");
                 int numTimes = 125;
                 if (numPerm==100) numTimes = 13;
	         for (int i = 0; i < numTimes;i++)
                 {
	            j.make_gsea(startIx,numPerm,phenotype,type,covars);
                    startIx = startIx + 8;
                 }
                j.make_gsea_submit_script(phenotype,type,numPerm,covars);
              }
            }
         }
         catch(Exception e)
         {
            System.out.println("Exception: " + e);
            e.printStackTrace();
         }
	
    }
    public void make_gsea(int startIx,int numPerm, String phenotype,String type,String covars)
    {
        try
        {
            String outfilename = ROOT_DIR+"jobfiles/gsea/run-gsea-"+phenotype+"-"+startIx+"-"+type+"-covars"+covars+".sh";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            bw.write("#!/bin/bash\n");
            bw.write("#SBATCH --nodes=1\n");
            bw.write("#SBATCH --cpus-per-task=1\n");
            bw.write("#SBATCH --ntasks=8\n");
            bw.write("#SBATCH --time=0:15:00\n");
            bw.write("#SBATCH --job-name run_gsea_"+phenotype+"-"+startIx+"-"+type+"-covars"+covars+"\n");
            bw.write("#SBATCH --output=run_gsea_"+phenotype+"-"+startIx+"-"+type+"-covars"+covars+"_%j.txt\n");
            bw.write("module load gnu-parallel/20180322\n");
            bw.write("module load python/2.7.14-anaconda5.1.0\n");
            bw.write("parallel -j 8 --joblog run_gsea_"+phenotype+"-"+startIx+"-"+type+"-covars"+covars+".log <<EOF\n");

           
            int endAt = startIx+8;
            if (startIx == 1)
            {
               if (type.equals("bin"))
               {
                  if (covars.equals("0")) {
                     bw.write("mkdir Gsea-"+phenotype+"_"+type+"-nocovar_0; cd Gsea-"+phenotype+"_"+type+"-nocovar_0; python "+SRC_DIR+"mygseaCalcESCentrality-perm-pvalues-out.py "+ROOT_DIR+"pvalues/perm"+numPerm+"/bin/"+covars+"/"+phenotype+"/"+phenotype+"_binAssoc.perm0.assoc.logistic.adjusted.nocovar.txt; echo \"job 0 finished\"\n");
                  } else {
                     bw.write("mkdir Gsea-"+phenotype+"_"+type+"-covars"+covars+"_0; cd Gsea-"+phenotype+"_"+type+"-covars"+covars+"_0; python "+SRC_DIR+"mygseaCalcESCentrality-perm-pvalues-out.py "+ROOT_DIR+"pvalues/perm"+numPerm+"/bin/"+covars+"/"+phenotype+"/"+phenotype+"_binAssoc.perm0.assoc.logistic.adjusted.covar"+covars+".txt; echo \"job 0 finished\"\n");
                  }
               }
            }
            int  numIt = 8;            
            if (endAt > numPerm)
               numIt = numPerm-startIx+1;  
              
            for (int i = 0;i< numIt;i++)
            {
               if (type.equals("bin")) {
                  if (covars.equals("0")) {
                     bw.write("mkdir Gsea-"+phenotype+"_"+type+"-nocovar_"+(startIx+i)+"; cd Gsea-"+phenotype+"_"+type+"-nocovar_"+(startIx+i) +"; python "+SRC_DIR+"mygseaCalcESCentrality-perm-pvalues-out.py "+ROOT_DIR+"pvalues/perm"+numPerm+"/bin/"+covars+"/"+phenotype+"/"+phenotype+"_binAssoc.perm"+(startIx+i)+".assoc.logistic.adjusted.nocovar.txt; echo \"job "+(startIx+i)+" finished\"\n");
                  } else {
                     bw.write("mkdir Gsea-"+phenotype+"_"+type+"-covars"+covars+"_"+(startIx+i)+"; cd Gsea-"+phenotype+"_"+type+"-covars"+covars+"_"+(startIx+i) +"; python "+SRC_DIR+ "mygseaCalcESCentrality-perm-pvalues-out.py "+ROOT_DIR+"pvalues/perm"+numPerm+"/bin/"+covars+"/"+phenotype+"/"+phenotype+"_binAssoc.perm"+(startIx+i)+".assoc.logistic.adjusted.covar"+covars+".txt; echo \"job "+(startIx+i)+" finished\"\n");
                  }
               }
            }
            bw.write("EOF");
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public void make_gsea_submit_script(String phenotype,String type,int numPerm, String covars)

    {
       try
       {
          String outfilename = ROOT_DIR + "submit-gsea-"+phenotype+"-"+type+"-covars"+covars+".sh";
          BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

          String dirName = ROOT_DIR + "jobfiles/gsea/";
          File dirFile = new File(dirName);
          File[] files =dirFile.listFiles();
          for (int i = 0 ; i < files.length;i++)
          {
             File file = files[i];

             String filename = file.getName();
             if (filename.endsWith("covars"+covars+".sh"))
	     {
                String [] splitfilename = filename.split("-");
                String filename0 = splitfilename[2];
       	        if (filename0.equals(phenotype))
                {  
 		   System.out.println(filename);
                   bw.write("sbatch "+ROOT_DIR+"jobfiles/gsea/"+filename+"\n");
	        }
             }
          }
          bw.close();
       }
       catch(Exception e)
       {
          System.out.println("Exception: " + e);
          e.printStackTrace();
       } 
    }

}
