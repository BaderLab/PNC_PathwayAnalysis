'''
This script performs enrichment score normalization given enrichment score output files, one for each permutation required.
Set parameters before running:
ROOT_DIR: The root directory
Parameters:
phenotype: The phenotype 
type: Set to be bin (for binary)
covars: Number of covariates used 
numPerm: Number of permutations
Authors: Shirley Hui, Asha Rostamianfar
'''
import util
import reader
import operator
import sys
import glob
import numpy as np
import copy
import os
minSize = 10
maxSize = 200

phenotype = str(sys.argv[1])
type=str(sys.argv[2])
covars=str(sys.argv[3])
numPerm=str(sys.argv[4])
print phenotype

ROOT_DIR = "/home/shirleyhui/Work/PNC/"
if type=='bin': 
   if covars=='0':
      originalPesScore = ROOT_DIR+'gsea/bin/gsea_'+phenotype+'_nocovar/obs/'+phenotype+'_binAssoc.perm0.assoc.logistic.adjusted.nocovar.txt.May012018_gmt-gsea.txt'
      inDir = ROOT_DIR+'gsea/bin/gsea_'+phenotype+'_nocovar/perm/*.txt'
   else:
      originalPesScore = ROOT_DIR+'gsea/bin/gsea_'+phenotype+'_covars'+covars+'/obs/'+phenotype+'_binAssoc.perm0.assoc.logistic.adjusted.covar'+covars+'.txt.May012018_gmt-gsea.txt'
      inDir = ROOT_DIR+'gsea/bin/gsea_'+phenotype+'_covars'+covars+'/perm/*.txt'

pesFiles = glob.glob(inDir)
print len(pesFiles)

def addPES(pesFile, pathwayToPes, pathwayToSize, firstTime = False):
  for line in open(pesFile):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    pes = float(parts[1])
    size = int(parts[2])
    if size < minSize or size > maxSize:
      continue
    if firstTime:
      pathwayToSize[pId] = size
      pathwayToPes[pId] = [pes]
    else:
      try:
         pathwayToPes[pId].append(pes)
      except KeyError,e:
         pass
def getAllPes():
  pathwayToSize = {}
  pathwayToPes = {}
  addPES(originalPesScore, pathwayToPes, pathwayToSize, True)
  i = 0
  for f in pesFiles:
    print 'Analyzing', f
    if i % 100 == 0:
      print i
      sys.stdout.flush()
    i += 1
    addPES(f, pathwayToPes, pathwayToSize, False)
  return pathwayToPes, pathwayToSize

def getNES(pathwayToPes, pIdx, pathwayToNES, firstTime = False):
  for p, pesScores in pathwayToPes.iteritems():
    #print p
    mean = np.mean(pesScores[1:])
    std = np.std(pesScores[1:])
    if pIdx >= len(pesScores):
       continue;
    if std == 0:
      nes = 0
    else:
      nes = ((pesScores[pIdx] - mean) / std)
    if firstTime:
      #print p, nes
      pathwayToNES[p] = nes
    else:
      pathwayToNES[p].append(nes)

pathwayToPes, pathwayToSize = getAllPes()
firstTime = True

outDir =( '%s/nes/%s/' % (GSEA_DIR,phenotype))
print 'Writing out to directory:', outDir
if not os.path.exists(outDir):
  print "Out dir doesn't exist. trying to make"
  os.makedirs(outDir)

stop = int(numPerm)+1
for i in range(0,stop):
  pathwayToNES = {}
  getNES(pathwayToPes, i, pathwayToNES, firstTime)
  outFileName =  outDir + ('nes-%s-%s-covars%s-%d.txt' % (phenotype,type,covars,i))
  util.printMapToFile(pathwayToNES, outFileName)

