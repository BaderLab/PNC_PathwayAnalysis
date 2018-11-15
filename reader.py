import sys
#import scipy.stats as st
#from scipy.stats import chisqprob
#import networkx as nx

ALL_STATUSES = ['0.1.777', '0.2.777', '0.4.777', '0.5.777', '0.888.777', '1.777.0', '1.777.1', '1.777.2', '1.777.888', '2.777.0', '2.777.1', '2.777.2', '2.777.888', '3.777.0', '3.777.1', '3.777.888']
ALLELE_MAP = {'AA': '0', 'AC': '1', 'AT': '2', 'AG': '3', 'CG': '4', 'CC': '5', 'TT':'6', 'GG':'7', 'ID':'8', 'DI':'9', 'II':'A', 'DD':'B', '00':'C'}
REVERSE_ALLELE_MAP = {}
for a, i in ALLELE_MAP.iteritems():
  REVERSE_ALLELE_MAP[i] = a

ETHMAP = {'European': '0', 'Asian': '1', 'African': '2'}
ERSTATUSMAP = {'0.00': '0', '1.00': '1', '' : '9', '888.00': '8'}

class PatientColumns():
  PID = 0
  STUDY = 2
  STATUS = 3
  CTRTYPE = 5
  STUDYTYPE = 6
  AGE = 10
  ETH_GENO = 12
  PC1_EURO = 20
  PC8_EURO = 27
  PC_LMBC = 28
  PC1_ASIAN = 29
  PC2_ASIAN = 30
  PC1_AFRICAN = 31
  PC2_AFRICAN = 32
  DIAG_AGE = 33
  ER = 36


def readFromMAF(fileName = 'data/snpCountsReference.txt.ALL_SNPS_CONTROLS_MAF.txt'):
  snpToMAF = {}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    snpToMAF[parts[0]] = float(parts[3])
  return snpToMAF


def getStudyNumbers(fileName = 'refdata/studies.txt'):
  studyToNumMap = {}
  lN = 0
  for line in open(fileName):
    studyToNumMap[line.replace('\n', '')] = lN
    lN += 1
  return studyToNumMap


def readFromPathwayWithGeneSums(pfile):
  pathwayToGeneSums = {}
  for line in open(pfile):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    pathwayToGeneSums[pId] = {}
    for ps in parts[1:]:
      if len(ps) == 0:
        continue
      pd = ps.split(':')
      geneId = pd[0]
      geneSum = pd[1]
      pathwayToGeneSums[pId][geneId] = float(geneSum)
  return pathwayToGeneSums


def readFromPatientFile(fileName = 'refdata/2013_simard/concept_263_andrulis_pheno.txt'):
  pIdMap = {}
  pStatusMap = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue
    line = line.replace('\n', '')
    parts = line.split('\t')
    statusId = parts[PatientColumns.STATUS] + "." + parts[PatientColumns.CTRTYPE] + "." + parts[PatientColumns.STUDYTYPE] + "." + ETHMAP[parts[PatientColumns.ETH_GENO]] + "." + ERSTATUSMAP[parts[PatientColumns.ER]]
    pList = []
    if statusId in pStatusMap:
      pList = pStatusMap[statusId]
    pList.append(parts[PatientColumns.PID])
    pStatusMap[statusId] = pList
    assert(parts[PatientColumns.PID] not in pIdMap)
    pIdMap[parts[PatientColumns.PID]] = statusId
  return pIdMap, pStatusMap


def readFromPatientFileDetails(fileName = 'refdata/concept_263_andrulis_pheno.txt'):
  pIdMap = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue
    line = line.replace('\n', '')
    parts = line.split('\t')
    pIdMap[parts[PatientColumns.PID]] = parts
  return pIdMap


def readPIdStatuses(fileName = '/home/g/gbader/asha/allPIDsAnnotated.txt'):
  statusIdToPIdsMap = {}
  #for s in ALL_STATUSES:
  #  statusIdToPIdsMap[s] = []
  pIdToStatusIdMap = {}
  for line in open(fileName):
    parts = line.split('\t')
    statusId = parts[0]
    if statusId not in statusIdToPIdsMap:
      statusIdToPIdsMap[statusId] = []
    pIds = parts[1].split(',')
    pIds[-1] = pIds[-1].replace('\n', '')
    for pId in pIds:
      pIdToStatusIdMap[pId] = statusId
      statusIdToPIdsMap[statusId].append(pId)
  return pIdToStatusIdMap, statusIdToPIdsMap


def readLineNumToPId(fileName = '/home/g/gbader/asha/lineNumToPid.txt'):
  lineNumToPid = {}
  for line in open(fileName):
    line = line.replace('\n', '')
    parts = line.split('\t')
    lineNumToPid[int(parts[0])] = parts[1]
  return lineNumToPid


def readPackedAlleleFile(pIdGroups, fileName = '/scratch/g/gbader/asha/packedBader.txt', shouldBreak = False):
  lineNumToPid = readLineNumToPId()
  snpCountsPerGroup = {}
  pIdToGroups = {}
  for group, pIds in pIdGroups.iteritems():
    for pId in pIds:
      pIdToGroups[pId] = group
    snpCountsPerGroup[group] = {}
  headings = []
  lineNum = 0
  for line in open(fileName):
    if lineNum % 10000 == 0:
      sys.stdout.write(str(lineNum) + '\n')
    lineNum += 1
    if lineNum == 1:
      headings = line.split('\t')
      headings[-1] = headings[-1].replace('\n', '')
      headings = headings[1:] # remove PId
      continue
    pId = lineNumToPid[lineNum]
    # check if PId is even in the list
    if pId not in pIdToGroups:
      continue
    group = pIdToGroups[pId]
    snpMap = snpCountsPerGroup[group]
    snpNum = -1
    line = line[:-1] # remove '\n'
    for a in line[10:]:
      snpNum += 1
      snp = headings[snpNum]
      if snp not in snpMap:
        snpMap[snp] = {}
      if a not in snpMap[snp]:
        snpMap[snp][a] = 0
      snpMap[snp][a] += 1
    if shouldBreak:
      break
  return snpCountsPerGroup


def readPackedAlleleFilePerSnp(pIds, snps, fileName = '/scratch/g/gbader/asha/packedBader.txt'):
  lineNumToPid = readLineNumToPId()
  pIdToSnps = {}
  headings = []
  lineNum = 0
  for line in open(fileName):
    if lineNum % 10000 == 0:
      sys.stdout.write(str(lineNum) + '\n')
    lineNum += 1
    if lineNum == 1:
      headings = line.split('\t')
      headings[-1] = headings[-1].replace('\n', '')
      headings = headings[1:] # remove PId
      continue
    pId = lineNumToPid[lineNum]
    # check if PId is even in the list
    if pId not in pIds:
      continue
    pIdToSnps[pId] = {}
    snpNum = -1
    line = line[:-1] # remove '\n'
    for a in line[10:]:
      snpNum += 1
      snp = headings[snpNum]
      if snp not in snps:
        continue
      pIdToSnps[pId][snp] = a
  return pIdToSnps


def readPackedAlleleFilePerSnpMaf(pIds, snps, snpMaj, pIdToListIdx, snpToListIdx, toFill, fileName = '/scratch/g/gbader/asha/packedBader.txt'):
  lineNumToPid = readLineNumToPId()
  #pIdToSnps = {}
  headings = []
  lineNum = 0
  for line in open(fileName):
    if lineNum % 10000 == 0:
      sys.stdout.write(str(lineNum) + '\n')
      sys.stdout.flush()
    lineNum += 1
    if lineNum == 1:
      headings = line.split('\t')
      headings[-1] = headings[-1].replace('\n', '')
      headings = headings[1:] # remove PId
      continue
    pId = lineNumToPid[lineNum]
    # check if PId is even in the list
    if pId not in pIds:
      continue
    #pIdToSnps[pId] = {}
    snpNum = -1
    line = line[:-1] # remove '\n'
    for a in line[10:]:
      snpNum += 1
      snp = headings[snpNum]
      if snp not in snps:
        continue
      actualA = REVERSE_ALLELE_MAP[a]
      maj = snpMaj[snp][0]
      acount = 0
      if actualA[0] == maj:
        acount += 1
      if actualA[1] == maj:
        acount += 1
      if actualA == '00':
        acount = 'NA'
      #pIdToSnps[pId][snp] = acount
      pIdx = pIdToListIdx[pId]
      snpIdx = snpToListIdx[snp]
      toFill[pIdx][snpIdx] = acount
  #return pIdToSnps


def readFromAlleleCountFile(snpFilters, fileName):
  allSnpInfo = {}
  lineNum = 0
  statusId = ''
  for line in open(fileName, 'r'):
    lineNum +=1
    line = line.replace('\n', '')
    if line[0] == '=':
      statusId = line.split(' ')[1]
      allSnpInfo[statusId] = {}
      continue
    if statusId == '':
      continue
    parts = line.split('\t')
    snpId = parts[0]
    if len(snpFilters) > 0 and snpId not in snpFilters:
      continue
    allSnpAlleles = allSnpInfo[statusId]
    alleleF = {}
    for a in REVERSE_ALLELE_MAP:
      alleleF[a] = 0  # initiate with 0 so that we can easily calculate the difference
    for p in range(len(parts)):
      if p == 0 or len(parts[p]) < 2:
        continue
      ac = parts[p].split(':')
      assert(ac[0] in alleleF)
      alleleF[ac[0]] = int(ac[1])
    allSnpAlleles[snpId] = alleleF
  return allSnpInfo


def readFromReferenceGenome(fileName):
  COLUMNS = {'CHR': 2, 'TX_START': 4, 'TX_END': 5, 'NAME': 12}
  geneInfo = {}
  geneInfoPerChr = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue # skip header
    parts = line.replace('\n', '').split('\t')
    geneName = parts[COLUMNS['NAME']]
    chrNum = parts[COLUMNS['CHR']]
    txStart = int(parts[COLUMNS['TX_START']])
    txEnd = int(parts[COLUMNS['TX_END']])
    if geneName in geneInfo:
      assert(geneInfo[geneName][0] == chrNum)
      txStart = min(txStart, geneInfo[geneName][1])
      txEnd = max(txEnd, geneInfo[geneName][2])
    geneInfo[geneName] = [chrNum, txStart, txEnd]
    if chrNum not in geneInfoPerChr:
      geneInfoPerChr[chrNum] = {}
    geneInfoPerChr[chrNum][geneName] = [chrNum, txStart, txEnd]
  return geneInfo, geneInfoPerChr


def readFromReferenceGenomeAnnovar(fileName):
  COLUMNS = {'CHR': 2, 'TX_START': 4, 'TX_END': 5, 'NAME': 12}
  geneInfo = {}
  geneInfoPerChr = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue # skip header
    parts = line.replace('\n', '').split('\t')
    geneName = parts[COLUMNS['NAME']].upper()
    chrNum = parts[COLUMNS['CHR']]
    if '_' in chrNum:
      continue
    txStart = int(parts[COLUMNS['TX_START']])
    txEnd = int(parts[COLUMNS['TX_END']])
    if geneName in geneInfo:
      if geneInfo[geneName][0] == chrNum:
        #if abs(geneInfo[geneName][1] - txStart) > 1000000:
        #  print '*********', geneName
        txStart = min(txStart, geneInfo[geneName][1])
        txEnd = max(txEnd, geneInfo[geneName][2])
        #print '**********', geneName
      else:
        pass
    geneInfo[geneName] = [chrNum, txStart, txEnd]
    if chrNum not in geneInfoPerChr:
      geneInfoPerChr[chrNum] = {}
    geneInfoPerChr[chrNum][geneName] = [chrNum, txStart, txEnd]
  return geneInfo, geneInfoPerChr

def readFromReferenceGenomeAnnovarSeparated(fileName):
  COLUMNS = {'CHR': 2, 'TX_START': 4, 'TX_END': 5, 'NAME': 12}
  geneInfo = {}
  geneInfoPerChr = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue # skip header
    parts = line.replace('\n', '').split('\t')
    geneName = parts[COLUMNS['NAME']].upper()
    chrNum = parts[COLUMNS['CHR']]
    if '_' in chrNum:
      continue
    txStart = int(parts[COLUMNS['TX_START']])
    txEnd = int(parts[COLUMNS['TX_END']])
    if geneName not in geneInfo:
      geneInfo[geneName] = []
    geneInfo[geneName].append([chrNum, txStart, txEnd])
    if chrNum not in geneInfoPerChr:
      geneInfoPerChr[chrNum] = {}
    if geneName not in geneInfoPerChr[chrNum]:
      geneInfoPerChr[chrNum][geneName] = []
    geneInfoPerChr[chrNum][geneName].append([chrNum, txStart, txEnd])
  return geneInfo, geneInfoPerChr


def readFromReferenceGenomeUCSCEntrez(fileName):
  COLUMNS = {'CHR': 1, 'TX_START': 3, 'TX_END': 4, 'SYMBOL': 7, 'ENTREZ': 9}
  geneInfo = {}
  geneInfoPerChr = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    lineNum += 1
    if lineNum == 1:
      continue # skip header
    parts = line.replace('\n', '').split('\t')
    geneSymbol = parts[COLUMNS['SYMBOL']]
    chrNum = parts[COLUMNS['CHR']]
    txStart = int(parts[COLUMNS['TX_START']])
    txEnd = int(parts[COLUMNS['TX_END']])
    entrezId = parts[COLUMNS['ENTREZ']]
    if entrezId == 'n/a':
      continue
    if entrezId in geneInfo:
      if geneInfo[entrezId][0] != chrNum:
        pass #print line.replace('\n', '')
      #assert(geneInfo[entrezId][0] == chrNum)
      txStart = min(txStart, geneInfo[entrezId][1])
      txEnd = max(txEnd, geneInfo[entrezId][2])
    geneInfo[entrezId] = [chrNum, txStart, txEnd]
    if chrNum not in geneInfoPerChr:
      geneInfoPerChr[chrNum] = {}
    geneInfoPerChr[chrNum][entrezId] = [chrNum, txStart, txEnd]
  return geneInfo, geneInfoPerChr


def readFromReferenceGenomeNCBIEntrezEntToEns(fileName):
  COLUMNS = {'ENTREZ': 1, 'SYMBOL': 2, 'DBXREFS': 5, 'CHR': 6}
  geneInfo = {}
  lineNum = 0
  for line in open(fileName, 'r'):
    line = line.replace('\n', '')
    lineNum += 1
    if lineNum == 1:
      continue # skip header
    parts = line.replace('\n', '').split('\t')
    geneSymbol = parts[COLUMNS['SYMBOL']]
    entrezId = parts[COLUMNS['ENTREZ']]
    dbxrefs = parts[COLUMNS['DBXREFS']]
    assert(entrezId not in geneInfo)
    if 'ENSG' not in dbxrefs:
      continue
    dbxparts = dbxrefs.split('|')
    for dbxpart in dbxparts:
      if 'ENSG' in dbxpart:
        ensId = dbxpart.split(':')[1]
        assert('ENSG' in ensId)
    geneInfo[entrezId] = ensId
  return geneInfo

def readFromAdditonalGeneInfo(fileName, geneInfo, geneInfoPerChr):
  if 'chrM' not in geneInfoPerChr:
    geneInfoPerChr['chrM'] = {}
  for line in open(fileName):
    line = line.replace('\n', '')
    parts = line.split('\t')
    geneId = parts[0]    
    geneName = parts[1].upper()
    chrN = parts[2]
    txStart = parts[3]
    txEnd = parts[4]
    assert(geneName not in geneInfo)
    geneInfo[geneName] = [chrN, txStart, txEnd]
    geneInfoPerChr[chrN][geneName] = [chrN, txStart, txEnd]

def readFromAdditonalGeneInfoSeparated(fileName, geneInfo, geneInfoPerChr):
  if 'chrM' not in geneInfoPerChr:
    geneInfoPerChr['chrM'] = {}
  for line in open(fileName):
    line = line.replace('\n', '')
    parts = line.split('\t')
    geneId = parts[0]    
    geneName = parts[1].upper()
    chrN = parts[2]
    txStart = parts[3]
    txEnd = parts[4]
    assert(geneName not in geneInfo)
    geneInfo[geneName] = [[chrN, txStart, txEnd]]
    geneInfoPerChr[chrN][geneName] = [[chrN, txStart, txEnd]]


def readSnpAnnotationsFromVar(fileName):
  COLUMNS = {'CHR': 2, 'POS': 3, 'SNP': 7}
  snps = {}   
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    snp = parts[COLUMNS['SNP']]
    pos = int(parts[COLUMNS['POS']])
    chrNum = parts[COLUMNS['CHR']]
    snps[snp] = [chrNum, pos]
  return snps


def readSnpAnnotationsFromCD(fileName):
  COLUMNS = {'CHR': 1, 'POS': 2, 'SNP': 0}
  snps = {}
  #i = 0   
  for line in open(fileName, 'r'):
   # print i
    
    parts = line.replace('\n', '').split('\t')
    snp = parts[COLUMNS['SNP']]
    if len(parts[2]) == 0:
      print line
    pos = int(parts[COLUMNS['POS']])
    chrNum = 'chr' + parts[COLUMNS['CHR']]
    snps[snp] = [chrNum, pos]
    #i = i+1
  return snps


def readSnpAnnotationsFromVarPlinkFormat(fileName):
  COLUMNS = {'CHR': 2, 'POS': 3, 'SNP': 7}
  snps = {}   
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    snp = parts[COLUMNS['SNP']]
    pos = int(parts[COLUMNS['POS']])
    chrNum = parts[COLUMNS['CHR']]
    chrNum = chrNum[3:]
    if chrNum == 'M':
      chrNum = 'MT'
    snps[snp] = [chrNum, pos]
  return snps

def readSnpAnnotationsFromVarComplete(fileName):
  COLUMNS = {'TYPE': 0, 'GENE': 1, 'CHR': 2, 'POS': 3, 'SNP': 7}
  snps = {}   
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    snp = parts[COLUMNS['SNP']]
    pos = int(parts[COLUMNS['POS']])
    chrNum = parts[COLUMNS['CHR']]
    snpType = parts[COLUMNS['TYPE']]
    geneName = parts[COLUMNS['GENE']]
    snps[snp] = [chrNum, pos, snpType, geneName]
  return snps

def readEntToEns(fileName = '../GeneIds/ens_to_ent2.txt'):
  entToEnsMap = {}
  for line in open(fileName, 'r'):
    line = line.upper()
    parts = line.replace('\n', '').split(',')
    if len(parts) < 2 or len(parts[1]) == 0 or '_' in parts[0]:
      continue
    ensId = parts[0]
    entId = parts[1]
    if entId not in entToEnsMap:
      entToEnsMap[entId] = [ensId]
    else:
      entToEnsMap[entId].append(ensId)
  return entToEnsMap

def readEnsToEnt(fileName):
  ensToEntMap = {}
  for line in open(fileName, 'r'):
    line = line.upper()
    parts = line.replace('\n', '').split(',')
    if len(parts) < 2 or len(parts[1]) == 0 or '_' in parts[0]:
      continue
    ensId = parts[0]
    entId = parts[1]
    if ensId not in ensToEntMap:
      ensToEntMap[ensId] = [entId]
    else:
      ensToEntMap[ensId].append(entId)
  return ensToEntMap


def readCombinedNetworkEns(fileName = '../Pathways/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt'):
  networkMap = {}
  reverseNetworkMap = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    if parts[0] == 'Gene_A':
      continue
    networkMap['%s_%s' % (parts[0], parts[1])] = float(parts[2])
    reverseNetworkMap[('%s_%s' % (parts[1], parts[0]))] = float(parts[2])
  return networkMap, reverseNetworkMap


def readFromPathwayFile(fileName = '../Pathways/Human_GO_AllPathways_no_GO_iea_January_14_2014_symbol.gmt'):
  pathwayIdMap = {}
  pathwayIdToName = {}
  geneToPathwayIds = {}
  for line in open(fileName, 'r'):
    line = line.replace('\n', '')
    parts = line.split('\t')
    while '' in parts:
      parts.remove('')
    for i in range(len(parts)):
      parts[i] = parts[i]
    parts[0] = parts[0].upper()
    assert(parts[0] not in pathwayIdMap)
    if len(parts) <= 2:
      continue
    pathwayIdMap[parts[0]] = parts[2:]
    for p in parts[2:]:
      pathwayList = []
      if p in geneToPathwayIds:
        pathwayList = geneToPathwayIds[p]
      pathwayList.append(parts[0])
      geneToPathwayIds[p] = pathwayList
    pathwayIdToName[parts[0]] = parts[1]
  return pathwayIdMap, pathwayIdToName, geneToPathwayIds


def readFromPathwayCentralityFile(fileName = ''):
  pathwayCentrality = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    pathwayId = parts[0]
    geneCentralityMap = {}
    geneCentralities = parts[1].split(',')
    for geneCentralityInfo in geneCentralities:
      if len(geneCentralityInfo) == 0:
        continue
      gParts = geneCentralityInfo.split(':')
      geneId = gParts[0]
      centralityScore = float(gParts[1])
      geneCentralityMap[geneId] = centralityScore
    pathwayCentrality[pathwayId] = geneCentralityMap
  return pathwayCentrality


def readFromGenePIterationFile(fileName):
  genePsAllIterations = {}
  for line in open(fileName):
    if 'ENSG' not in line:
      continue
    parts = line.replace('\n', '').split('\t')
    assert('ENSG' in parts[0])
    fParts = []
    for p in parts[1:]:
      fParts.append(float(p))
    genePsAllIterations[parts[0]] = fParts
  return genePsAllIterations


def readFromGenePWithSnpInfo(fileName):
  genePs = {}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    geneId = parts[0]
    pVal = float(parts[1])
    genePs[geneId] = pVal
  return genePs

def getGeneNameToIdMapFromRef(fileName = 'refdata/all_gene_names_and_ids.txt', includeOtherNames = True):
  geneNameToIdMap = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    officialGeneName = parts[1].split('~')[0]
    otherGeneNames = parts[4].replace(' ', '').split(',')
    synGeneNames = parts[5].replace(' ', '').split(',')
    entrezId = parts[8]
    allGeneNames = [officialGeneName]
    if includeOtherNames:
      allGeneNames.extend(otherGeneNames)
      allGeneNames.extend(synGeneNames)
    while '' in allGeneNames:
      allGeneNames.remove('')
    for geneName in allGeneNames:
      geneName = geneName.upper()
      ids = []
      if geneName in geneNameToIdMap:
        ids = geneNameToIdMap[geneName]
      if len(entrezId) == 0:
        continue
      ids.append(entrezId)
      geneNameToIdMap[geneName] = ids
  return geneNameToIdMap


def getGeneIdToNameMap(geneNameToIdMap):
  geneIdToNameMap = {}
  for gname, gIds in geneNameToIdMap.iteritems():
    for gId in gIds:
      if gId not in geneIdToNameMap:
        geneIdToNameMap[gId] = []
      geneIdToNameMap[gId].append(gname)
  return geneIdToNameMap


def getGeneNameToIdMapFromExtra(fileName = 'refdata/extra_ids_2.txt', geneNameToIdMap = {}, replace = False):
  for line in open(fileName, 'r'):
    line = line.upper()
    parts = line.replace('\n', '').split('\t')
    if len(parts) == 1:
      print line
    geneName = parts[0]
    entrezId = parts[1]
    if len(entrezId) == 0:
      continue
    #print line
    ids = []
    if not replace:
      if geneName in geneNameToIdMap:
        ids = geneNameToIdMap[geneName]
      ids.append(entrezId)
      geneNameToIdMap[geneName] = ids
    else:
      geneNameToIdMap[geneName] = [entrezId]
  return geneNameToIdMap


def getEnsToEntMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  entToEnsMap = getEntToEnsMapGeneMania(fileName)
  ensToEntMap = {}
  for ent, ensIds in entToEnsMap.iteritems():
    for ensId in ensIds:
      if ensId not in ensToEntMap:
        ensToEntMap[ensId] = []
      ensToEntMap[ensId].append(ent)
  return ensToEntMap


def getEnsToSymbolMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  ensToSymMap = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    if 'Gene Name' not in parts[2]:
      continue
    #if parts[0] not in ensToEntMap:
    #  ensToEntMap[parts[0]] = []
    assert(parts[0] not in ensToSymMap)
    ensToSymMap[parts[0]] = parts[1]
  return ensToSymMap


def getEnsToEnsMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  ensToEnsMap = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    if 'Ensembl Gene ID' not in parts[2]:
      continue
    ensId1 = parts[0]
    ensId2 = parts[1]
    if ensId1 not in ensToEnsMap:
      ensToEnsMap[ensId1] = []
    ensToEnsMap[ensId1].append(ensId2)
  return ensToEnsMap


def getEntToEnsMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  ensToEnsMap = getEnsToEnsMapGeneMania(fileName)
  entToEnsMap = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    if 'Entrez' not in parts[2]:
      continue
    ensId = parts[0]
    entId = parts[1]
    if entId not in entToEnsMap:
      entToEnsMap[entId] = []
    entToEnsMap[entId].append(ensId)
    for e in ensToEnsMap[ensId]:
      if e not in entToEnsMap[entId]:
        entToEnsMap[entId].append(e)
  return entToEnsMap


def getEntToSymbolMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  entToSymMap = {}
  ensToSymMap = getEnsToSymbolMapGeneMania(fileName)
  ensToEntMap = getEnsToEntMapGeneMania(fileName)
  for e in ensToEntMap:
    if e not in ensToSymMap:
      continue
    for entId in ensToEntMap[e]:
      assert(entId not in entToSymMap)
      entToSymMap[entId] = ensToSymMap[e]
  return entToSymMap


def getSymbolToEntMapGeneMania(fileName = 'refdata/GeneMania/identifier_mappings.txt'):
  entToSymMap = getEntToSymbolMapGeneMania(fileName)
  symToEntMap = {}
  for e, s in entToSymMap.iteritems():
    if s not in symToEntMap:
      symToEntMap[s] = []
    #assert(s not in symToEntMap)
    symToEntMap[s].append(e)
  return symToEntMap

def readFromiCOGS(fileName):
  snpInfo = {}
  COLUMNS = {'SNP': 0, 'CHI2': 8, 'CHI2ERp': 13, 'CHI2ERn': 18}
  for line in open(fileName):
    if line[7] == 'C': # header
      continue
    line = line.replace('\n', '')
    parts = line.split(',')
    snpId = parts[COLUMNS['SNP']].replace('"', '')
    chi2 = parts[COLUMNS['CHI2ERn']]
    if chi2 == 'NA':
      chi2 = 0 #continue
    snpInfo[snpId] = float(chi2)
  return snpInfo

def readFromiCOGSAllele(fileName):
  snpInfo = {}
  COLUMNS = {'SNP': 0, 'ALLELE': 3}
  for line in open(fileName):
    if line[1] == 'S':
      continue
    line = line.replace('\n', '')
    parts = line.split(',')
    snpId = parts[COLUMNS['SNP']].replace('"', '')
    allele = parts[COLUMNS['ALLELE']].replace('"', '')
    if chi2 == 'NA':
      continue
    snpInfo[snpId] = float(chi2)
  return snpInfo

def readFromSnpsToGenesMapGenGen(fileName, snpFilters = None, DIST = 500000):
  snpsToGenes = {}
  geneToSnps = {}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    snpId = parts[0]
    if snpFilters is not None and snpId not in snpFilters:
      continue
    geneIds = parts[1].split(',')
    dist = parts[2]
    if int(dist) > DIST:
      continue
    snpsToGenes[snpId] = [geneIds, dist]
    for g in geneIds:
      if g not in geneToSnps:
        geneToSnps[g] = []
      geneToSnps[g].append(snpId)
  #print '-----------------'
  return snpsToGenes, geneToSnps


def readFromGenePsFile(fileName):
  geneDetails = {}
  rank = 0
  for line in open(fileName):
    rank += 1
    line = line.replace('\n', '')
    parts = line.split('\t')
    geneName = parts[0]
    geneP = float(parts[1])
    geneDetails[geneName] = [rank, geneP, parts[2]]
  return geneDetails

def readFromMyGseaResults(fileName):
  pathwayToDetails = {}
  rank = 0
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    size = int(parts[1])
    es = float(parts[2])
    nes = float(parts[3])
    pVal = 1.0 - st.norm.cdf(nes)
    fdr = parts[4]
    pathwayToDetails[pId] = [pVal, rank, nes, es, size, fdr]
    rank += 1
  return pathwayToDetails

def readFromMyGseaResultsCombined(fileName):
  pathwayToDetails = {}
  rank = 0
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    size = int(parts[1])
    nes = float(parts[2])
    fdr = parts[3]
    pVal = float(parts[4])

    pathwayToDetails[pId] = [pVal, rank, nes, 0, size, fdr]
    rank += 1
  return pathwayToDetails

def readFromMyGseaResultsCytoscape(fileName):
  pathwayToDetails = {}
  rank = 0
  for line in open(fileName):
    if line[0] == '#':
      continue
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    pName = parts[1]
    p = float(parts[2])
    q = parts[3]
    if q == 'NA':
      q = 1.0
    else:
      q = float(q)
    pathwayToDetails[pId] = [p, rank, q]
    rank += 1
  return pathwayToDetails

def readFromMyGseaResultsOnlyPathways(fileName):
  pathwayToDetails = {}
  rank = 0
  for line in open(fileName):
    if line[0] == '#':
      continue
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    try:
      p = float(parts[1])
    except:
      p = 1.0
    pathwayToDetails[pId] = [p, rank]
    rank += 1
  return pathwayToDetails


def readFromPathwayCommons5(fileName, interactions):
  geneInteractions = {}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    m1 = parts[0]
    m2 = parts[2]
    if not m1.isupper() or not m2.isupper():
      continue
    iType = parts[1]
    if len(interactions) > 0 and iType not in interactions:
      continue
    m12 = '%s&%s' % (m1, m2)
    if m12 not in geneInteractions:
      geneInteractions[m12] = []
    geneInteractions[m12].append(iType)
    if iType == 'neighbor-of' or iType == 'interacts-with' or iType == 'in-complex-with':
      m21 = '%s&%s' % (m2, m1)
      if m21 not in geneInteractions:
        geneInteractions[m21] = []
      geneInteractions[m21].append(iType)
  return geneInteractions

def readFromMAFDet(fileName):
  snpInfo = {}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    snpInfo[parts[0]] = [float(parts[1]), float(parts[2]), parts[3], parts[4]]
  return snpInfo
    

def readFromMyKEGG(fileName):
  entToSymMap = getEntToSymbolMapGeneMania()
  geneNameToIdMapFull = getGeneNameToIdMapFromRef(includeOtherNames = False)
  getGeneNameToIdMapFromExtra(fileName = 'refdata/my_extra_ids.txt', geneNameToIdMap = geneNameToIdMapFull, replace = False)
  geneIdToNameMap = getGeneIdToNameMap(geneNameToIdMapFull)

  allGeneIdToSymMap = {}

  pIdToGraph = {}
  pId = ''
  g = nx.DiGraph()
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    if parts[0] == '*':
      if len(pId) > 0:
        pIdToGraph[pId] = g
      pId = ''
      g = nx.DiGraph()
      continue
    if pId == '':
      pId = parts[0]
      nodes = parts[1:]
      nodesSym = []
      for n in nodes:
        geneSym = ''
        if n in entToSymMap:
          geneSym = entToSymMap[n]
        elif n in geneIdToNameMap:
          geneSym = geneIdToNameMap[n][0]
        else:
          geneSym = 'LOC' + n
        allGeneIdToSymMap[n] = geneSym
        nodesSym.append(allGeneIdToSymMap[n])
      g.add_nodes_from(nodesSym)
      continue
    if len(parts) > 1 and len(parts[1]) > 0:
      g.add_edge(allGeneIdToSymMap[parts[0]], allGeneIdToSymMap[parts[1]])
  pIdToGraph[pId] = g 
  return pIdToGraph

def readFromMappedHg19(fileName):
  hgPositions = {}
  COLUMNS = {'SOURCE': 1, 'MAPPED': 2, 'SOURCEC': 3, 'MAPPEDC': 4, 'SOURCEP': 7, 'MAPPEDP': 12}
  for line in open(fileName):
    parts = line.replace('\n', '').split('\t')
    if not (parts[COLUMNS['SOURCE']] == '1' and parts[COLUMNS['MAPPED']] == '1'):
      continue
    sourceFull = parts[COLUMNS['SOURCEC']] + ':' + parts[COLUMNS['SOURCEP']]
    mappedFull = parts[COLUMNS['MAPPEDC']] + ':' + parts[COLUMNS['MAPPEDP']]
    hgPositions[sourceFull] = mappedFull
  return hgPositions

def readFromExcludedSnpsWTCCC(fileName):
  excl = {}
  for line in open(fileName):
    if line[0] == '#':
      continue
    parts = line.replace('\n', '').split('\t')
    snpId = parts[1]
    excl[snpId] = 1
  return excl


