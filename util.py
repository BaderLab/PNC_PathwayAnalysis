import operator
import random
import gzip

def printSortedMap(regularMap):
  sortedMap = sorted(regularMap.iteritems(), key=operator.itemgetter(1))
  for s in sortedMap:
    print s[0] + '\t' + str(s[1])


def printSortedMapToFile(regularMap, fileName):
  sortedMap = sorted(regularMap.iteritems(), key=operator.itemgetter(1))
  f = open(fileName, 'w')
  for s in sortedMap:
    f.write(s[0] + '\t' + str(s[1]) + '\n')
  f.close() 

def printMapToFile(regularMap, fileName):
  f = open(fileName, 'w')
  for k, v in regularMap.iteritems():
    f.write(k + '\t' + str(v) + '\n')
  f.close()

def printMapWithList(mapWithList):
  for key, lst in mapWithList.iteritems():
    print key + '\t' + '\t'.join(str(x) for x in lst)

def printMapWithListToFile(mapWithList, fileName):
  f = open(fileName, 'w')
  for key, lst in mapWithList.iteritems():
   # print key, len(lst)
    f.write(key + '\t' + '\t'.join(str(x) for x in lst) + '\n')
  f.close()

def print2DList(tdList):
  for r in tdList:
    cs = ''
    for c in r:
      cs += str(c) + '\t'
    print cs

def printSnpAlleleCounts(s):
  for snpId in s:
    cs = ''
    for j, v in s[snpId].iteritems():
      if v > 0:
        cs += j + ':' + str(v) + '\t'
    print snpId, cs

def readMapWithListFromFile(fileName):
  mapWithList = {}
  for line in open(fileName, 'r'):
    parts = line.replace('\n', '').split('\t')
    key = parts[0]
    elements = parts[1:]
    mapWithList[key] = elements
  return mapWithList

def readMapWithOneFloat(fileName, ignoreInvalid = True):
  mapWithFloat = {}
  for line in open(fileName, 'r'):
    if len(line) == 0:
      continue
    parts = line.replace('\n', '').split('\t')
    if len(parts) < 2:
      print line
      continue
    if (parts[1] == 'NA' or parts[1] == 'NaN') and ignoreInvalid:
      parts[1] = 0.0
      #continue
    mapWithFloat[parts[0].strip()] = float(parts[1])
  return mapWithFloat

def readMapWithOneFloatGzip(fileName, ignoreInvalid = True):
  mapWithFloat = {}
  for line in gzip.open(fileName, 'r'):
    if len(line) == 0:
      continue
    parts = line.replace('\n', '').split('\t')
    if len(parts) < 2:
      print line
      continue
    if (parts[1] == 'NA' or parts[1] == 'NaN') and ignoreInvalid:
      parts[1] = 0.0
      #continue
    mapWithFloat[parts[0].strip()] = float(parts[1])
  return mapWithFloat


def readListAsMapFromFile(fileName):
  lst = {} 
  for line in open(fileName, 'r'):
    lst[line.replace('\n', '')] = 1
  return lst

def combineCounts(allSnpInfos):
  combinedSnpInfo = {}
  for snpInfo in allSnpInfos:
    for snpId, counts in snpInfo.iteritems():
      if snpId not in combinedSnpInfo:
        combinedSnpInfo[snpId] = counts
      else:
        for c in counts:
          combinedSnpInfo[snpId][c] += counts[c] 
  return combinedSnpInfo

def invertMap(originalMapWithList):
  invertedMap = {}
  for k, lst in originalMapWithList.iteritems():
    for l in lst:
      if l not in invertedMap:
        invertedMap[l] = [k]
      else:
        invertedMap[l].append(k)
  return invertedMap

def sortAndReturnIndex(lst):
  newMap = {}
  for i in range(len(lst)):
    newMap[i] = lst[i]
  sortedMap = sorted(newMap.iteritems(), key=operator.itemgetter(1))
  sortedMapIndex = {}
  for s in range(len(sortedMap)):
    furthestIndex = s
    for s2 in range(s, len(sortedMap)):
      if sortedMap[s][1] == sortedMap[s2][1]:
        furthestIndex = s2
      else:
        break
    sortedMapIndex[sortedMap[s][0]] = furthestIndex 
  return sortedMapIndex

def sortAndReturnIndexBig(lst):
  newMap = {}
  for i in range(len(lst)):
    newMap[i] = lst[i]
  sortedMap = sorted(newMap.iteritems(), key=operator.itemgetter(1))
  sortedMapIndex = {}
  for s in range(len(sortedMap)):
    furthestIndex = s
    for s2 in reversed(range(0, s + 1)):
      if sortedMap[s][1] == sortedMap[s2][1]:
        furthestIndex = s2
      else:
        break
    sortedMapIndex[sortedMap[s][0]] = furthestIndex
  return sortedMapIndex

def getFriendlyPathwayName(pFull):
  parts = pFull.split('%')
  return (parts[0], parts[1], parts[2])

def getGroupNames(groupIds, expr):
  res = []
  eps = expr.split('.')
  for g in groupIds:
    gps = g.split('.')
    assert(len(gps) == len(eps))
    dontInclude = False
    for i in range(len(gps)):
      hasMatch = False
      if eps[i] == 'X':
        hasMatch = True
      epsp = eps[i].split(',')
      for epspi in epsp:
        if epspi == gps[i]:
          hasMatch = True
      if not hasMatch:
        dontInclude = True
        break
    if not dontInclude:
      res.append(g)
  return res

def shuffleValues(m):
  vs = m.values()
  random.shuffle(vs)
  i = 0
  rndM = {}
  for k in m:
    rndM[k] = vs[i]
    i += 1
  return rndM 

      

