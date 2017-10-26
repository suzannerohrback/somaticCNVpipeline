#!/usr/bin/python
import numpy as np





#OVERALL NOTES - created during code development#

###DETERMINING QC CUTOFFS###	
##noticed strong clustering behavior of samples with good QC values (high reads, low noise, high CS) when I made a 3D plot, wanted to use this to define the cutoff boundaries
##so I wrote a sort of contraction algorithm to define which samples are in this cluster
	#range-scale all statistics so the minimum value is 0 and the maximum is 1 (so that everything has equal weight, which otherwise gets screwed up by the reads)
	#calculate a distance matrix, using euclidian distance of all three QC variables
	#find the minimum pairwise distance among all comparisons and connect those points in a cluster
	#to shorten/simplify processing I am not looking for multiple connections within a cluster, once there is a path between 2 samples, their direct distance measurement is removed from the pool of consideration
		#ex, round 1: A connects to B, round 2: B connects to C, round 3: A to C is smaller than C to D but C to D expands the cluster and is selected
	#once 95% of samples are connected, iteration stops

	#Notes: 
		#I'm not actually paying attention to whether or not the samples were identified because the clustering is so clear and showed a natural bias against the unidentified lymphocytes
		#when using a 100% sample clustering cutoff outliers start to be drawn into the good cluster, so I cut off at 95%
		#expanding intervals slightly to account for sample size limitations (disproportionate # of WGA4 samples get cut out if I go exact)
			#reads-round down to nearest 100000 for minimum (622705 -> 600000)
			#MAPD-round up to nearest 0.5 multiple (0.38 -> 0.4)
			#CS-round down to nearest 0.5 multiple (0.81 -> 0.80)
			#actually only includes an extra 18 samples but I feel more comfortable with being able to generalize with these values










#calculate euclidean distance between two datapoints#
def getDist(data1, data2):
	thisSum = 0
	for i in ['Reads', 'MAPD', 'CS']:
		thisSum += (data1[i] - data2[i]) ** 2
	eucDist = thisSum ** 0.5
	return eucDist










#determine appropriate QC cutoffs#
  #data is a structured array (entries following the .txt file in this repo)
  
def qcCompare(data, name, numBins):

  #set up variables
  testData = np.copy(data)
	currentDtype = testData.dtype.descr
	readLoc = [x for x,y in enumerate(currentDtype) if y[0] == 'Reads'][0]
	readEntryFix = ('Reads', 'float')
	currentDtype[readLoc] = readEntryFix
	fixDtype = np.dtype(currentDtype)
	testData = testData.astype(fixDtype)





  #range scale data and calculate distance metric
  for i in ['Reads', 'MAPD', 'CS']:
    testData[i] = ( testData[i] - min(testData[i]) ) / ( max(testData[i]) - min(testData[i]) )

	allDistDict = {x['Sample']: {y['Sample']: getDist(x, y) for y in testData if y['Sample'] != x['Sample']} for x in testData}





  #iterate through contraction algorithm to merge similar data points
	connectionDict = {}
	connectedPoints = []
	roundCount = 0
	while len(connectedPoints) <= len(allDistDict) * 0.95:
		minPairs = {x: [(y, allDistDict[x][y]) for y in allDistDict[x] if allDistDict[x][y] == min(allDistDict[x].values())][0] for x in allDistDict}
		nameDists = {x: minPairs[x][1] for x in minPairs}
        
		overallMin = [(x, minPairs[x][0]) for x in minPairs if nameDists[x] == min(nameDists.values())][0]
		minData1 = data[data['Sample'] == overallMin[0]]
		minData2 = data[data['Sample'] == overallMin[1]]
                
		if overallMin[0] in connectionDict:
			connectionDict[overallMin[0]].append(overallMin[1])
		else:
			connectionDict[overallMin[0]] = [overallMin[1]]
            
		if overallMin[1] in connectionDict:
			connectionDict[overallMin[1]].append(overallMin[0])
		else:
			connectionDict[overallMin[1]] = [overallMin[0]]
            
		connectedPoints = list(set(sum([connectionDict[x] for x in connectionDict], [])))
        
		for i in connectionDict:
			for j in connectionDict[i]:
				otherConnections = [i] + [x for x in connectionDict[i] if x != j]
				connectionDict[j] = list(set( connectionDict[j] + otherConnections ))
                
				if i in allDistDict.keys() and j in allDistDict[i].keys():
					del allDistDict[i][j]

		roundCount += 1
		if roundCount % 10 == 0:
			print '\t', roundCount, 'iterations complete', len(connectedPoints), 'connected points out of', len(allDistDict)
		if roundCount >= 500:
			print 'ERROR: why are so many contraction iterations needed?'
			raise SystemExit





  #Interpret clusters to obtain cutoff values
	clusterSizes = list(set([len(connectionDict[x]) for x in connectionDict]))
	goodCluster = [[x] + connectionDict[x] for x in connectionDict if len(connectionDict[x]) == max(clusterSizes)][0]
	print 'The largest cluster contains', len(goodCluster), 'samples',
	print sum([1 for x in data if x['Sample'] in goodCluster and x['Type'] in ['B', 'T']]), 'of', len(IDdata), 'are identified',
	print sum([1 for x in data if x['Sample'] in goodCluster and x['Type'] not in ['B', 'T']]), 'of', len(UNdata), 'are not identified'

	minReads = min([data[data['Sample'] == x]['Reads'][0] for x in goodCluster])
	maxReads = max([data[data['Sample'] == x]['Reads'][0] for x in goodCluster])

	minMAPD = min([data[data['Sample'] == x]['MAPD'][0] for x in goodCluster])
	maxMAPD = max([data[data['Sample'] == x]['MAPD'][0] for x in goodCluster])

	minCS = min([data[data['Sample'] == x]['CS'][0] for x in goodCluster])
	maxCS = max([data[data['Sample'] == x]['CS'][0] for x in goodCluster])
  
  

	readsCutoff = np.arange(0, minReads, 100000)[-1]
	mapdCutoff = np.arange(0, maxMAPD+0.05, 0.05)[-1]
	csCutoff = np.arange(0, minCS, 0.05)[-1]

	cutoffDict = {'Reads': readsCutoff, 'MAPD': mapdCutoff, 'CS': csCutoff}
	print '\nOptimal QC cutoffs:\n', cutoffDict, '( from', minCS, minReads, maxMAPD, ')\n'
