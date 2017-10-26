#!/usr/bin/python
import random
import numpy as np
from sklearn import svm





#OVERALL NOTES--created during code development so some may be irrelevant#

###I used one-class support vector machines to get cutoffs for CNVs based on the size (in bins) and integer distance of a segment
	##building two separate models
		#'Small' segments, using the immune loci recombination events
		#'Large' segments, using the sex chromosome and euploid segments
	##to compensate for the limited sample size of the immune recombination and potential skewing by outliers I am using a bootstrapping approach
		#run a large number of iterations
		#each time randomly sample with replacement from appropriate data source to build a dataset half as large as the original one
		#build a one class SVM for the large and small test groups
		#save the information about which size/intdist values pass and fail
	##combine the pass/fail results from all iterations
		#This lets me find a sort of confidence interval -- which values pass in >/= 95% of iterations

	##The overall border had issues after doing this
		#SVM seems to require a sort of symmetry where some samples are excluded in every (2D) direction
		#this makes for biologically inaccurate cases such as a size/int dist of 10 and 0.1 passes but 10 and 0.01 do not
		#so I needed to draw my own boundary

	##drawing boundaries strategy
		#minimum sizes (small or large) are the smallest number in the relevant 95% CI
			#therefore the max size of small segments is one less than the minimum of the large
				#at 25k this is only slightly larger than the largest observed immune recombination 
		#find maximal int dists at each size in the relevant 95% CI
			#for small segments these are filtered down to only the pairs where the intdist is >= that at smaller sizes
			#for large segments these are filtered to start at the size with the maximal intdist and then the first size/intdist pair where no larger size has a higher intdist
				#ignores that the distribution briefly increases at the smaller end sizes 
				#removes local dips like 0.2->0.19->0.2
				#taking the 'underside' of the curve, slightly more conserverative
		#the maximal size for large segments is the largest observed segment (so the size of chromosome 1)
		#the int dist at the max size is calculated by extending the slope of the nearest datapoints
		#the int dist at the small-large crossing point is the maximum of either
			#the cutoff for large if extended to the switching point using the slope of the nearest datapoints
			#the maximum passing intdist of small segments
	##cutoffs are saved as txt files for consistant future usage
  
###THIS CODE DOES PROVIDE SLIGHTLY DIFFERENT OUTPUTS UPON EACH USAGE DUE TO THE RANDOMIZATION COMPONENT
  ##10,000 iterations gave the most minimal variability for processing time required
  ##the variability that does occur does not appear to alter biological analysis










#perform one iteration of SVM on the input (bootstrapped) dataset#
def runIteration(iterationDict, parameterDict, xxx, yyy):
  Zdict = {}
  for i in iterationDict:
    thisModel = svm.OneClassSVM(kernel=parameterDict['kernel'], gamma=parameterDict['gamma'], nu=parameterDict['nu'])
    thisModel.fit(iterationDict[i])

    thisZ = thisModel.predict( np.c_[xxx.ravel(), yyy.ravel()] )
    Zdict[i] = thisZ.reshape(xxx.shape)

  return Zdict










#used when determining which combinations of size and integer distance are at the boundary between good and bad#
def CheckIfEdge(i1, i2, Zdata):
  edgeOptions =	[
      (i1, i2-1),
      (i1, i2+1),
      (i1-1, i2-1), 
      (i1-1, i2), 
      (i1-1, i2+1), 
      (i1+1, i2-1), 
      (i1+1, i2), 
      (i1+1, i2+1)
      ]
  edges = [x for x in edgeOptions if x[0] >= 0 and x[1] >= 0]

  thisValue = Zdata[i1, i2]
  surroundingData = [Zdata[x[0], x[1]] for x in edges]
  if any(x for x in surroundingData if x < thisValue):
    return True
  else:
    return False









#the master code that performs SVM and calculates the cutoff boundary#
  #dataDict is a dictionary of structured numpy arrays 
    #one array for each category of segments
  #rawOutputFilename is a path and name of a file 
    #include if you would like to save the SVM Frequency output to
  #same concept for smallCutoffFile and largeCutoffFFile
    #need paths if you want to save the information to a text file
    
def getCNVcutoffs(dataDict, rawOutputFilename=False, smallCutoffFile=False, largeCutoffFile=False):

	###variable setup
  
	variables = ['bins', 'intdist']
	parameterDict = {'kernel': 'rbf', 'gamma': 10, 'nu': 0.125}
  iterations = 10000
 # numBins = 25
  
	maxDict = {x: max([max(dataDict[y][x]) for y in dataDict]) for x in variables}
	intdistVals = np.arange(0.00, 0.51, 0.01)
	sizeVals = np.array([0.9] + list(np.arange(1, 10001, 1)))
	xx, yy = np.meshgrid( np.log10(sizeVals), intdistVals )
	xxx = xx/np.log10(maxDict['bins'])
	yyy = yy/maxDict['intdist']

	testDataDict =	{
			'Immune': [[np.log10(x[y])/np.log10(maxDict[y]) if y == 'bins' else x[y]/maxDict[y] for y in variables] for x in dataDict['Immune']],
			'Sex Chrom': [[np.log10(x[y])/np.log10(maxDict[y]) if y == 'bins' else x[y]/maxDict[y] for y in variables] for x in dataDict['Sex Chrom']],
			'Euploid': [[np.log10(x[y])/np.log10(maxDict[y]) if y == 'bins' else x[y]/maxDict[y] for y in variables] for x in dataDict['Euploid']],
			'Putative CNV': [[np.log10(x[y])/np.log10(maxDict[y]) if y == 'bins' else x[y]/maxDict[y] for y in variables] for x in dataDict['Putative CNV']]
			}

	downsampleSizeDict = { x: len(testDataDict[x]) / 2 for x in testDataDict }

	maxSizeEup = np.log10(max(dataDict['Euploid']['bins']))





	###Run bootstrapping of fitting a one-class SVM model to large and small alterations	
	allZdict = {'small': np.zeros(xxx.shape), 'large': np.zeros(xxx.shape)}

  #note: doing this in multiple iterations to avoid running out of memory
	for i in range(iterations/250):

		argList = []
		for j in range(250):
			iterationData = { x: [random.choice(testDataDict[x]) for y in range(downsampleSizeDict[x])] for x in testDataDict if x != 'Putative CNV'}
	
			thisDict =	{
					'small': np.array(iterationData['Immune']),
					'large': np.array( iterationData['Sex Chrom'] + iterationData['Euploid'] )
					}
			argList.append( (thisDict, parameterDict, xxx, yyy, ) )


		startArg = i * 250
		endArg = min([iterations, (i+1)*250])
		print i+1, 'round of iteration multiprocessing, will go from', startArg, 'to', endArg-1
#!#!#!#!#!#!#!#!#!#!#!#!#!#!# I NEED TO FIX THIS DAEMON THING
		results = daemon( runIteration, argList, 0.5, True, 'SVM bootstrapping round ' + str(i+1) )

		for j in results:
			allZdict = { x: allZdict[x] + j[x] for x in allZdict }
		print '\n'



	allZdictNorm = { x: (allZdict[x] / float(iterations)) for x in allZdict }

	thisShape = allZdictNorm[allZdictNorm.keys()[0]].shape
	overallZ = np.zeros(thisShape)
	for i in range(thisShape[0]):
		for j in range(thisShape[1]):
			thisMax = max([ allZdictNorm[x][i,j] for x in allZdictNorm ])
			overallZ[i,j] = thisMax

  if rawOutputFilename:
	  np.savetxt(rawOutputFilename, overallZ, delimiter='\t')





	###Interpret bootstrapping results
	cutoffDtypes = {'names': ('bins', 'intdist'), 'formats': ('int', 'float')}

	Zdict95 = {x: [index for index,value in np.ndenumerate(allZdictNorm[x]) if value >= 0.95] for x in allZdictNorm}
	Zdict95edge = {x: [y for y in Zdict95[x] if CheckIfEdge(y[0], y[1], allZdictNorm[x])] for x in Zdict95}
	Zdict95edge = { x: np.sort( np.array([ (10**xx[y[0], y[1]], yy[y[0], y[1]]) for y in Zdict95edge[x] ], dtype=cutoffDtypes), order='bins' ) for x in Zdict95 }

	maxSizePassingLym = min(Zdict95edge['large']['bins']) - 1



	filteredZedgeDict = {}
	for i in Zdict95edge:
		filterData = []
		thisData = Zdict95edge[i]
		for j in np.unique(thisData['bins']):
			useTest = True

			sizeData = thisData[thisData['bins']==j]
			sizeData = np.sort(sizeData, order='intdist')

			if i == 'small':
				if len(filterData) != 0 and sizeData['intdist'][-1] <= filterData[-1][1]:
					useTest = False
			else:
				if (len(filterData) == 0 and max(sizeData['intdist']) < max(thisData['intdist'])):
					useTest = False
				elif len(filterData) != 0 and np.round(max(sizeData['intdist']), 3) == np.round(filterData[-1][1], 3):
					useTest = False
				elif j < maxSizePassingLym:
					useTest = False

				elif len(filterData) > 0 and max(sizeData['intdist']) > filterData[-1][1]:
					useTest = False
					while max(sizeData['intdist']) > filterData[-1][1]:
						filterData.pop()

				if len(sizeData) > 0 and max(sizeData['intdist']) != sizeData['intdist'][-1]:
					print 'ERROR'


			if useTest:
				filterData.append( (j, sizeData['intdist'][-1]) )



		if i == 'small':
			filterData.insert(0, (filterData[0][0], 0.) )
			filterData.append( (maxSizePassingLym, filterData[-1][1]) )
			filterData.append( (maxSizePassingLym, 0.) )

		else:
			smallDiff = filterData[0][0] - maxSizePassingLym
			diffDown = filterData[0][0] + smallDiff
			diffEdge = [(x, abs(y[0] - diffDown)) for x,y in enumerate(filterData)]
			diffMin = min([x[1] for x in diffEdge])
			diffEdge = [filterData[x[0]] for x in diffEdge if x[1] == diffMin]
			nextIntDist = max([x[1] for x in diffEdge])
			diffIntDist = filterData[0][1] - nextIntDist
			filterDiff = filterData[0][1] + diffIntDist
			if filterDiff < filteredZedgeDict['small'][-2][1]:
				filterDiff = filteredZedgeDict['small'][-2][1]
			filterData.insert( 0, (filteredZedgeDict['small'][-2][0], filterDiff) )



			bigDiff = 10**maxSizeEup - filterData[-1][0]
			diffUp = filterData[-1][0] - bigDiff
			diffEdge = [(x, abs(y[0] - diffUp)) for x,y in enumerate(filterData)]
			diffMin = min([x[1] for x in diffEdge])
			diffEdge = [filterData[x[0]] for x in diffEdge if x[1] == diffMin]

			prevIntDist = max([x[1] for x in diffEdge])
			diffIntDist = prevIntDist - filterData[-1][1]
			filterDiff = filterData[-1][1] - diffIntDist

			filterData.append( ((10 ** maxSizeEup), filterDiff) )
			filterData.append( ((10 ** maxSizeEup), 0.) )

		filterData = np.array(filterData, dtype=cutoffDtypes)
		filteredZedgeDict[i] = filterData



	if max(filteredZedgeDict['large']['intdist']) > max(filteredZedgeDict['small']['intdist']):
		filteredZedgeDict['small']['intdist'][-2] = max(filteredZedgeDict['large']['intdist'])





	







	if smallCutoffFile:
    OUT = open(smallCutoffFile, 'w')
    OUT.write('Size\tIntegerDistance\n')

    smallSizeOptions = np.arange(filteredZedgeDict['small']['bins'][0], filteredZedgeDict['small']['bins'][-1]+1)
    smallThresholdDict = {}
    smallData = filteredZedgeDict['small'][1:-1]
    for i in smallSizeOptions:
      if i in smallData['bins']:
        cutoff = [x['intdist'] for x in smallData if x['bins'] == i]
        if len(cutoff) > 1:
          print 'ERROR: why are there multiple cutoffs? Fix edge assignment code', cutoff
        cutoff = cutoff[0]

      else:
        distances = [(x, i-y) for x,y in enumerate(smallData['bins'])]
        bigger = [x for x in distances if x[1] < 0][0][0]
        smaller = [x for x in distances if x[1] > 0][-1][0]

        big = smallData[bigger]
        small = smallData[smaller]

        slope = (big['intdist'] - small['intdist']) / float(big['bins'] - small['bins'])
        cutoff = ( slope * (i - small['bins']) ) + small['intdist']

      smallThresholdDict[int(i)] = np.round(cutoff, 2)
      OUT.write(str(int(i)))
      OUT.write('\t')
      OUT.write(str(np.round(cutoff, 2)))
      OUT.write('\n')

    OUT.close()





	if largeCutoffFile:
    OUT = open(largeCutoffFile, 'w')
    OUT.write('Size\tIntegerDistance\n')

    largeSizeOptions = np.arange(smallSizeOptions[-1]+1, filteredZedgeDict['large']['bins'][-1]+1)
    largeThresholdDict = {}
    largeData = filteredZedgeDict['large'][:-1]
    for i in largeSizeOptions:
      if i in largeData['bins']:
        cutoff = [x['intdist'] for x in largeData if x['bins'] == i]
        if len(cutoff) > 1:
          print 'ERROR:', cutoff
        cutoff = cutoff[0]

      else:
        distances = [(x, i-y) for x,y in enumerate(largeData['bins'])]
        bigger = [x for x in distances if x[1] < 0][0][0]
        smaller = [x for x in distances if x[1] > 0][-1][0]
        big = largeData[bigger]
        small = largeData[smaller]

        slope = (big['intdist'] - small['intdist']) / float(big['bins'] - small['bins'])
        cutoff = ( slope * (i - small['bins']) ) + small['intdist']

      largeThresholdDict[int(i)] = np.round(cutoff, 2)
      OUT.write(str(int(i)))
      OUT.write('\t')
      OUT.write(str(np.round(cutoff, 2)))
      OUT.write('\n')

    OUT.close()
    
    
    
