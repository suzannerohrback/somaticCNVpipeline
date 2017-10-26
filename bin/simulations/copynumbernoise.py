#!/usr/bin/python
import numpy as np
import scipy.stats as stats











#Variable explanation
  #I believe maleSamples and femaleSamples are meant to be lists of the required samples to use for processing
  #numBins is likely archaic to my system but included in case compatability is needed
  
def checkXcounts(maleSamples, femaleSamples, numBins=25):
	def basicCheck(sampleNames, expected, allCounts, gender):
		IQRs = []
		for j,i in enumerate(sampleNames):
      ##########Add the correct filepath
			#filename = bincountsFolder + i + '.' + str(numBins) + 'k.bincounts.txt'
			binCounts = np.loadtxt(filename, usecols=[3])

			if expected == 0.5:
				allCounts[j, :] = binCounts

			normCounts = np.log2( binCounts / np.median(binCounts) )
			normXbins = normCounts[xLoc]
			thisIQR = np.percentile(normXbins, 75) - np.percentile(normXbins, 25)
			IQRs.append(thisIQR)

		if expected == 0.5:
			return IQRs, allCounts
		else:
			return IQRs





##########Add the correct filepath for your system
#	chromList = list(np.loadtxt('filepath/example/' + cfg.species + '/' + cfg.species + '.varbin.gc.content.' + str(numBins) + 'k.bowtie.k36.txt', dtype='S5', usecols=[0]))
	xLoc = [x for x,y in enumerate(chromList) if y == 'chrX']
	autoLoc = range(xLoc[0])

	allCounts = np.zeros([len(maleSamples), len(chromList)])


	maleIQR, allCounts = basicCheck(maleSamples, 0.5, allCounts, 'male')
	femaleIQR = basicCheck(femaleSamples, 1, False, 'female')


	sums = [sum(allCounts[x,:]) for x in range(allCounts.shape[0])]
	rankData = map(int,list(stats.rankdata(sums,method='ordinal')))
	rankDict = {rankData[x]:x for x,y in enumerate(sums)}

	male2IQR = []
	for i in range(1,len(rankData)):
		testData = allCounts[rankDict[i], :]
		addData = allCounts[rankDict[i+1], :]
		testData[xLoc] = testData[xLoc] + addData[xLoc]

		normCounts = np.log2( testData / np.median(testData) )
		normXbins = normCounts[xLoc]
		thisIQR = np.percentile(normXbins, 75) - np.percentile(normXbins, 25)
		male2IQR.append(thisIQR)




  #print out a brief stats comparison report#
	data = [maleIQR, femaleIQR, male2IQR]
	medians = [np.median(maleIQR), np.median(femaleIQR), np.median(male2IQR)]
	lowIQR = [medians[0] - np.percentile(maleIQR, 25), medians[1] - np.percentile(femaleIQR, 25), medians[2] - np.percentile(male2IQR, 25)]
	highIQR = [np.percentile(maleIQR, 75) - medians[0], np.percentile(femaleIQR, 75) - medians[1], np.percentile(male2IQR, 75) - medians[2]]
	names = ['Male', 'Female', 'Male+Male']
	
	diffDict = {}

	for i,j in enumerate(medians):
		if names[i] != names[-1]:
			for k,l in enumerate(names[i+1:]):
				thisMWU = 2 * stats.mannwhitneyu(data[i], data[k+i+1])[1]
				diffDict[(i, k+i+1)] = thisMWU

	###correction to p-value for multiple comparisons being made###
	numDiff = sum([1. for x in diffDict if diffDict[x] < 0.05])
	diffDict = {x: numDiff*diffDict[x] for x in diffDict} 

	for i,j in enumerate(medians):
		print '\tAssessing:', names[i]
		print '\t\tMedian (IQR):', j, '(', np.percentile(data[i], 25) , '-', np.percentile(data[i], 75), ')', lowIQR[i], highIQR[i]
		if names[i] != names[-1]:
			for k,l in enumerate(names[i+1:]):
				print '\t\tMWU vs', l, ':', diffDict[(i, k+i+1)]

