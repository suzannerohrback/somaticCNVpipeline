#!/usr/bin/python
import sys
import os
import inspect
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg

normVars = cfg.Segment()
#import statsmodels
from statsmodels.nonparametric.smoothers_lowess import lowess as lowess









#Lowess regression to correct bin counts for gc-content bias#
def runLowess(counts, gc):
	counts = counts + 1
	counts = counts / np.median(counts)
	counts = np.log(counts)
	
#	lowessModel = statsmodels.nonparametric.smoothers_lowess.lowess(counts, gc, frac=0.05, return_sorted=False)
	lowessModel = lowess(counts, gc, frac=0.05, return_sorted=False)

	lowessData = counts - lowessModel
	lowessData = np.exp(lowessData)
	lowessData = np.log2(lowessData)
	
	return lowessData










#Adjust the sex chromosome method reference values to -> 2 (otherwise everything will appear to be 2X/2Y)#
def adjustSexChroms(data, xLocs, yLocs, numSamples):
	euploidMedian = np.median(data[:xLocs[0]])
	xMedian = np.median(data[xLocs])
	yMedian = np.median(data[yLocs])
	
	xCorrection = euploidMedian / xMedian
	yCorrection = euploidMedian / yMedian
	
	#every male should generate (euploid / 2) counts...so only one male would generate (euploid / (2 x numSamples) ) counts
	maleTest = True
	maleMin = euploidMedian / (2 * numSamples)
	if yMedian < maleMin / 2:	#Is this overly conservative? Need to test
		maleTest = False
		
	data[xLocs] = data[xLocs] * xCorrection
	if maleTest:
		data[yLocs] = data[yLocs] * yCorrection
		
	return data, maleTest





#Make a bincount reference to use for normalizing single cell data#
	#There must be at least 10 samples with at least 500,000 counted reads to run
def runMakeMethodRef(species, sampleList, methodName, lowessDir):
	if len(sampleList) < 10:
		return False
		

	
	binArray = common.importInfoFile(normVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)
	xLocs = [x for x,y in enumerate(binArray['chrom']) if y == 'chrX']
	yLocs = [x for x,y in enumerate(binArray['chrom']) if y == 'chrY']
	
	mergeArray = np.zeros(len(binArray), dtype='int')
	
	sampleCount = 0
	for i in sampleList:
		data = np.loadtxt(i, usecols=[3], dtype='int')
		
		if sum(data) < 500000:
			continue
			
		mergeArray = mergeArray + data
		sampleCount += 1
	
	if sampleCount < 10:
		return False
		
	
	
	mergeArray, maleTest = adjustSexChroms(mergeArray, xLocs, yLocs, sampleCount)
	
	lowessData = runLowess(mergeArray, binArray['gc'])
	
	if not maleTest:
		lowessData[yLocs] = len(yLocs) * [0.]

	np.savetxt(lowessDir + methodName + '.methodRef.lowess.txt', lowessData)
	
	
	
	printText = '\tAmplification method reference for ' + methodName + ' has been generated from ' + str(sampleCount) + ' samples '
	if maleTest:
		printText += 'which included at least one male\n'
	else:
		printText += 'with included no male samples\n'
	print(printText)

	return lowessData
	
	








#Run normalization on a single bincount sample#
def runNormalizeOne(species, infile, methodRef, outfile):
#	normVars = cfg.Segment()
	
	binArray = common.importInfoFile(normVars.binDict[species], [0, 2, 6], 'normref', skiprows=1)

	data = np.loadtxt(infile, usecols=[3], dtype='int')
	lowessData = runLowess(data, binArray['gc'])
	
	if methodRef:
		lowessData = lowessData - methodRef
		
	np.savetxt(outfile, lowessData)
	



