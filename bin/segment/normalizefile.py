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










#Lowess regression to correct bin counts for gc-content bias#
def runLowess(counts, gc, saveName):
	counts = counts + 1
	counts = counts / np.median(counts)
	counts = np.log(counts)
	
	lowessModel = statsmodels.nonparametric.smoothers_lowess.lowess(counts, gc, frac=0.05, return_sorted=False)

	lowessData = counts - lowessModel
	lowessData = np.exp(lowessData)
	lowessData = np.log2(lowessData)
	
	np.savetxt(saveName, lowessData)

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
	#Sample must have at least 500,000 counted reads to be used
def runMakeMethodRef(species, sampleList, methodName, lowessDir):
	if len(sampleList) < 10:
		return False
		
	normVars = cfg.Segment()
	
	binArray = common.importInfoFile(normVars.binDict[species], [0, 2, 6], 'normref', skiprows=1)
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
	
	mergeArray, maleTest = adjustSexChroms(mergeArray, xLocs, yLocs, sampleCount)
	
	lowessData = runLowess(mergeArray, binArray['gc'], lowessDir + methodName + '.methodRef.lowess.txt')
	
	if not maleTest:
		lowessData[yLocs] = len(yLocs) * [0.]

	return lowessData
	
	








#Run normalization on a single bincount sample#
def runNormalizeOne(species):
	normVars = cfg.Segment()
	
	binArray = common.importInfoFile(normVars.binDict[species], [0, 2, 6], 'normref', skiprows=1)

	return 0











