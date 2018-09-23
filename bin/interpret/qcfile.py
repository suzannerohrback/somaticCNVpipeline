#!/usr/bin/python
import os, sys, inspect
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg










def calcReads(sample, statsDir):
	infile = statsDir + common.findInfile(sample, statsDir, ext='.bincount.stats.txt')
	
	with open(infile, 'r') as IN:
		data = IN.readline()
		readCount = data.rstrip().split('\t')[1]

	return int(readCount)

	
	
	
	
	
def calcMAPD(sample, lowessDir):
	infile = lowessDir + common.findInfile(sample, lowessDir)
	
	data = np.loadtxt(infile)
	MAPD = np.median(abs(data[1:] - data[:-1]))

	return MAPD, data
	
	
	
	
	
def calcCS(data):
	CS = 1 - ( 2 * np.median( abs( data - np.round(data) ) ) )
	return CS
	
	
	
	
	
def getPloidy(segData, sample, plotDir, cutoff):
	ploidyTestValues = np.arange(1.25, 2.76, 0.01)
	CSarray = np.zeros(len(ploidyTestValues))
	peakPloidy = 1.
	peakCS = 0.
	
	for i,j in enumerate(ploidyTestValues):
		testData = np.copy(segData) * j
		CSarray[i] = calcCS(testData)
		
		if CSarray[i] > peakCS:
			peakCS = CSarray[i]
			peakPloidy = j

			
	xTicks = np.arange(1.25, 2.76, 0.25)
	yTicks = np.arange(0, 1.1, 0.2)
	
	fig, ax = plt.subplots()		
	
	ax.plot(ploidyTestValues, CSarray, color='#2000b1', lw=3)
	ax.plot([2, 2], [0, 1.1], color='#6b7c85', lw=0.5, zorder=0)
	ax.plot([1., 3], [cutoff, cutoff], color='#6b7c85', lw=0.5, zorder=0)
	
	ax.set_xticks(xTicks)
	ax.set_xticklabels(xTicks)
	ax.set_xlabel('Ploidy Value', labelpad=5)
	ax.set_xlim(1.25, 2.75)

	ax.set_yticks(yTicks)
	ax.set_yticklabels(yTicks)
	ax.set_ylabel('Confidence Score', labelpad=5)
	ax.set_ylim(-0.02, 1.02)
	
	ax.tick_params(direction='out', which='both', pad=0., length=3, top='off', right='off')
	
	fig.set_size_inches(4, 4, forward=True)
	plt.subplots_adjust(left=0.15, right=0.98, bottom=0.12, top=0.95)
	plt.savefig(plotDir + sample + '.ploidyDeterminationPlot.png', dpi=666)
	plt.close()
	
	return peakCS, peakPloidy





def getGender(data, chroms, ploidy):
	binCN = 2 ** data * ploidy
	
	xData = [y for x,y in enumerate(binCN) if chroms[x] == 'chrX']
	xCN = np.median(xData)
	xGender = 'NA'
	if np.round(xCN) >= 2.:
		xGender = 'F'
	elif np.round(xCN) <= 1.:
		xGender = 'M'

	yData = [y for x,y in enumerate(binCN) if chroms[x] == 'chrY']
	yCN = np.median(yData)
	yGender = 'NA'
	if np.round(yCN) < 1.:
		yGender = 'F'
	elif np.round(yCN) >= 1.:
		yGender = 'M'

	
	
	if xGender == yGender and xGender != 'NA':
		gender = xGender

	elif (yCN >= 0.25 and xGender == 'M') or np.round(np.percentile(yData, 75)) == 1.:
			gender = 'M'
		
	###Note: this is strange, but gender can be indistinguishable in this case, 
		#so going the more conservative route (for CNV calls) and assuming it is male
	elif (yGender == 'F' and np.round(yCN) == 0.) and (xGender == 'M' and np.round(xCN) == 1.):
		gender = 'M'		

	else:
		gender = yGender

	return gender










def runQCone(sample, species, statsDir, lowessDir, segmentDir, QCdir, plotDir):
	#import config info#
	interpretVars = cfg.Interpret()
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)

	
	
	#determine read counts from the countstats file#
	readCount = calcReads(sample, statsDir)
	
	
	
	#determine MAPD from lowess bins#
	MAPD, lowessData = calcMAPD(sample, lowessDir)
	
	
	
	#determine the optimal ploidy value and CS#
	segData, segArray = common.importSegData(sample, segmentDir, binArray)
	CS, ploidy = getPloidy(segArray, sample, plotDir, interpretVars.QCdict['CS'])
		
	
	
	#determine gender
	gender = getGender(lowessData, binArray['chrom'], ploidy)
	
	
	
	#determine if the sample is good
	qcPass = False
	if readCount >= interpretVars.QCdict['Reads'] and MAPD <= interpretVars.QCdict['MAPD'] and CS >= interpretVars.QCdict['CS']:
		qcPass = True
	
	
	
	#write an intermediate file with name, reads, MAPD, CS, ploidy, gender, good/bad
	outfile = QCdir + sample + '.qcTEMP.txt'
	OUT = open(outfile, 'w')
	outData = [sample, readCount, MAPD, CS, ploidy, gender, qcPass]
	outText = '\t'.join(map(str, outData)) + '\n'
	OUT.write(outText)
	OUT.close()
	
	
	
	printText = '\t\tFinished assessing QC for ' + sample
	print(printText)
	
	
	
	
	
	
	
	
	
	
