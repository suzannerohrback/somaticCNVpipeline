#!/usr/bin/python
import os,sys,inspect
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import config as cfg
import common





def plotProfile(sample, outDir, lowessData, cnvData, refArray):
	xVals = [x['abspos'] + (x['size']/2) for x in refArray]

	chromStarts = [refArray[refArray['chrom'] == x]['abspos'][0] for x in np.unique(refArray['chrom'])]
	print chromStarts
	
#	chromStarts = [min([y['start'] for y in refArray[refArray['chrom'] == x]]) for x in np.unique(refArray['chrom'])]
	chromEnds = [max([y['start'] + y['size'] for y in refArray[refArray['chrom'] == x]]) for x in np.unique(refArray['chrom'])]
	chromEdges = chromEnds[:-1]
	xTicks = [np.mean([chromStarts[x], chromEnds[x]]) for x,y in enumerate(np.unique(refArray['chrom'])[:-1])]

	fig, ax = plt.subplots()

	#ideally make these prettier colors like were used in the paper#
	ax.scatter(xVals, binData, color='b', marker='d', s=3, linewidths=0)
	ax.plot(xVals, cnvData, color='r', lw=1, ls='steps')

	for j in chromEdges:
		ax.plot([j, j], [-1, 5], lw=1, ls='-', color='gray', zorder=0)

	ax.set_xticks(xTicks)
	ax.set_xticklabels(xUse, fontsize=pFonts.fontSizeDict['tick'], fontname=fontType, rotation=45)
	ax.set_xlabel('Genome Location (Chrom)', fontsize=pFonts.fontSizeDict['axis'], fontname=fontType, labelpad=1)
	ax.set_xlim(0, xVals[-1])

	yTicks = [0, 1, 2, 3, 4]
	ax.set_yticks(yTicks)
	ax.set_yticklabels(yTicks, fontsize=pFonts.fontSizeDict['tick'], fontname=fontType)
	ax.set_ylabel('Copy Number', fontsize=pFonts.fontSizeDict['axis'], fontname=fontType, labelpad=1)
	ax.set_ylim(-0.1, 4.6)

	ax.tick_params(direction='out', which='both', pad=0., length=3, top='off', right='off')

	fig.set_size_inches(1.4, 0.9, forward=True)
	plt.subplots_adjust(left=0.13, right=0.98, bottom=0.3, top=0.91)

	plt.savefig(outDir + sample + 'copyNumberProfile.png', dpi=1200, transparent=True)
	plt.savefig(plotNameAlt, dpi=333)

	plt.close()


#function to create cell summary file (one row per sample)

	
	
def analyzeOne(sample, species, cnvDir, lowessDir, plotDir, ploidy, gender):

	interpretVars = cfg.Interpret()
	
	#load reference data#
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)
	xBins = [x for x,y in enumerate(binArray) if y['chrom'] == 'chrX']
	yBins = [x for x,y in enumerate(binArray) if y['chrom'] == 'chrY']
	
	
	
	#load lowess counts and convert to CN state#
	binData = np.loadtxt(lowessDir + common.findInfile(sample, lowessDir))
	binData = (2 ** binData) * ploidy

	
	
	#load CNV data and convert to array form#
	cnvData = np.array( len(binArray) * [2], dtype='int' )
	if gender == 'M':
		cnvData[xBins] = len(xBins) * [1]
		cnvData[yBins] = len(yBins) * [1]
	else:
		cnvData[yBins] = len(yBins) * [0]

	listFile = cnvDir + common.findInfile(sample, cnvDir)
   	listDtype = {'names': ('chrom', 'start', 'end', 'CN'), 'formats': ('S10', 'int', 'int', 'int')}
	if os.stat(listFile).st_size > 32:
		cnvs = np.loadtxt(listFile, skiprows=1, dtype=listDtype)
		cnvs = np.atleast_1d(cnvs)
		for j in cnvs:
			cnvData[j['start']:j['end']] = j['CN']
			
			
	plotProfile(sample, plotDir, binData, cnvData, binArray)






