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





def plotProfile(sample, outDir, lowessData, cnvData, refArray, chromList):
	xVals = [x['abspos'] + (x['size']/2) for x in refArray]

	chromStarts = [refArray[refArray['chrom'] == x]['abspos'][0] for x in chromList]
	chromEnds = [refArray[refArray['chrom'] == x]['abspos'][-1] + refArray[refArray['chrom'] == x]['size'][-1] for x in chromList]

	chromEdges = chromEnds[:-1]
	xTicks = [np.mean([chromStarts[x], chromEnds[x]]) for x in range(len(chromList))]
	yTicks = [0, 1, 2, 3, 4]


	
	fig, ax = plt.subplots()

	#ideally make these prettier colors like were used in the paper#
	ax.scatter(xVals, lowessData, color='#247afd', marker='d', s=3)
	ax.plot(xVals, cnvData, color='#980002', lw=1, ls='steps')

	for j in chromEdges:
		ax.plot([j, j], [-1, 5], lw=1, ls='-', color='#6b7c85', zorder=0)

	ax.set_xticks(xTicks)
	ax.set_xticklabels(chromList, rotation=45)
	ax.set_xlim(-10000000, xVals[-1]+10000000)

	ax.set_yticks(yTicks)
	ax.set_yticklabels(yTicks)
	ax.set_ylabel('Copy Number', labelpad=5)
	ax.set_ylim(-0.1, 4.6)

	ax.tick_params(direction='out', which='both', pad=0., length=3, top='off', right='off')

	fig.set_size_inches(8, 4, forward=True)
	plt.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.95)
	plt.savefig(outDir + sample + '.copyNumberProfile.png', dpi=666)

	plt.close()
	
	
	
	
#plot estimated overall chromosome ploidy states--CALCULATED FROM LOWESS DATA, NOT SEGMENT DATA#
def plotChroms(sample, outDir, lowessData, refArray, chromList):
	graphData = [np.mean(lowessData[refArray['chrom'] == x]) for x in chromList]
	graphErr = [np.std(lowessData[refArray['chrom'] == x]) for x in chromList]
	graphErr = [abs(graphData[x] - graphErr[x]) for x in range(len(chromList))]

	xTicks = np.arange(0, len(chromList))
	yTicks = [0, 1, 2, 3, 4]

	fig, ax = plt.subplots()
	
	ax.bar(xTicks, graphData, width=0.8, color='#c0022f', align='center', yerr=graphErr, ecolor='k')
	
	ax.set_xticks(xTicks)
	ax.set_xticklabels(chromList, rotation=45)
	ax.set_xlim(-1, xTicks[-1]+1)

	ax.set_yticks(yTicks)
	ax.set_yticklabels(yTicks)
	ax.set_ylabel('Copy Number', labelpad=5)
	ax.set_ylim(0, 4.6)
	
	ax.tick_params(direction='out', which='both', pad=0., length=3, top='off', right='off')
	
	fig.set_size_inches(8, 4, forward=True)
	plt.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.95)
	plt.savefig(outDir + sample + '.chromosomeCopyNumber.png', dpi=666)
	plt.close()

	
	
	
	
#function to calculate cell summary stats
def getSummaryStats(cnvs, gender, chromList, chromSizes,noCNVs=False):
	cellStats = {
		'delCount': 0,
		'ampCount': 0,
		'delMB': 0.,
		'ampMB': 0.,
	}
	chromAmp = {x: 0. for x in chromList if x != 'chrY'}
	chromDel = {x: 0. for x in chromList if x != 'chrY'}
	if not noCNVs:
		for i in cnvs:
			normalCN = common.getNormalCN(i['chrom'], gender)
		
			if i['chrom'] == 'chrY':
				break
		
			if i['CN'] < normalCN:
				cellStats['delCount'] += 1
				cellStats['delMB'] += float(abs(normalCN - i['CN']) * (i['end'] - i['start'] + 1)) / 1e6
				chromDel[i['chrom']] += float(i['end'] - i['start'] + 1)
			else:
				cellStats['ampCount'] += 1
				cellStats['ampMB'] += float(abs(normalCN - i['CN']) * (i['end'] - i['start'] + 1)) / 1e6
				chromAmp[i['chrom']] += float(i['end'] - i['start'] + 1)
			
		chromAmp = {y: (100. * chromAmp[y]) / float(chromSizes[x]) for x,y in enumerate(chromList) if y != 'chrY'}
		chromDel = {y: (100. * chromDel[y]) / float(chromSizes[x]) for x,y in enumerate(chromList) if y != 'chrY'}

		thisResult = {'chroms': [x for x in chromList if x != 'chrY'], 'cellStats': cellStats, 'chromAmp': chromAmp, 'chromDel': chromDel}
	return thisResult
	
	
	
	
	
	
	
	
	
	
def analyzeOne(sample, species, cnvDir, lowessDir, plotDir, chromPlotDir, ploidy, gender):

	interpretVars = cfg.Interpret()
	
	#load reference data#
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)
	xBins = [x for x,y in enumerate(binArray) if y['chrom'] == 'chrX']
	yBins = [x for x,y in enumerate(binArray) if y['chrom'] == 'chrY']
	
	chromList = ['chr1'] + [y for x,y in enumerate(binArray['chrom'][1:]) if y != binArray['chrom'][x]]
	chromSizes = [binArray[binArray['chrom'] == x]['chrStart'][-1] + binArray[binArray['chrom'] == x]['size'][-1] for x in chromList]
		#this might need an extra -1
	
	
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
			startBin = [x for x,y in enumerate(binArray) if y['chrom'] == j['chrom'] and y['chrStart'] == j['start']][0]
			endBin = [x for x,y in enumerate(binArray) if y['chrom'] == j['chrom'] and y['chrStart'] + y['size'] - 1 == j['end']][0]
			cnvData[startBin:endBin] = j['CN']
		summaryStats = getSummaryStats(cnvs, gender, chromList, chromSizes)			
	else:
		cnvs=0
		summaryStats = getSummaryStats(cnvs, gender, chromList, chromSizes, noCNVs=True)
	plotProfile(sample, plotDir, binData, cnvData, binArray, chromList)
	plotChroms(sample, chromPlotDir, binData, binArray, chromList)
	#summaryStats = getSummaryStats(cnvs, gender, chromList, chromSizes)
	
	return summaryStats




