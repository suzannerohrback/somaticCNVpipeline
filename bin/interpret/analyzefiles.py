#!/usr/bin/python
import os,sys,inspect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import config as cfg





def plotProfile(sample, outDir, lowessData, cnvData, refArray):
	xVals = [x['abspos'] + (x['size']/2) for x in refArray]

	chromStarts = [min([y['start'] for y in refArray[refArray['chrom'] == x]]) for x in np.unique(refArray['chrom'])]
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

#function to creat CNV summary file (one line per CNV)
	
	
def analyzeOne(sample, segDir, lowessDir, outDir, ploidy, gender):

	interpretVars = cfg.Interpret()
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)

	#this is probably the wrong filename
	binData = np.loadtxt(lowessDir + sample + '.lowessBinCounts.txt')
	binData = (2 ** binData[:len(refData)]) * ploidyDict[i]

	
	
	listFile = './data/' + i + '.CNVlist.bed'
	cnvData = np.array( len(refData) * [2], dtype='int' )
	if genderDict[i] == 'M':
		cnvData[xBins] = len(xBins) * [1.]

	#THIS NUMBER ALMOST CERTAINLY NEEDS TO BE SMALLER
   	listDtype = {'names': ('chrom', 'start', 'end', 'CN', 'type'), 'formats': ('S10', 'int', 'int', 'int', 'S10')}
	if os.stat(listFile).st_size > 47:
		cnvs = np.loadtxt(listFile, skiprows=1, dtype=listDtype)
		cnvs = np.atleast_1d(cnvs)
		for j in cnvs:
			cnvData[j['start']:j['end']] = j['CN']
			
			
	plotProfile(sample, outDir, binData, cnvDat, binArray)






