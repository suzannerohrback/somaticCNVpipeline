#!/usr/bin/python
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg





def getNormalCN(chrom, gender):
	normalCN = 2.
	if chrom in ['chrX', 'chrY'] and gender == 'M':
		normalCN = 1.
	elif chrom == 'chrY' and gender == 'F':
		normalCN = 0.
	return normalCN





def mergeCNinitial(dataDict):
	newData = [dataDict[0]]

	for i in dataDict[1:]:
		addData = i

		if newData[-1]['chrom'] == i['chrom'] and np.round(newData[-1]['CN']) == np.round(i['CN']):
			currentWeight = newData[-1]['end'] - newData[-1]['start']
			nextWeight = i['end'] - i['start']
			weightedAverageCN = np.average([newData[-1]['CN'], i['CN']], weights = [currentWeight, nextWeight])
			mergeIntdist = abs(np.round(weightedAverageCN) - weightedAverageCN)

			currentIntdist = abs(np.round(newData[-1]['CN']) - newData[-1]['CN'])
			nextIntdist = abs(np.round(i['CN']) - i['CN'])
			weightedAverageIntdist = np.average([currentIntdist, nextIntdist], weights = [currentWeight, nextWeight])

			if mergeIntdist <= weightedAverageIntdist:
				prevMerge = newData.pop()
				addData =	{
					'chrom': i['chrom'], 
					'start': prevMerge['start'], 
					'end': i['end'], 
					'CN': weightedAverageCN
					}
		newData.append(addData)

	return newData





def FUnC(dataDict, refArray, cutoffDict, gender):
	#make dict of ref array so I can calculate out bin number
	binDict = {x: y for x,y in enumerate(refArray['abspos'])}
	binDict[len(refArray)] = refArray[-1]['abspos'] + refArray[-1]['size']
	
	#for each entry, calc bin size, compare to cutoffDict
	for i in dataDict:
		normalCN = getNormalCN(i['chrom'], gender)
		###I AM HERE###
	
	#if passes, keep CN, if not, convert to euploid (TAKE GENDER INTO ACCOUNT)
	#return dict
	
	return 0





def mergeCNfinal():
	#merging
	#convert to np array at the end (condensed format)
	return 0





def runFUNCone(sample, species, segmentDir, CNVdir, ploidy, gender):
	
	#import config info#
	interpretVars = cfg.Interpret()
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 2, 4, 5], 'normref', skiprows=1)
	cutoffArray = np.loadtxt(interpretVars.cutoffFile, skiprows=1, dtype={'names': ('bins', 'intD'), formats: ('int', 'float')})
	cutoffDict = {x['bins']: x['intD'] for x in cutoffArray}
	
	#load segment data#
	segData, segArray = common.importSegData(sample, segmentDir, binArray)
	segData['CN'] = segData['CN'] * ploidy
	dataDict = [ {y: x[y] for y in segData.dtype['names']} for x in segData ]

	#merge adjacent segments that have the same copy number, when merging improves intD#
	mergeDataDict = mergeCNinitial(dataDict)
	
	#run FUnC#
	
	
	#second merge#
	
	
	#write output file#
	
	
	
	printText = '\t\tFinished performing FUnC on ' + sample
	print(printText)

	
	
