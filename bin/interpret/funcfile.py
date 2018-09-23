#!/usr/bin/python
import os,sys,inspect
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





def mergeSegCN(seg1, seg2, intD=False):
	currentWeight = seg1['end'] - seg1['start']
	nextWeight = seg2['end'] - seg2['start']
	weightedAverageCN = np.average([seg1['CN'], seg2['CN']], weights = [currentWeight, nextWeight])
	
	if intD:
		currentIntdist = abs(np.round(seg1['CN']) - seg1['CN'])
		nextIntdist = abs(np.round(seg2['CN']) - seg2['CN'])
		weightedAverageIntdist = np.average([currentIntdist, nextIntdist], weights = [currentWeight, nextWeight])
		mergeIntdist = abs(np.round(weightedAverageCN) - weightedAverageCN)
		
		return weightedAverageCN, weightedAverageIntdist, mergeIntdist

	else:
		return weightedAverageCN
	
	
	
	

def mergeCNinitial(dataDict, gender):
	newData = [dataDict[0]]

	for i in dataDict[1:]:
		addData = i

		if newData[-1]['chrom'] == i['chrom'] and np.round(newData[-1]['CN']) == np.round(i['CN']):
			weightedAverageCN, weightedAverageIntdist, mergeIntdist = mergeSegCN(newData[-1], i, intD=True)

			#previously was <= but that's stupid b/c it can really only improve
			if mergeIntdist < weightedAverageIntdist or np.round(i['CN']) == getNormalCN(i['chrom'], gender):
			#	print 'passed', weightedAverageIntdist, mergeIntdist
				prevMerge = newData.pop()
				addData =	{
					'chrom': i['chrom'], 
					'start': prevMerge['start'], 
					'end': i['end'], 
					'CN': weightedAverageCN
					}
			#else:
			#	print 'failed', weightedAverageIntdist, mergeIntdist

		newData.append(addData)

	return newData





def FUnC(dataDict, binDict, cutoffDict, gender):
	dataDict[-1]['end'] = max(binDict.keys())-1

	#for each entry, calc bin size, compare to cutoffDict
	for i in range(len(dataDict)):
		normalCN = getNormalCN(dataDict[i]['chrom'], gender)
		thisSize = binDict[dataDict[i]['end']+1] - binDict[dataDict[i]['start']]
		dataDict[i]['bins'] = thisSize
		
		if np.round(dataDict[i]['CN']) == normalCN:
			dataDict[i]['pass'] = 'eup'
		elif abs(np.round(dataDict[i]['CN']) - dataDict[i]['CN']) <= cutoffDict[thisSize]:
			dataDict[i]['pass'] = 'cnv'
		else:
			dataDict[i]['pass'] = 'no'
			
	return dataDict





def mergePassing(funcDict):
	mergePass = [funcDict[0]]
	for i, j in enumerate(funcDict[1:]):
		thisEntry = j
		if j['pass'] == funcDict[i]['pass'] == 'cnv': #both entries passed FUnC
			if j['chrom'] == funcDict[i]['chrom']: #both entries on the same chromosome
				if np.round(j['CN']) == np.round(funcDict[i]['CN']): #both entries have the same copy number
					prev = mergePass.pop()
					thisEntry = {
						'chrom': j['chrom'],
						'start': prev['start'], 
						'end': j['end'],
						'CN': mergeSegCN(prev, j),
						'bins': prev['bins'] + j['bins'],
						'pass': 'cnv'
					}	
		mergePass.append(thisEntry)
		
	return mergePass




def mergeCNfinal(funcDict):
	
	print 'starting with', len(funcDict), 'segments'
	
	#pass 1: merge passing CNVs if same CN and same chrom 
	mergePass = mergePassing(funcDict)
	print 'after first merge of passing CNVs', len(mergePass), 'segments remain'
		
	#pass 2: merge small segments with most similar adjacent (passing or euploid) neighbor#
	mergeSmall = []
	skipTest = False
	for i,j in enumerate(mergePass):
		if skipTest:
			skipTest = False
			continue
		thisEntry = j
		mergeTest = False

		if j['bins'] < 3: #segment is unreliably small
	#		print 'SMALL SEG NEEDED MERGING'

			#situations where you automatically merge with the previous segment
			if (i == len(mergePass) - 1) or (j['chrom'] != mergePass[i+1]['chrom']) or (mergePass[i+1]['pass'] == 'no' != mergeSmall[-1]):
		#		print 'auto merge with previous segment'
				mergeTest = i-1
				
			#situations where you automatically merge with the next segment
			elif (i == 0) or (j['chrom'] != mergeSmall[-1]['chrom']) or (mergeSmall[-1]['pass'] == 'no' != mergePass[i+1]):
		#		print 'auto merge with next segment'
				mergeTest = i+1
			
			#situations where you need to figure out which segment to merge with
			else:
				if mergeSmall[-1]['pass'] == mergePass[i+1]['pass'] == 'no':
					if mergeSmall[-1]['bins'] > mergePass[i+1]['bins'] or (mergeSmall[-1]['bins'] == mergePass[i+1]['bins'] and abs(j['CN'] - mergeSmall[-1]['CN']) < abs(j['CN'] - mergePass[i+1]['CN'])):
				#		print 'merge with non passing but larger or CN matching prev seg'
						mergeTest = i-1
					else:
				#		print 'merge with non passing but larger or CN matching next seg'
						mergeTest = i+1
				elif abs(j['CN'] - mergeSmall[-1]['CN']) < abs(j['CN'] - mergePass[i+1]['CN']):
					mergeTest = i-1
				else:
					mergeTest = i+1

		#	print mergeSmall[-1]
		#	print j
		#	print mergePass[i+1]
			
			#actually do the merging
			if mergeTest == i-1:
				parent = mergeSmall.pop()
				thisEntry = {
					'start': parent['start'],
					'end': j['end'],
				}
			elif mergeTest == i+1:
				skipTest = True
				parent = mergePass[mergeTest]
				thisEntry = {
					'start': j['start'],
					'end': parent['end'],
				}

			else:
				print('ERROR: why was no segment found for merging?')
			
			thisEntry['chrom'] = j['chrom']
			thisEntry['CN'] = mergeSegCN(j, parent)
			thisEntry['bins'] = j['bins'] + parent['bins']
			thisEntry['pass'] = parent['pass']

		#	print mergeTest - i, skipTest
		#	print thisEntry
		#	print '\n'
			
		mergeSmall.append(thisEntry)
		
	print 'after merging excessively small segments', len(mergeSmall), 'segments remain'

	#pass 3: merge passing CNVs again if same CN and same chrom because merging small segments can bring them into contact
	remergePass = mergePassing(mergeSmall)
	print 'after second merge of passing CNVs', len(remergePass), 'segments remain'
	

	"""		
	Third, small segments (<= 25 bins) that fail as CNV calls are either
		Combined with neighboring CNVs IF 
			Both are on the same chromosome
			Both have the same copy number
			Both pass FUnC
			Both are > 25 bins
			(because this indicates the baseline CN in this locus is not euploid)
		Otherwise, converted to euploid
	Large segments (>25 bins) are automatically treated as euploid too
	"""
	
	#pass 4: <=25 bin segments that fail as CNV calls are merged with neighboring CNVs if they show a baseline CN shift
	baselineMerge = []
	skipTest = False
	for i,j in enumerate(remergePass):
		if skipTest:
			skipTest = False
			continue
		thisEntry = j

		if j['pass'] == 'no' and j['bins'] <= 25: #segment is a candidate for additional merging
			
			if i == 0 or j['chrom'] != baselineMerge[-1]['chrom']: #segment at 5' end of chromosome#
				if remergePass[i+1]['pass'] == 'cnv' and remergePass[i+1]['bins'] > 25:
					if np.round(remergePass[i+1]['CN']) == np.round(mergeSegCN(j, remergePass[i+1])): #extra check b/c at chrom edge
			#			print 'At the 5prime chrom end and passes merging requirements'
			#			print j
			#			print remergePass[i+1]
						thisEntry = remergePass[i+1]
						thisEntry['bins'] = thisEntry['bins'] + j['bins']
						thisEntry['start'] = j['start']
						skipTest = True	
			#			print thisEntry, '\n'
			
			if i == len(remergePass)-1 or j['chrom'] != remergePass[i+1]['chrom']: #segment at 3' end of chromosome#
				if baselineMerge[-1]['pass'] == 'cnv' and baselineMerge[-1]['bins'] > 25:
					if np.round(baselineMerge[-1]['CN']) == np.round(mergeSegCN(j, baselineMerge[-1])): #extra check b/c at chrom edge
			#			print 'At the 3prime chrom end and passes merging requirements'
			#			print baselineMerge[-1]
			#			print j
						thisEntry = baselineMerge.pop()
						thisEntry['bins'] = thisEntry['bins'] + j['bins']
						thisEntry['end'] = j['end']
			#			print thisEntry, '\n'

			else: #segment in middle of chromosome#
				if baselineMerge[-1]['pass'] == remergePass[i+1]['pass'] == 'cnv': #both are passing CNVs#
					if np.round(baselineMerge[-1]['CN']) == np.round(remergePass[i+1]['CN']): #both are the same CN#
						if baselineMerge[-1]['bins'] > 25 and remergePass[i+1]['pass'] > 25: #both are more than 25 bins#
							thisEntry = baselineMerge.pop()
							thisEntry['bins'] = thisEntry['bins'] + j['bins'] + remergePass[i+1]['bins']
							thisEntry['end'] = remergePass[i+1]['end']
							thisEntry['CN'] = mergeSegCN(thisEntry, remergePass[i+1])
							skipTest = True
			
		baselineMerge.append(thisEntry)
	
	
	print 'after merging  small segments surrounded by baseline copy number shift', len(baselineMerge), 'segments remain'
	
	cnvData = [x for x in baselineMerge if x['pass'] == 'cnv']
	return cnvData










def FUnCone(sample, species, segmentDir, CNVdir, ploidy, gender):
	
	#import config info#
	interpretVars = cfg.Interpret()
	refArray = common.importInfoFile(interpretVars.binDict[species], [0, 1, 2, 4, 5], 'normref', skiprows=1)

	binDict = {y:x for x,y in enumerate(refArray['abspos'])} #make dict of ref array so I can calculate bin number
	binDict[refArray[-1]['abspos'] + refArray[-1]['size']+1] = len(refArray) #fix for last bin not having the correct end position#
	
	cutoffArray = np.loadtxt(interpretVars.cutoffFile, skiprows=1, dtype={'names': ('bins', 'intD'), 'formats': ('int', 'float')})
	cutoffDict = {x['bins']: x['intD'] for x in cutoffArray}
	
	
	
	#load segment data#
	segData, segArray = common.importSegData(sample, segmentDir, refArray)
	segData['CN'] = segData['CN'] * ploidy
	dataDict = [ {y: x[y] for y in segData.dtype.names} for x in segData ]

	
	
	#merge adjacent segments that have the same copy number, when merging improves intD#
	mergeDataDict = mergeCNinitial(dataDict, gender)
	
	
	
	#run FUnC#
	funcDataDict = FUnC(mergeDataDict, binDict, cutoffDict, gender)
	numFail = [x['pass'] for x in funcDataDict].count('no')
	numPass = [x['pass'] for x in funcDataDict].count('cnv')
	

	
	#second merge#
	cnvData = mergeCNfinal(funcDataDict)

	
	
	#write output file#
	outfile = CNVdir + sample + 'CNVlist.bed'
	OUT = open(outfile, 'w')
	OUT.write('Chromosome\tStart\tEnd\tCopyNumber\n')
	for i in cnvData:
		print i
		OUT.write(i['chrom'])
		OUT.write('\t')
		OUT.write(str(refArray[binDict[i['start']]]['chrStart']))
		OUT.write('\t')
		OUT.write(str(refArray[binDict[i['end']+1]]['chrStart']-1))
		OUT.write('\t')
		OUT.write(str(np.round(i['CN'])))
		OUT.write('\n')
	OUT.close

	
	
	printText = '\t\tFinished performing FUnC on ' + sample + ', removed ' + str(numFail) + ' of ' + str(numFail + numPass) + ' putative CNV segments'
	printText += '\n\t\t\tFor a total of ' + str(len(cnvData)) + ' CNVs after merging'
	print(printText)
	raise SystemExit

	
	
	
	
	
