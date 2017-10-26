#!/usr/bin/python
import random
import os
import subprocess as sub
import shlex
import multiprocessing as mp
import numpy as np

##########Modules specific to SR's environment
#import config as cfg
#import parallel





#THIS CODE WILL TAKE SOME EFFORT TO MODIFY FOR A NEW ENVIRONMENT#





#Notes created during development--no guarantee they are all relevant
#####use the X chromosome of males and females to assess FDR and FNR, particularly using my cutoffs

	#DONE!!!##import reference information
		#DONE#sample names, QC info, (ploidy?), genders
		#DONE#cnv cutoffs
		#DONE#determine which QC-passing samples have no putative CNVs on the X chromosome
			#DONEat the same time, build a information about CNVs on the X chromosome of female cells
				#DONEconsider both all putative CNVs, DONE and filtered CNVs

	###FDR will be assessed by adding together the data of 2 male X chromosomes and seeing if CNVs are called (expect 0!)
		#DONE#using samples where no CNVs called originally

		#DONE#don't want to merge samples with extremely different QC metrics...(or do I???)
			#DONEcalculate range-scaled distance matrix of QC values
			#DONEfor each sample, find the 10 others that are most similar to it
			#DONEfilter down the list to only contain unique combinations
				###THIS GIVES 582 COMPARISONS
				###Will this be enough for convincing comparisons? Might need to expand based on the output I get

		#DONE#for each pair
			#generate a new lowess file by adding together the X chrom lowess bin counts
			#run CBS
			#convert the segment data to linear form 
				###2 ** ratio * 2 should be fine here because there's only one chromosome and it was euploid to start...but I might want to confirm that
			#look for the presence of any CNVs 
				#merge if needed
			#assess if any CNVs would be removed by my filtering
			#return the number of CNVs with and without filtering

		#DONE#statistics to assess
			#how many CNVs does filtering remove (both for females and simulated females)
			#FDR
				#FP rate can be the total number of CNVs from simFem runs / number of runs (bad calls per cell)
				#DONEthen the TP+FP is the number of CNVs from realFems / number of realFems  (calls made per cell)
				#FDR is FP/(TP+FP)
				#look at for both putative and filtered CNVs, and compare the difference between the two

	###FNR will be assessed by combining portions of male and female X chromosomes
		#DONE#using samples where no CNVs called originally

		#DONE#general
			#losses simulated by replacing a random region in the female X with a random set of contiguous male X bins
			#gains simulated by adding together the randomly selected contiguous regions
			#use lowess bin counts so differences in sequencing depth shouldn't be an issue (adjust both for ploidy, then use the female ploidy throughout)

		#DONE#again don't want to merge samples with extremely different QC metrics...(or do I???)
			#DONEcalculate range-scaled distance matrix of QC values
			#DONEfor each female sample, find the 10 males that are most similar to it
			#DONEfilter down the list to only contain unique combinations
				###THIS GIVES 1160 COMPARISONS
				###Will this be enough for convincing comparisons? Might need to expand based on the output I get

		#DONE#for each pair
			#randomly select a cnv size, change direction, and start location
				#will need to come up with a list of CNV sizes of interest
				#change could only be gain or loss
				#start could be anywhere that lets the whole length fit
			#mash in male data however is appropriate (see above)
			#save a lowess file and run CBS
			#assess if the CNV was identified
				#at least 50% reciprocal overlap
				#correct copy number
				#assess if it would be filtered out with my cutoffs or not
			#save information: size, copy number, putativeID, filteredID
			#repeat these steps some number of times to assess a range of size/location options with the combination

		##statistics to assess
			#DONEoverall FNR
				#number of times sample was not ID'd / number of trials
				#for both putative and filtered ID rates
			#check for a size dependency of FNR
				#if yes, I might want to adjust calculations for reporting purposes...or adjust filtering...we will see...
			









###calculate euclidian distance of QC variables between two samples###
def getSampleDist(data1, data2):
	thisSum = 0
	for i in data1.keys():
		thisSum += (data1[i] - data2[i]) ** 2
	eucDist = thisSum ** 0.5
	return eucDist










def getRefInfo(name, xStart, xEnd, absPosDict):
	print '\n\n\nPreparing to assess FDR and FNR using the X chromosome of male and female samples'
	segDtype = {'names': ('chrom', 'start', 'end'), 'formats': ('S15', 'S15', 'S15')}

	###get a list of single cell samples###NEEDS UPDATING TO UNIVERSAL SYSTEM
	infile = folders.getInfoFile(name)
	data = np.loadtxt(infile, skiprows=1, usecols=[0,7,10], dtype={'names': ('name', 'cells', 'method'), 'formats': ('S10', 'int', 'S10')})
	useCells = data[data['cells'] == 1]
	useCells = useCells[useCells['method'] == 'VeriSeq']
	useCells = list(useCells['name'])





	###import currently used CNV filtering cutoffs###
	cutoffDict = {}
  ##########need to change filepaths
	#cutoffFile1 = '/filepath/mm10.smallThresholdCutoffs.25k.bowtie.k36.txt'
	#cutoffFile2 = '/filepath/mm10.largeThresholdCutoffs.25k.bowtie.k36.txt'
	for i in [cutoffFile1, cutoffFile2]:
		data = np.loadtxt( i, skiprows=1, dtype={'names': ('size', 'intdist'), 'formats': ('int', 'float')} )
		for j in data:
			cutoffDict[j['size']] = j['intdist']
	for i in range(1, min(cutoffDict.keys())):
		cutoffDict[i] = 0.0





	###Import QC data for all samples######NEEDS UPDATING TO UNIVERSAL SYSTEM
	QCfile = folders.getQCfile(name)
	qcData = np.loadtxt(QCfile, dtype={'names': ('Sample', 'Reads', 'MAPD', 'CS', 'Ploidy', 'Gender'), 'formats': ('S10', 'float', 'float', 'float', 'float', 'S1')}, skiprows=1, usecols=[0, 1, 2, 3, 4, 5])
	qcData = qcData[qcData['Reads'] >= 600000]
	qcData = qcData[qcData['MAPD'] <= 0.38]
	qcData = qcData[qcData['CS'] >= 0.83]
	ploidyDict = {x['Sample']: x['Ploidy'] for x in qcData}
	sampleList = [x for x in qcData['Sample'] if x in useCells]





	###Assess which samples have no CNVs on the X chromosome, and save references of the CNVs on Female X chromosomes###
	maleCleanList = []
	femaleCleanList = []
	femaleCells = [x['Sample'] for x in qcData if x['Gender'] == 'F' and x['Sample'] in sampleList]
	putativeFemaleCNVs = 0
	putFemLarge = 0
	putFemSmall = 0
	filteredFemaleCNVs = 0
	filtFemLarge = 0
	filtFemSmall = 0

	for i in sampleList:
		normalCN = 1.
		if i in femaleCells:
			normalCN = 2.





  ##########need to update filepaths
			#listFile = /filepath/cnvLists/' + i + '.CNVlist.txt'
			data = np.loadtxt(listFile, dtype=segDtype, usecols=[0, 1, 2,])#, skiprows=1)
			data = np.atleast_1d(data)
			data = data[data['chrom'] == 'chrX']
			filteredFemaleCNVs += len(data)
			filtFemLarge += sum([1 for x in data if int(x['end']) - int(x['start']) > sizeBoundary])
			filtFemSmall += sum([1 for x in data if int(x['end']) - int(x['start']) < sizeBoundary])

	


  ##########need to update filepaths
		#segFile = '/filepath/' + i + 'k25.segments.txt'
		data = np.loadtxt(segFile, dtype={'names': ('start', 'end', 'CN'), 'formats': (int, int, float)})
		usableSegs = [x for x,y in enumerate(data) if y['end'] > y['start'] and xStart <= y['start'] <=xEnd]
		data = data[usableSegs]
		data['CN'] = [2.0 if np.isinf(x) else (2 ** x) * ploidyDict[i] for x in data['CN']]
		cnvSegs = data[np.round(data['CN']) != normalCN]

		if i not in femaleCells:
			if len(cnvSegs) == 0:
				maleCleanList.append(i)

		else:
			if len(cnvSegs) == 0:
				femaleCleanList.append(i)
			else:
				mergeData = [cnvSegs[0]]
 				putativeFemaleCNVs += len(mergeData)
				putFemLarge += sum([1 for x in mergeData if absPosDict[x[1]+1] - absPosDict[x[0]] > sizeBoundary])
				putFemSmall += sum([1 for x in mergeData if absPosDict[x[1]+1] - absPosDict[x[0]] < sizeBoundary])

	putativeFemaleRate = float(putativeFemaleCNVs) / float(len(femaleCells))
	putLargeRate = float(putFemLarge) / float(len(femaleCells))
	putSmallRate = float(putFemSmall) / float(len(femaleCells))
	filteredFemaleRate = float(filteredFemaleCNVs) / float(len(femaleCells))
	filtLargeRate = float(filtFemLarge) / float(len(femaleCells))
	filtSmallRate = float(filtFemSmall) / float(len(femaleCells))

	print '\t', len(femaleCells), 'female cells have'
	print '\t\t', putativeFemaleCNVs, 'putative and', filteredFemaleCNVs, 'filtered X chromosome CNVs'
	print '\t\tfor respective rates of', np.round(putativeFemaleRate, 3), 'and', np.round(filteredFemaleRate, 3), 'per cell'

	print '\tThere are', len(maleCleanList), 'males and', len(femaleCleanList),  'females with no X chromosome CNVs'





	###Determine which samples to combine based on their difference in QC variables###
		#I am calculating all pairwise differences in QC statistics, then combining paris of samples with differences in the bottom 15% for males, and in the bottom 5% for females#
	qcCompareDict = { x['Sample']: { y: (x[y] - min(qcData[y])) / (max(qcData[y]) - min(qcData[y])) for y in ['Reads', 'MAPD', 'CS'] } for x in qcData }
	malePercent = 15
	femalePercent = 5

	allComparisons = sum([ [(y, z, getSampleDist(qcCompareDict[y], qcCompareDict[z])) for z in maleCleanList[x+1:]] for x,y in enumerate(maleCleanList[:-1]) ], [])
	allComparisons = np.array(allComparisons, dtype={'names': ('cell1', 'cell2', 'dist'), 'formats': ('S10', 'S10', 'float')})
	FDRcutoff = np.percentile(allComparisons['dist'], malePercent)
	FDRcomparisons = allComparisons[allComparisons['dist'] <= FDRcutoff]
	print '\t\t', len(FDRcomparisons), 'male X combinations for FDR testing with a distance range of', min(FDRcomparisons['dist']), np.mean(FDRcomparisons['dist']), FDRcutoff


	allComparisons = sum([ [(y, z, getSampleDist(qcCompareDict[y], qcCompareDict[z])) for z in maleCleanList] for x,y in enumerate(femaleCleanList) ], [])
	allComparisons = np.array(allComparisons, dtype={'names': ('cell1', 'cell2', 'dist'), 'formats': ('S10', 'S10', 'float')})
	FNRcutoff = np.percentile(allComparisons['dist'], femalePercent)
	FNRcomparisons = allComparisons[allComparisons['dist'] <= FNRcutoff]
	print '\t\t', len(FNRcomparisons), 'female X combinations for FNR testing with a distance range of', min(FNRcomparisons['dist']), np.mean(FNRcomparisons['dist']), FNRcutoff


	femaleRateDict =	{
			'putative':	{
					'all': putativeFemaleRate,
					'large': putLargeRate,
					'small': putSmallRate,
					},
			'filtered':	{
					'all': filteredFemaleRate,
					'large': filtLargeRate,
					'small': filtSmallRate,
					},
			}

	return cutoffDict, ploidyDict, FDRcomparisons, FNRcomparisons, femaleRateDict










#loads lowess bin count data, adjusts it to copy number form, and returns it#
def loadCountData(name, cell, numBins, binLocs, ploidy):
	file = folders.getLowessFile(name, cell, numBins)
	data = np.loadtxt(file)
	data = (2 ** data[binLocs]) * ploidy
	return data










def loadSegData(segFile, usePloidy):
	segDtype = dtype={'names': ('start', 'end', 'CN'), 'formats': ('int', 'int', 'float')}

	segData = np.loadtxt(segFile, dtype=segDtype)
	segData = segData[segData['CN'] != np.inf]
	segData = segData[segData['CN'] != -np.inf]
	segData['CN'] = (2 ** segData['CN']) *usePloidy

	if len(segData) > 1:
		mergeData = []
		mergeData.append(segData[0])
		for i in segData[1:]:
			newData = i
			if np.round(i['CN']) == np.round(mergeData[-1][2]):
				old = mergeData.pop()
				newCN = np.average([old[2], i['CN']], weights=[old[1]-old[0], i['end'] - i['start']])
				newData = (old[0], i['end'], newCN)
			mergeData.append(newData)

		segData = np.array(mergeData, dtype=segDtype)

	return segData













def runOneFDR(name, numBins, xBins, absPosDict, overwrite, cell1, cell2, ploidy1, ploidy2):
	matlabName = cell1 + 'V' + cell2
  ##########Probably need to update filename
#	refFile = folders.getOutDir(name, 'prep') + 'xChromRef.txt'
	segDtype = dtype={'names': ('start', 'end', 'CN'), 'formats': ('int', 'int', 'float')}
	segFile = folders.getOutDir(name, 'XchromFDR') + matlabName + '.segments.txt'



	#only combine samples and run CBS if necessary#
	if overwrite or not os.path.exists(segFile):

		###Merge X chromosome bin counts and write to file###
		data1 = loadCountData(name, cell1, numBins, xBins, ploidy1)
		data2 = loadCountData(name, cell2, numBins, xBins, ploidy2)

		mergedXdata = data1 + data2
		testXdata = np.log2(mergedXdata / np.median(mergedXdata))

  ##########Probably need to update filename
#		dataFile = folders.getOutDir(name, 'XchromFDR') + matlabName + '.lowessBinCounts.txt'
		np.savetxt(dataFile, testXdata)



		###write Matlab CBS script to file###
  ##########Probably need to update filename
#		scriptFile = folders.getOutDir(name, 'XchromFDR') + matlabName + '.m'
		SCRIPT = open(scriptFile, 'w')

		importLines = "sampleName = genvarname('" + matlabName + "');\n"
		importLines += "sampleData = importdata('" + dataFile + "');\n"
		importLines += "ref = importdata('" + refFile + "');\n"

		refLines = "chroms = zeros(length(ref.textdata), 1);\n\n"
		refLines += "for i = 1:length(ref.textdata)\n    chrom = textscan(ref.textdata{i},'%s %s','delimiter','r');\n    chroms(i) = 20;\nend\n\n"
		refLines += "bins = ref.data;\nlogLowess = sampleData;\nclear ref sampleData\n\n"

		CBSlines = "cbsInput = [chroms, bins, logLowess];\ncbsOutput = cghcbs(cbsInput, 'ALPHA', " + str(0.01) + ", 'PERMUTATIONS', 10000, 'STOPPINGRULE', false, 'SHOWPLOT', false);\n\n"
		CBSlines += "numOfSegments = 0;\nfor i = 1:length(cbsOutput.SegmentData)\n    for j = 1:length(cbsOutput.SegmentData(1, i).Mean)\n        numOfSegments = numOfSegments + 1;\n    end\nend\n\n"
		CBSlines += "segments = zeros(length(numOfSegments), 3);\nk = 0;\nfor i = 1:length(cbsOutput.SegmentData)\n    for j = 1:length(cbsOutput.SegmentData(1, i).Mean)\n        k = k + 1;\n        segments(k, 1) = cbsOutput.SegmentData(1, i).Start(j);\n        segments(k, 2) = cbsOutput.SegmentData(1, i).End(j);\n        segments(k, 3) = cbsOutput.SegmentData(1, i).Mean(j);\n    end\nend\n\n"
		CBSlines += "segFile = fopen([sampleName '.segments.txt'], 'w');\nfor i = 1:numOfSegments\n    fprintf(segFile,'%d\t',segments(i,1));\n    fprintf(segFile,'%d\t',segments(i,2));\n    fprintf(segFile,'%.15f\\n',segments(i,3));\nend\n\n"

		lastLine = "quit\n"

		allLines = '\n'.join([importLines, refLines, CBSlines, lastLine])
		SCRIPT.write(allLines)
		SCRIPT.close()



		###run CBS###
    ##########THIS MAY BE HARD TO RUN ON MANY SYSTEMS, need access to Matlab from the command line
		matlabCommand = "matlab -nodisplay -r " + matlabName
  ##########Probably need to update filename
#		stderrfile = folders.getOutDir(name, 'XchromFDR') + matlabName + '.stderr'

		cmd = shlex.split(matlabCommand)
		stdout = open(stderrfile, 'w')

		p = sub.Popen(cmd, stdout=stdout, stderr=sub.STDOUT)
		p.wait()

		os.remove('./' + matlabName + '.m')
		os.remove('./' + matlabName + '.stderr')



	###assess the presence of false positive CNVs###
	segData = loadSegData(segFile, 2.)
	segCNVs = segData[np.round(segData['CN']) != 2.]


	segInfo = []
	for i in segCNVs:
		binStart = absPosDict[i['start']]
		binEnd = absPosDict[i['end']]
		size = binEnd - binStart
		CN = np.round(i['CN'])
		intdist = abs(i['CN'] - CN)
		segInfo.append( (matlabName, binStart, binEnd, CN, size, intdist) )

#	print '\t', segInfo

	return segInfo





def runFDR(name, overwrite, FDRcomparisons, femaleRateDict, ploidyDict, xBins, numBins, absPosDict, cutoffDict):
	print '\n\n\nAssessing FDR by combining male X chromosome bin counts'
  ##########probably need to update filename
#	os.chdir(folders.getOutDir(name, 'XchromFDR'))
	resultDtype = {'names': ('cells', 'start', 'end', 'CN', 'size', 'intdist'), 'formats': ('S20', 'int', 'int', 'float', 'int', 'float')}
	distanceDict = {x[0] + 'V' + x[1]: x[2] for x in FDRcomparisons}



	argList = [(name, numBins, xBins, absPosDict, overwrite, x[0], x[1], ploidyDict[x[0]], ploidyDict[x[1]]) for x in FDRcomparisons]
  
  pool = mp.Pool(mp.cpu_count())
	processes = [ pool.apply_async(runOneFDR, args=x) for x in argList ]
  pool.close()
  
  results = len(argList) * [ [] ]
	for i,j in enumerate(processes):
		j.wait()
		results[i] = j.get()

  

	print '\nFinished segmenting merged male X chromosome bin counts'
	resultsArray = np.array(sum([x for x in results if len(x) > 0], []), dtype=resultDtype)

	putativeCNVnum = len(resultsArray)
	putativeCNVrate = float(putativeCNVnum) / float(len(FDRcomparisons))
	putLargeRate = sum([1. for x in resultsArray if x['size'] > sizeBoundary]) / float(len(FDRcomparisons))
	putSmallRate = sum([1. for x in resultsArray if x['size'] < sizeBoundary]) / float(len(FDRcomparisons))
	putFDRall = putativeCNVrate / femaleRateDict['putative']['all']
	putFDRlarge = putLargeRate / femaleRateDict['putative']['large']
	putFDRsmall = putSmallRate / femaleRateDict['putative']['small']

	filteredResultsArray = np.array([x for x in resultsArray if cutoffDict[x['size']] >= x['intdist']], dtype=resultsArray.dtype)
	filteredCNVnum = len(filteredResultsArray)
	filteredCNVrate = float(filteredCNVnum) / float(len(FDRcomparisons))
	filtLargeRate = sum([1. for x in filteredResultsArray if x['size'] > sizeBoundary]) / float(len(FDRcomparisons))
	filtSmallRate = sum([1. for x in filteredResultsArray if x['size'] < sizeBoundary]) / float(len(FDRcomparisons))
	filtFDRall = filteredCNVrate / femaleRateDict['filtered']['all']
	filtFDRlarge = filtLargeRate / femaleRateDict['filtered']['large']
	filtFDRsmall = filtSmallRate / femaleRateDict['filtered']['small']

	percentDiffAll = (filtFDRall - putFDRall) / putFDRall
	percentDiffLarge = (filtFDRlarge - putFDRlarge) / putFDRlarge
	percentDiffSmall = (filtFDRsmall - putFDRsmall) / putFDRsmall

	foldDiffAll = putFDRall / filtFDRall
	try:
		foldDiffLarge = putFDRlarge / filtFDRlarge
	except ZeroDivisionError:
		foldDiffLarge = 'infinity'
	foldDiffSmall = putFDRsmall / filtFDRsmall

	print '\tAfter', len(FDRcomparisons), 'iterations'
	print '\t\t', putativeCNVnum, 'total putative CNVs were identified for a rate of', putativeCNVrate, 'CNVs per cell'
	print '\t\t\tComparing to', femaleRateDict['putative']['all'], 'this results in an FDR of', np.round(100 * putFDRall, 3), 'percent of all putative diploid X chrom calls'
	print '\t\t\tComparing', putLargeRate, 'to', femaleRateDict['putative']['large'], 'this results in an FDR of', np.round(100 * putFDRlarge, 3), 'percent of large putative diploid X chrom calls'
	print '\t\t\tComparing', putSmallRate, 'to', femaleRateDict['putative']['small'], 'this results in an FDR of', np.round(100 * putFDRsmall, 3), 'percent of small putative diploid X chrom calls'	

	print '\t\t', filteredCNVnum, 'CNVs remain after filtering for a rate of', filteredCNVrate, 'CNVs per cell'
	print '\t\t\tComparing to', femaleRateDict['filtered']['all'], 'results in an FDR of', np.round(100 * filtFDRall, 3), 'percent of all filtered diploid X chrom calls'
	print '\t\t\tComparing', filtLargeRate, 'to', femaleRateDict['filtered']['large'], 'results in an FDR of', np.round(100 * filtFDRlarge, 3), 'percent of large filtered diploid X chrom calls'
	print '\t\t\tComparing', filtSmallRate, 'to', femaleRateDict['filtered']['small'], 'results in an FDR of', np.round(100 * filtFDRsmall, 3), 'percent of small filtered diploid X chrom calls'	

	print '\t\tPercent change and fold difference of FDR from filtering'
	print '\t\t\tAll CNVs\t', 100 * percentDiffAll, '\t', foldDiffAll
	print '\t\t\tLarge CNVs\t', 100 * percentDiffLarge, '\t', foldDiffLarge
	print '\t\t\tSmallCNVs\t', 100 * percentDiffSmall, '\t', foldDiffSmall
	









def runOneFNR(female, male, distance, femalePloidy, malePloidy, repeats, name, numBins, xBins, absPosDict, cutoffDict):
	###misc imports and variables###
  ##########Probably need to update filename
  #refFile = folders.getOutDir(name, 'prep') + 'xChromRef.txt'

	#parameterOptions#
	CNoptions = [1, 3]
#	sizeOptions = [np.arange(3, 11), np.arange(11, 28), np.arange(28, 101), np.arange(101, 1315)]
	sizeOptions = [[3]]
	while sizeOptions[-1][-1] != 1314:
		prevSize = len(sizeOptions[-1])
		thisSize = prevSize * 2
		thisStart = sizeOptions[-1][-1]+1
		thisEnd = thisStart + thisSize
		theseOptions = np.arange(thisStart, min(thisEnd, 1315))
		sizeOptions.append(theseOptions)


	comboName = female + 'V' + male
	femaleData = loadCountData(name, female, numBins, xBins, femalePloidy)
	maleData = loadCountData(name, male, numBins, xBins, malePloidy)


	#NOTE this end index is one greater than the modified bins for ease of python indexing strategy#
		#starts and ends are in the xChrom bin location, starting from 0 on the chromosome#
	resultDtype = {'names': ('name', 'distance', 'rep', 
			'start', 'end', 'size', 'CN', 
			'segStart', 'segEnd', 'segSize', 'segCN', 'segIntDist', 'putativeID', 'filteredID'), 
		'formats': ('S20', 'float', 'int',
			'int', 'int', 'int', 'int',
			'int', 'int', 'int', 'float', 'float', 'bool', 'bool')}
	resultsArray = np.zeros(repeats, dtype = resultDtype)

  ##########Probably need to update filename
#	refoutfile = folders.getOutDir(name, 'XchromFNR') + comboName + '.reference.txt'
	if os.path.exists(refoutfile):
		REF = open(refoutfile, 'a+')
		allData = REF.read()
	
		allDataList = allData.split('\n')[1:]
		for i,j in enumerate(allDataList):
			if i == len(resultsArray):
				break
			repData = j.split('\t')
			resultsArray[i]['name'] = comboName
			resultsArray[i]['distance'] = distance
			resultsArray[i]['rep'] = int(repData[0])
			resultsArray[i]['start'] = int(repData[1])
			resultsArray[i]['end'] = int(repData[2])
			resultsArray[i]['size'] = int(repData[2]) - int(repData[1])
			resultsArray[i]['CN'] = int(repData[3])

		REF.flush()
	else:
		REF = open(refoutfile, 'w')
		REF.write('Repeat\tStart\tEnd\tCN')



	for i,j in enumerate(resultsArray):
		thisRep = i+1
		matlabName = comboName + 'R' + str(thisRep)
      ##########Probably need to update filename
#		segFile = folders.getOutDir(name, 'XchromFNR') + matlabName + '.segments.txt'

		###processing only if this simulation has not already been run (will need to delete all files if I change strategies again)###
		if j['rep'] != thisRep:
			testData = np.copy(femaleData)
			j['name'] = comboName
			j['distance'] = distance
			j['rep'] = thisRep


			###determine simulation iteration variables (start, stop, ((size)), copy number)###
			j['CN'] = random.choice(CNoptions)
			sizeGroup = random.choice(sizeOptions)
			j['size'] = random.choice(sizeGroup)
			j['start'] = random.randint(0, len(xBins) - j['size'])
			j['end'] = j['start'] + j['size']

			REF.write('\n')
			REF.write(   str(  '\t'.join( map(str, [thisRep, j['start'], j['end'], j['CN']]) )  )   )



			###merge bin data, re-log-correct and save it###
				#losses simulated by replacing a random region in the female X with a random set of contiguous male X bins
				#gains simulated by adding together the randomly selected contiguous regions

			if j['CN'] == 1:
				testData[j['start']:j['end']] = maleData[j['start']:j['end']]
			elif j['CN'] == 3:
				testData[j['start']:j['end']] = femaleData[j['start']:j['end']] + maleData[j['start']:j['end']]
			testData = np.log2(testData / femalePloidy)
	
  ##########Probably need to update filename
#			dataFile = folders.getOutDir(name, 'XchromFNR') + matlabName + '.lowessBinCounts.txt'
			np.savetxt(dataFile, testData)


		
			###write and run matlab script###
			importLines = "sampleName = genvarname('" + matlabName + "');\n"
			importLines += "sampleData = importdata('" + dataFile + "');\n"
			importLines += "ref = importdata('" + refFile + "');\n"

			refLines = "chroms = zeros(length(ref.textdata), 1);\n\n"
			refLines += "for i = 1:length(ref.textdata)\n    chrom = textscan(ref.textdata{i},'%s %s','delimiter','r');\n    chroms(i) = 20;\nend\n\n"
			refLines += "bins = ref.data;\nlogLowess = sampleData;\nclear ref sampleData\n\n"

			CBSlines = "cbsInput = [chroms, bins, logLowess];\ncbsOutput = cghcbs(cbsInput, 'ALPHA', " + str(0.01) + ", 'PERMUTATIONS', 10000, 'STOPPINGRULE', false, 'SHOWPLOT', false);\n\n"
			CBSlines += "numOfSegments = 0;\nfor i = 1:length(cbsOutput.SegmentData)\n    for j = 1:length(cbsOutput.SegmentData(1, i).Mean)\n        numOfSegments = numOfSegments + 1;\n    end\nend\n\n"
			CBSlines += "segments = zeros(length(numOfSegments), 3);\nk = 0;\nfor i = 1:length(cbsOutput.SegmentData)\n    for j = 1:length(cbsOutput.SegmentData(1, i).Mean)\n        k = k + 1;\n        segments(k, 1) = cbsOutput.SegmentData(1, i).Start(j);\n        segments(k, 2) = cbsOutput.SegmentData(1, i).End(j);\n        segments(k, 3) = cbsOutput.SegmentData(1, i).Mean(j);\n    end\nend\n\n"
			CBSlines += "segFile = fopen([sampleName '.segments.txt'], 'w');\nfor i = 1:numOfSegments\n    fprintf(segFile,'%d\t',segments(i,1));\n    fprintf(segFile,'%d\t',segments(i,2));\n    fprintf(segFile,'%.15f\\n',segments(i,3));\nend\n\n"

			lastLine = "quit\n"

			allLines = '\n'.join([importLines, refLines, CBSlines, lastLine])
	
	  ##########Probably need to update filename
#		scriptFile = folders.getOutDir(name, 'XchromFNR') + matlabName + '.m'
			SCRIPT = open(scriptFile, 'w')
			SCRIPT.write(allLines)
			SCRIPT.close()

  ##########Probably need to update filename
#			stderrfile = folders.getOutDir(name, 'XchromFNR') + matlabName + '.stderr'
			stdout = open(stderrfile, 'w')

			matlabCommand = "matlab -nodisplay -r " + matlabName
			cmd = shlex.split(matlabCommand)
			p = sub.Popen(cmd, stdout=stdout, stderr=sub.STDOUT)
			p.wait()

			os.remove('./' + matlabName + '.m')
			os.remove('./' + matlabName + '.stderr')



		###assess if each CNV was identified and return information###
			# (cells, simNumber, start, stop, size, copy number, segStart, segStop, segSize, segCN, segIntdist, putativeID, filteredID)#
		segData = loadSegData(segFile, femalePloidy)
		segData['start'] = [absPosDict[x] - xBins[0] for x in segData['start']]
		segData['end'] = [absPosDict[x] - xBins[0] for x in segData['end']]

		cnvSeg = segData[np.argmax([min(x['end'], j['end']) - max(x['start'], j['start']) for x in segData])]
		j['segStart'] = cnvSeg['start']
		j['segEnd'] = cnvSeg['end']
		j['segSize'] = cnvSeg['end'] - cnvSeg['start']
		j['segCN'] = cnvSeg['CN']
		j['segIntDist'] = abs(cnvSeg['CN'] - np.round(cnvSeg['CN']))

		if 0.5 <= float(j['segSize']) / float(j['size']) <= 2 and np.round(j['segCN']) == j['CN']:
			j['putativeID'] = True

			if cutoffDict[j['segSize']] >= j['segIntDist']:
				j['filteredID'] = True

	REF.close()

	return list(resultsArray)










def runFNR(name, numBins, repeats, overwrite, ploidyDict, xBins, absPosDict, cutoffDict, FNRcomparisons):
	print '\n\n\nAssessing FNR by spiking male X chromosome bin counts to the female X chromosome'
  ##########Probably need to update filename
#	os.chdir(folders.getOutDir(name, 'XchromFNR'))
	resultDtype = {'names': ('name', 'distance', 'rep', 
			'start', 'end', 'size', 'CN', 
			'segStart', 'segEnd', 'segSize', 'segCN', 'segIntDist', 'putativeID', 'filteredID'), 
		'formats': ('S20', 'float', 'int',
			'int', 'int', 'int', 'int',
			'int', 'int', 'int', 'float', 'float', 'bool', 'bool')}
  ##########Probably need to update filename
#	outfilename = folders.getOutDir(name, 'stats') + 'xChrom.FNRtestData.txt'
	thisHeader = 'Comparison\tDistance\tRepeat\tSpikeinStart\tSpikeinEnd\tSpikeinSize\tSpikeinCN\tSegStart\tSegEnd\tSegSize\tSegCN\tSegIntDist\tSegPutativeID\tSegFileredID'



	if overwrite or not os.path.exists(outfilename):
		argList = [(x[0], x[1], x[2], ploidyDict[x[0]], ploidyDict[x[1]], repeats, name, numBins, xBins, absPosDict, cutoffDict,) for x in FNRcomparisons]

    
    pool = mp.Pool(mp.cpu_count())
    processes = [ pool.apply_async(runOneFDR, args=x) for x in argList ]
    pool.close()

    results = len(argList) * [ [] ]
    for i,j in enumerate(processes):
      j.wait()
      results[i] = j.get()


		print '\nFinished segmenting male X spike-ins'
		resultsArray = np.array(sum(results, []), dtype=resultDtype)

		np.savetxt(outfilename, resultsArray, fmt='%s', delimiter='\t', comments='', header=thisHeader)


		CNVnum = float(len(resultsArray))
		putativeID = resultsArray[resultsArray['putativeID'] == True]
		filteredID = resultsArray[resultsArray['filteredID'] == True]
		putativeFNR = 1. - (float(len(putativeID)) / CNVnum)
		filteredFNR = 1. - (float(len(filteredID)) / CNVnum)

	else:
		resultsArray = np.genfromtxt(outfilename, skiprows=1, dtype=resultDtype)










	def getFilteredID(putative, bins, intdist):
		thisTest = False
		if putative and cutoffDict[bins] >= intdist:
			thisTest = True
		return thisTest

	resultsArray['filteredID'] = [ getFilteredID(x['putativeID'], x['segSize'], x['segIntDist']) for x in resultsArray]
	np.savetxt(outfilename, resultsArray, fmt='%s', delimiter='\t', comments='', header=thisHeader)






	failData = resultsArray[resultsArray['putativeID'] == False]
	IDdata = resultsArray[resultsArray['putativeID'] == True]
	filtData = resultsArray[resultsArray['filteredID'] == True]
	putFNRall = 100 * (1. - (float(len(IDdata)) / float(len(resultsArray))))
	filtFNRall = 100 * (1. - (float(len(filtData)) / float(len(resultsArray))))

	largeData = resultsArray[resultsArray['size'] > sizeBoundary]
	largePutID = largeData[largeData['putativeID'] == True]
	largeFiltID = largeData[largeData['filteredID'] == True]
	putFNRlarge = 100 * (1. - (float(len(largePutID)) / float(len(largeData))))
	filtFNRlarge = 100 * (1. - (float(len(largeFiltID)) / float(len(largeData))))

	smallData = resultsArray[resultsArray['size'] < sizeBoundary]
	smallPutID = smallData[smallData['putativeID'] == True]
	smallFiltID = smallData[smallData['filteredID'] == True]
	putFNRsmall = 100 * (1. - (float(len(smallPutID)) / float(len(smallData))))
	filtFNRsmall = 100 * (1. - (float(len(smallFiltID)) / float(len(smallData))))



	print '\tAssessing all', len(resultsArray), 'CNVs'
	print '\t\t', len(IDdata), 'putative CNVs were identified for a FNR of', putFNRall, '%'
	print '\t\t', len(filtData), 'filtered CNVs were identified for a FNR of', filtFNRall, '%'
	print '\t\t', np.round(100 * ((filtFNRall - putFNRall) / putFNRall), 3), '% increase and', np.round(filtFNRall / putFNRall, 3), 'fold change after filtering'
	

	print '\tAssessing', len(largeData), 'large CNVs'
	print '\t\t', len(largePutID), 'putative CNVs were identified for a FNR of', putFNRlarge, '%'
	print '\t\t', len(largeFiltID), 'filtered CNVs were identified for a FNR of', filtFNRlarge, '%'
	print '\t\t', np.round(100 * ((filtFNRlarge - putFNRlarge) / putFNRlarge), 3), '% increase and', np.round(filtFNRlarge / putFNRlarge, 3), 'fold change after filtering'

	print '\tAssessing', len(smallData), 'small CNVs'
	print '\t\t', len(smallPutID), 'putative CNVs were identified for a FNR of', putFNRsmall, '%'
	print '\t\t', len(smallFiltID), 'filtered CNVs were identified for a FNR of', filtFNRsmall, '%'
	print '\t\t', np.round(100 * ((filtFNRsmall - putFNRsmall) / putFNRsmall), 3), '% increase and', np.round(filtFNRsmall / putFNRsmall, 3), 'fold change after filtering'










def runAll(name, numBins, repeats, xStart, xEnd, xBins, absPosDict, overwrite):

	cutoffDict, ploidyDict, FDRcomparisons, FNRcomparisons, femaleRateDict = getRefInfo(name, xStart, xEnd, absPosDict)

	runFDR(name, overwrite, FDRcomparisons, femaleRateDict, ploidyDict, xBins, numBins, absPosDict, cutoffDict)

	runFNR(name, numBins, repeats, overwrite, ploidyDict, xBins, absPosDict, cutoffDict, FNRcomparisons)






sizeBoundary = 25.5#####


##########LEFTOVER FROM ORIGINAL CODE, will hopefully be replaced by fixing comments above
folders = cfg.Folders()






