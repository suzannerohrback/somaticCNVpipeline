#!/usr/bin/python
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg










def calcReads(sample, statsDir):
	infile = statsDir + sample + '.bincount.stats.txt'
	
	with open(infile, 'r') as IN:
		data = IN.readline()
		readCount = data.rstrip().split('\t')[1]
		
	return int(readCount)
	
	
	
	
	
	
	
	
	
	
def calcMAPD(sample, lowessDir):
	infile = lowessDir + sample + '.lowess.txt'
	
	data = np.loadtxt(infile)
	MAPD = np.median(abs(data[1:] - data[:-1]))
	
	return MAPD, data
	
	
	
	
	
	
	
	
	
	
def calcCS(data):
	CS = 1 - ( 2 * np.median( abs( data - np.round(data) ) ) )
	return CS
	
	
	
	
	
	
	
	
	
	
def getPloidy(segData):
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

	#DO I WANT THIS FUNCTION? DOESN'T YET EXIST#
#	thisPlots.makeQCplot(ploidyTestValues, CSarray, CScutoff, peakCS, peakPloidy)

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










def runQCone(sample, species, statsDir, lowessDir, segmentDir, QCdir):
	#import config info#
	interpretVars = cfg.Interpret()
	binArray = common.importInfoFile(interpretVars.binDict[species], [0, 2, 4, 5], 'normref', skiprows=1)

	
	
	#determine read counts from the countstats file#
	readCount = calcReads(sample, statsDir)
	
	
	
	#determine MAPD from lowess bins#
	MAPD, lowessData = calcMAPD(sample, lowessDir)
	
	
	
	#determine the optimal ploidy value and CS#
	segData, segArray = common.importSegData(sample, segmentDir, binArray)
	CS, ploidy = getPloidy(segArray)
		
	
	
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
	
	
	
	
	
	
	
	
	
	
