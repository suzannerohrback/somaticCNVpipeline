#!/usr/bin/python
import gzip
import numpy as np
import sys, os, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg










def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")
	rd = dict()
	
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		
		if rd.has_key(id):
			#rd[id].append(arow)
			print "duplicate knowngene id = " + id
			print "arow =   " + str(arow)
			print "rd[id] = " + str(rd[id])
			
		else:
			rd[id] = arow	
			
	input.close()
	return(rd)










def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")
	ra = []
	
	for i in range(skipFirst):
		input.readline()
		
	for x in input:
		arow = x.rstrip().split("\t")		
		ra.append(arow)
		
	input.close()
	return(ra)










def countBins(samFile, countFile, statFile, sizeRef, binRef):
	if samFile[-2:] == 'gz':
		INFILE = gzip.open(samFile, "rb")
	else:
		INFILE = open(samFile, "r")
		
	OUTFILE = open(countFile, "w")
	STATFILE = open(statFile, "w")

	chrominfo = fileToDictionary(sizeRef, 0)
	bins = fileToArray(binRef, 0)

	

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	counter = 0
	totalReads = 0
	prevChrompos = ""
	for x in INFILE:
		arow = x.rstrip().split("\t")
		if arow[0][0] == '@':
			continue
		thisChrom = arow[2]
		thisChrompos = arow[3]
		if thisChrom.find("_") > -1:
			continue
		if thisChrom == "chrM":
			continue
		if thisChrom == "":
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue

		totalReads += 1

		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
		
		counter += 1
		
		indexUp = len(bins) - 1
		indexDown = 0
		indexMid = int((indexUp - indexDown) / 2.0)

		while True:
			if thisAbspos >= long(bins[indexMid][2]):
				indexDown = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexMid
			else:
				indexUp = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexDown

			if indexUp - indexDown < 2:
				break

		#####I ADDED IN THIS IF/ELSE TEST -- ORIGINALLY ONLY THE SECOND STATEMENT WAS HERE BUT THAT PREVENTS GETTING COUNTS IN THE LAST BIN#####
		if thisAbspos >= long(bins[indexUp][2]):
			binCounts[indexUp] += 1
		else:
			binCounts[indexDown] += 1

		prevChrompos = thisChrompos
		
		
		
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

		
		
	binCounts.sort()

	STATFILE.write('Reads\t')
	STATFILE.write(str(totalReads))
	STATFILE.write('\n')
	
	STATFILE.write('AverageCount\t')
	STATFILE.write(str(np.mean(binCounts)))
	STATFILE.write('\n')
		
	STATFILE.write('MedianCount\t')
	STATFILE.write(str(binCounts[len(bins)/2]))
	STATFILE.write('\n')

	
	
	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


	
	
	
	
	
	
	
		
def runOne(samFile, countDir, statsDir, species):
  
	#get environment prepared#
	countVars = cfg.Count()

	if samFile[-2:] == 'gz':
		sampleName = os.path.basename(samFile)[:-13]
	else:
		sampleName = os.path.basename(samFile)[:-11]
  
	statFile = statsDir + sampleName + '.bincount.stats.txt'
	countFile = countDir + sampleName + '.bincounts.txt'

	
	
	countBins(samFile, countFile, statFile, countVars.chromDict[species], countVars.binDict[species])
	
	
	
	printText = '\t\tFinished counting reads for ' + os.path.basename(samFile)
	print(printText)

	
	
	
	
