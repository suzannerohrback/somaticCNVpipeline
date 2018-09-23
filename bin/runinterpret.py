#!usr/bin/python
import os
import numpy as np

import common
from interpret import qcfile, funcfile, analyzefiles









def runAll(args):

	print('\n\n\nYou have requested to analyze CNV call data')
	print('\tWARNING:')
	print('\t\tIF USING ANY REFERENCES OTHER THAN THOSE I PROVIDE I CANNOT GUARANTEE RESULT ACCURACY')
	print('\n')


	
	
	
	#Set up environment#
	args.AnalysisDirectory = common.fixDirName(args.AnalysisDirectory)
	
	
	
	folderDict = {'LowessBinCounts': args.lowess, 
		    'Segments': args.segments, 
		    'PipelineStats': args.countstats}
	
	for i in list(folderDict.keys()):
		if not folderDict[i]:
			folderDict[i] = args.AnalysisDirectory + i + '/'
		else:
			folderDict[i] = common.fixDirName(folderDict[i])
	
	
	
	QCdir = args.AnalysisDirectory + 'QC/'
	CNVdir = args.AnalysisDirectory + 'CNVlists/'
	summaryDir = args.AnalysisDirectory + 'SummaryFiles/'
	PloidyPlotDir = args.AnalysisDirectory + 'PloidyDeterminationPlots/'
	CNplotDir = args.AnalysisDirectory + 'CopyNumberProfilePlots/'
	ChromPlotDir = args.AnalysisDirectory + 'ChromosomeCopyNumberPlots/'
	
	for i in [args.AnalysisDirectory, QCdir, CNVdir, summaryDir, PloidyPlotDir, CNplotDir, ChromPlotDir]:#
		common.makeDir(i)
	
	
	
	#get list of samples to process 
		#will involve checking infofile (if present) and whether required input files exist
	sampleFiles = common.getSampleList(folderDict['Segments'], args.samples, 'segments')
	sampleNames = [x.split('/')[-1].split('.')[0] for x in sampleFiles]

#	info = common.importInfoFile(args.infofile, args.columns, 'interpret')
#	if args.infofile:
#		refArray = info
#	else:
#		thisDtype = info
#		refArray = np.array(
#			[ (x, 1, 'unk',) for x in sampleNames],
#			dtype=thisDtype)
		

	
	
	
	#QC assessment#
#	qcfile.runQCone(sampleNames[0], args.species, folderDict['PipelineStats'], folderDict['LowessBinCounts'], folderDict['Segments'], QCdir, PloidyPlotDir)
	argList = [(x, args.species, folderDict['PipelineStats'], folderDict['LowessBinCounts'], folderDict['Segments'], QCdir, PloidyPlotDir) for x in sampleNames]
	common.daemon(qcfile.runQCone, argList, 'assess sample quality')

	analysisSamples = []
	ploidyDict = {}
	genderDict = {}
	
	mergeQCfile = summaryDir + 'QCmetrics.txt'
	OUT = open(mergeQCfile, 'w')
	OUT.write('Name\tReads\tMAPD\tCS\tPloidy\tGender\tPASS\n')
	
	for i in sampleNames:
		IN = open(QCdir + i + '.qcTEMP.txt', 'r')
		data = IN.readline()	
		OUT.write(data)
		
		data = data.rstrip().split('\t')
		if data[-1] == 'True':
			analysisSamples.append(i)
			ploidyDict[i] = float(data[4])
			genderDict[i] = data[-2]
		
		IN.close()
		os.remove(QCdir + i + '.qcTEMP.txt')
		
	OUT.close()
	os.rmdir(QCdir)
	
	
	
	#FUnC: CNV filtering#
	if args.nofilter:
		print '\nFURTHER CODE IS ONLY DEVELOPED FOR WHEN FUnC IS IMPLEMENTED, EXITING NOW\n\n\n'
		raise SystemExit
		
#	funcfile.FUnCone(analysisSamples[0], args.species, folderDict['Segments'], CNVdir, 
#			 ploidyDict[analysisSamples[0]], genderDict[analysisSamples[0]])
	argList = [(x, args.species, folderDict['Segments'], CNVdir, ploidyDict[x], genderDict[x]) for x in analysisSamples]
	common.daemon(funcfile.FUnCone, argList, 'remove unreliable CNV calls')
	
	
	
	#CNV analysis#
#	summaryStats = analyzefiles.analyzeOne(analysisSamples[0], args.species, CNVdir, folderDict['LowessBinCounts'], CNplotDir, ChromPlotDir, ploidyDict[analysisSamples[0]], genderDict[analysisSamples[0]])
#	summaryStats = [summaryStats]
	argList = [(x, args.species, CNVdir, folderDict['LowessBinCounts'], CNplotDir, ChromPlotDir, ploidyDict[x], genderDict[x]) for x in analysisSamples]
	summaryStats = common.daemon(analyzefiles.analyzeOne, argList, 'create summary files')
	
	cellStatsFile = summaryDir + 'CellStats.txt'
	chromAmpFile = summaryDir + 'ChromosomeAmplifiedPercent.txt'
	chromDelFile = summaryDir + 'ChromosomeDeletedPercent.txt'
	
	#write summary statistics files#
	with open(cellStatsFile, 'w') as CELL, open(chromAmpFile, 'w') as AMP, open(chromDelFile, 'w') as DEL:
		CELL.write('Sample\tDeletionNumber\tAmplificationNumber\tTotalCNVnumber\tDeletedMB\tAmplifiedMB\tNetDNAalterdMB\n')
		chromHeader = 'Sample\t' + '\t'.join(summaryStats[0]['chroms']) + '\n'
		AMP.write(chromHeader)
		DEL.write(chromHeader)
		
		for i,j in enumerate(analysisSamples):
			CELL.write(str(j + '\t'))
			cellOut = [summaryStats[i]['cellStats']['delCount'],
				   summaryStats[i]['cellStats']['ampCount'],
				   summaryStats[i]['cellStats']['delCount'] + summaryStats[i]['cellStats']['ampCount'],
				   np.round(summaryStats[i]['cellStats']['delMB'], 3),
				   np.round(summaryStats[i]['cellStats']['ampMB'], 3),
				   np.round(summaryStats[i]['cellStats']['ampMB'] - summaryStats[i]['cellStats']['delMB'], 3)]
			cellOut = '\t'.join(map(str, cellOut)) + '\n'
			CELL.write(cellOut)
			
			AMP.write(str(j + '\t'))
			ampOut = [np.round(summaryStats[i]['chromAmp'][x], 3) for x in summaryStats[0]['chroms']]
			ampOut = '\t'.join(map(str, ampOut)) + '\n'
			AMP.write(ampOut)
					   
			DEL.write(str(j + '\t'))
			delOut = [np.round(summaryStats[i]['chromDel'][x], 3) for x in summaryStats[0]['chroms']]
			delOut = '\t'.join(map(str, delOut)) + '\n'
			DEL.write(delOut)
			
	
	
	
	print('\nCNV analysis complete\n\n\n')
	
	
	
	
	
	
