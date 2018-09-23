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
	summaryDir = args.AnalysisDirectory + 'CNVsummary/'
	CNplotDir = args.AnalysisDirectory + 'CopyNumberProfilePlots/'
	#ChromPlotDir = args.AnalysisDirectory + 'ChromosomeCopyNumberPlots/'
	#summaryPlotDir = args.AnalysisDirectory + 'CombinedSamplesPlots/'
	
	for i in [args.AnalysisDirectory, QCdir, CNVdir, summaryDir, CNplotDir]:#, ChromPlotDir]#, summaryPlotDir]:
		common.makeDir(i)
	
	
	
	#get list of samples to process 
		#will involve checking infofile (if present) and whether required input files exist
	sampleFiles = common.getSampleList(folderDict['Segments'], args.samples, 'segments')
	sampleNames = [x.split('/')[-1].split('.')[0] for x in sampleFiles]

	info = common.importInfoFile(args.infofile, args.columns, 'interpret')

	if args.infofile:
		refArray = info
	else:
		thisDtype = info
		refArray = np.array(
			[ (x, 1, 'unk',) for x in sampleNames],
			dtype=thisDtype)
		

	
	
	
	#QC assessment#
	argList = [(x, args.species, folderDict['PipelineStats'], folderDict['LowessBinCounts'], folderDict['Segments'], QCdir) for x in sampleNames]
	common.daemon(qcfile.runQCone, argList, 'assess sample quality')

	analysisSamples = []
	ploidyDict = {}
	genderDict = {}
	
	mergeQCfile = QCdir + 'ALL_SAMPLES.QC.txt'
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
	
	
	
	#FUnC: CNV filtering#
#	funcfile.FUnCone(analysisSamples[0], args.species, folderDict['Segments'], CNVdir, 
#			 ploidyDict[analysisSamples[0]], genderDict[analysisSamples[0]])
	argList = [(x, args.species, folderDict['Segments'], CNVdir, ploidyDict[x], genderDict[x]) for x in analysisSamples]
	common.daemon(funcfile.FUnCone, argList, ' remove unreliable CNV calls')
	
	
	
	
	#CNV analysis#
	analyzefiles.analyzeOne(analysisSamples[0], folderDict['Segments'], folderDict['LowessBinCounts'], CNplotDir, ploidyDict[analysisSamples[0]], genderDict[analysisSamples[0]])
#	argList = [(x, folderDict['Segments'], folderDict['LowessBinCounts'], CNplotDir, ploidyDict[x], genderDict[x]) for x in analysisSamples]
#	common.daemon(analyzefiles.analyzeOne, argList, ' create summary file(s)')
	
	raise SystemExit
	
	print('\nCNV analysis complete\n\n\n')
	
	
	
	
	
	
