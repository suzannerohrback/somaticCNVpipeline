#!usr/bin/python
import os

import common
from interpret import qcfile, funcfile, analyzefiles









def runAll(args):

	print('\n\n\nYou have requested to analyze CNV call data')
	print('\tWARNING:')
	print('\t\tIF USING ANY REFERENCES OTHER THAN THOSE I PROVIDE I CANNOT GUARANTEE RESULT ACCURACY')
	print('\n')


	
	
	
	#Set up environment#
	args.AnalysisDirectory = common.fixDirName(args.AnalysisDirectory)
	
	
	
	folderDict = {'Lowess': args.lowess, 
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
	ChromPlotDir = args.AnalysisDirectory + 'ChromosomeCopyNumberPlots/'
	summaryPlotDir = args.AnalysisDirectory + 'CombinedSamplesPlots/'
	
	for i in [args.AnalysisDirectory, QCdir, CNVdir, summaryDir, CNplotDir, ChromPlotDir, summaryPlotDir]:
		common.makeDir(i)
	
	
	
	#get list of samples to process 
		#will involve checking infofile (if present) and whether required input files exist
	sampleFiles = common.getSampleList(folderDict['Segments'], args.samples, 'segments')

	info = common.importInfoFile(args.infofile, args.columns, 'interpret')

	if args.infofile:
		refArray = info
	else:
		thisDtype = info
		refArray = np.array(
			[ (basename(x)[:-13], 1, 'unk',) for x in sampleFiles],
			dtype=thisDtype)
		

	
	
	errorText = 'SORRY, THIS FUNCTION IS STILL BEING WRITTEN, TRY AGAIN LATER\n\n\n'
	print(errorText)
	raise SystemExit
	
	
	
	#QC assessment#
	
	#will be in bin/interpret/qcfile.py
	#read number, MAPD calculation, determine gender, ploidy, CS calculation
	#ploidy determination plots made here
	#save summary file and get a list of high quality samples
	
	
	
	
	
	
	#CNV filtering#
	
	#will be in bin/interpret/funcfile.py
	#for high quality samples, read through the segments.txt files
	#caclulate size and int dist of each call
	#merge adjacent segments with the same CN on the same chrom
	#if not args.nofiler run full FUnC, otherwise remove CNVs if < 3 bins
	#for each file save a CNV location summary file
	
	
	
	
	#CNV analysis#
	
	#will be in bin/interpret/analyzefiles.py
	#make CellStats, CNVstats, ChromosomeStats summary files
	#make Copy number profile plot and Chromosome copy number plot for each sample
	#make Heatmap of CNV locations and plot of genomic locus CNV frequency for all samples (group separation)
	#any sort of statistical testing????? (aneuploidy, CNVs/cell, ampVdel, DNA altered, CNV size...?)
	
	
	
	print('\nCNV analysis complete\n\n\n')

	
	
	############
#	parser.add_argument('species', choices=['hg38', 'mm10'], 
#		help = 'The genome build of the species being assessed')
#	
#	#optional arguments#
#	parser.add_arugment('-f', '--nofilter', action='store_true'
#		help = 'Set this flag if you do not want to perform FUnC filtering of low-quality CNV calls')
#	parser.add_argument('-i', '--infofile', metavar='/path/to/sample.info.txt', default=False,
#		help='Path to a .txt file containing information about the samples to be processed (unique name, number of cells, group)\n\tIf not all are identical. This file should not have a header row')
#	parser.add_argument('-c', '--columns', metavar='X X X', default=[0, 1, 2], type=int, nargs=3,
#		help='The zero-indexed locations of the columns to import from the infofile in the order: name, cell number, group (if not the first 3 columns)')
#	parser.add_argument('-s', '--samples', metavar='/path/to/sample_list.txt', default=False,
#		help='Path to a file containing a list of sample names to be processed\n\tno path or file extension needed')

	
	
	
	
	
	
	
