#!/usr/bin/python
import os
import numpy as np

import common
from segment import normalizefile, segmentfile













def runAll(args):

	print('\n\n\nYou have requested to normalize and segment bincounts files')
	print('\tWARNING:')
	print('\t\tIF USING ANY REFERENCES OTHER THAN THOSE I PROVIDE I CANNOT GUARANTEE RESULT ACCURACY')
	print('\n')
 


	#Set up environment#
	#make folders, get list of sample names, etc
	args.CountDirectory = common.fixDirName(args.CountDirectory)
	
	lowessDir = os.path.dirname(args.CountDirectory[:-1]) + '/LowessBinCounts/'
	segmentDir = os.path.dirname(args.CountDirectory[:-1]) + '/Segments/'
	tempDir = os.path.dirname(args.CountDirectory[:-1]) + '/Temp/'

	if args.output:
		lowessDir = common.fixDirName(args.output) + 'LowessBinCounts/'
		segmentDir = common.fixDirName(args.output) + 'Segments/'

	common.makeDir(lowessDir)
	if not args.normalizeonly:
		common.makeDir(segmentDir)
		common.makeDir(tempDir)

	sampleFiles = common.getSampleList(args.CountDirectory, args.samples, 'bincounts')
		
		
		
		
		
	#Use multiprocessing daemon to run normalization for all samples#
	
	info = common.importInfoFile(args.infofile, args.columns, 'normalize')
	
	if args.infofile:
		refArray = info
	else:
		thisDtype = info
		refArray = np.array(
			[ (basename(x)[:-14], 'unk', 1,) for x in sampleFiles],
			dtype=thisDtype)
		
	sampleDict = {x: [y for y in sampleFiles if x in y][0] for x in refArray['name']}
		
		
	methodDict = {x: False for x in np.unique(refArray['method'])}	
	methodDict['NA'] = False
	sampleNormMethodDict = {x['name']: 'NA' for x in methodDict}
	
	if not args.gconly:
		for i in methodDict:
			refSlice = refArray[(refArray['method'] == i) & (refArray['cells'] == 1)]
			methodSamples = [sampleDict[x] for x in refSlice['name']]
			
			methodDict[i] = normalizefile.runMakeMethodRef(args.species, methodSamples, i, lowessDir)
		
			if methodDict[i] != False:
				for j in refSlice['name']:
					sampleNormMethodDict[j] = i
		
		
		
		
		
	#run multiprocessing for gc (+ method) correction
	argList = [(args.species, sampleDict[x], methodDict[sampleNormMethodDict[x]], lowessDir + x + '.lowess.txt') for x in sampleDict]
	common.daemon(normalizefile.runNormalizeOne, argList, 'normalize bincount files')

	
	
	print('\nNormalization complete\n\n\n')
	
	if args.normalizeonly:
		return 0
		
		
		
		
		
	#Run segmentation#
	#write matlab script
	#run matlab script
	#run for all samples in parallel with multiprocessing daemon
	
	
	
	
	
	printText = '\nTHIS FUNCTION IS NOT YET COMPLETE, SORRY\n'
	print(printText)
	raise SystemExit
	

	
	print('\nSegmentation complete\n\n\n')

	
	
	
